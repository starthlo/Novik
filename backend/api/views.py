import csv
import logging
import os
import random
import re
import string
import time

import nltk
import openai
from Bio import Entrez
from django.contrib.gis.geoip2 import GeoIP2
from django.db.models import F
from django.http import HttpResponse
from django.utils import timezone
from dotenv import load_dotenv
from drf_spectacular.utils import OpenApiParameter, OpenApiResponse, extend_schema
from rest_framework import status, viewsets
from rest_framework.decorators import action, api_view, permission_classes
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from .models import Banner, BannerStat, CustomUser, PatientContext
from .schemas.openai_schemas import (
    ErrorResponseSerializer,
    OpenAIPDFRequestSerializer,
    OpenAIResponseSerializer,
)
from .schemas.user_schemas import (
    ClearSessionsResponseSerializer,
    ToggleUserStatusRequestSerializer,
    ToggleUserStatusResponseSerializer,
    UserDeleteRequestSerializer,
    UserDeleteResponseSerializer,
)
from .serializers import (
    BannerSerializer,
    BannerStatSerializer,
)

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")

load_dotenv()

logger = logging.getLogger(__name__)

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
ASSISTANT_ID = os.environ.get("OPENAI_ASSISTANT_ID")

Entrez.email = os.environ.get("ENTREZ_EMAIL")

client = openai.OpenAI(api_key=OPENAI_API_KEY)


@extend_schema(
    tags=["Dashboard"],
    description="Process a text message or PDF file with message through OpenAI",
    request=OpenAIPDFRequestSerializer,
    responses={
        200: OpenApiResponse(
            response=OpenAIResponseSerializer, description="Successful response"
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
        500: OpenApiResponse(
            response=ErrorResponseSerializer, description="Server error"
        ),
    },
)
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def patient_assistant_view(request):
    """
    Process a text message or PDF file with message through OpenAI.

    This endpoint can handle:
    1. Text-only messages - just provide the 'message' parameter
    2. PDF files with messages - provide both 'pdf' file and 'message' parameter

    The function extracts text from PDF (if provided), combines with the message,
    sends it to OpenAI, and returns an enhanced response with additional context
    and references if available.
    """
    try:
        content = request.data.get("message")
        pdf_file = request.FILES.get("pdf")

        if not content and not pdf_file:
            return Response(
                {"error": "Either message content or PDF file is required"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # If PDF is provided, process it and combine with message content
        if pdf_file:
            import PyPDF2

            reader = PyPDF2.PdfReader(pdf_file)

            text = ""
            for page in reader.pages:
                text += page.extract_text() or ""

            if not text.strip():
                return Response(
                    {"error": "Could not extract text from PDF"},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if content:
                content_to_process = text + "\r\n" + content
            else:
                content_to_process = text
        else:
            if not content:
                return Response(
                    {"error": "Message content is required when no PDF is provided"},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            content_to_process = content

        enhanced_ai_response = openai_response_shared_code(request, content_to_process)

        if isinstance(enhanced_ai_response, Response):
            return enhanced_ai_response

        return Response({"message": enhanced_ai_response}, status=status.HTTP_200_OK)

    except Exception as e:
        logger.error(f"patient_assistant_view error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


def openai_response_shared_code(request, content):
    try:
        sessionId = request.data.get("session_id")
        createdBy = request.data.get("createdBy", 0)

        PatientContext.objects.create(
            session_id=sessionId,
            content=content,
            created_by=createdBy,
            created_at=timezone.now().date(),
        )

        ## retrieve context i.e. previous & current chat
        chat_context = ""
        pcs = list(PatientContext.objects.filter(session_id=sessionId).order_by("id"))
        for i, pc in enumerate(pcs):
            suffix = " [LANGUAGE]" if i == len(pcs) - 1 else ""
            chat_context += "\r\n" + pc.content.strip() + suffix

        print(f"chat_context:\n {chat_context}")

        thread = client.beta.threads.create()
        client.beta.threads.messages.create(
            thread_id=thread.id, role="user", content=chat_context
        )
        run = client.beta.threads.runs.create(
            thread_id=thread.id, assistant_id=ASSISTANT_ID
        )

        while True:
            run_status = client.beta.threads.runs.retrieve(
                thread_id=thread.id, run_id=run.id
            )
            if run_status.status == "completed":
                break
            time.sleep(1)

        messages = client.beta.threads.messages.list(thread_id=thread.id)
        m_d_c_t = messages.data[0].content[0].text
        ai_response = m_d_c_t.value
        ai_annotations = m_d_c_t.annotations

        print("RAW ai_response:\r\n ******* \r\n %s \r\n *******" % ai_response)

        anchor_random_string = "".join(
            random.choices(string.ascii_letters + string.digits, k=5)
        )

        # Iterate over the annotations and add footnotes
        citations = []
        for index, annotation in enumerate(ai_annotations):
            anchor_id = "%s%s" % (index + 1, anchor_random_string)
            ## Replace the text with a footnote
            ai_response = ai_response.replace(
                annotation.text,
                f" [^{index + 1}^](#{anchor_id})",
            )

            # Gather citations based on annotation attributes
            if file_citation := getattr(annotation, "file_citation", None):
                cited_file = client.files.retrieve(file_citation.file_id)
                citations.append(f"[{index + 1}] {cited_file.filename}")
            elif file_path := getattr(annotation, "file_path", None):
                cited_file = client.files.retrieve(file_path.file_id)
                citations.append(
                    f"[{index + 1}] Click <here> to download {cited_file.filename}"
                )
                # Note: File download functionality not implemented above for brevity

        if citations:
            ref_index = ai_response.find("### References")
            if ref_index != -1:
                ai_response = ai_response[:ref_index]
            else:
                ref_index = ai_response.find("**References")
                if ref_index != -1:
                    ai_response = ai_response[:ref_index]

            ai_response += (
                "\r\n"
                + "### References:  "
                + "\r\n"
                + "  \r\n".join(citations)
                + "\r\n"
            )

        ai_response = re.sub(
            r"^.*https?:\/\/(?:www\.|go\.)?drugbank\.(?:com|ca)\/drugs\/DB\d+.*\n?",
            "",
            ai_response,
            flags=re.MULTILINE,
        )

        updated_response = enrich_response(ai_response, anchor_random_string)

        return updated_response

    except Exception as e:
        logger.error(f"OpenAI API error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


def enrich_response(text, anchor_random_string):
    medications = {
        "Amoxicillin": "DB01060",
        "Clindamycin": "DB01190",
        "Metronidazole": "DB00916",
        "Azithromycin": "DB00207",
        "Erythromycin": "DB00199",
        "Cephalexin": "DB00567",
        "Cefuroxime": "DB01112",
        "Doxycycline": "DB00254",
        "Tetracycline": "DB00759",
        "Ciprofloxacin": "DB00537",
        "Levofloxacin": "DB01137",
        "Spiramycin": "DB06145",
        "Vancomycin": "DB00512",
        "Rifampicin": "DB01045",
        "Paracetamol": "DB00316",
        "Ibuprofen": "DB01050",
        "Naproxen": "DB00788",
        "Ketorolac": "DB00465",
        "Diclofenac": "DB00586",
        "Metamizole": "DB04817",
        "Celecoxib": "DB00482",
        "Etoricoxib": "DB01628",
        "Tramadol": "DB00193",
        "Codeine": "DB00318",
        "Morphine": "DB00295",
        "Acetylsalicylic acid": "DB00945",
        "Meloxicam": "DB00814",
        "Piroxicam": "DB00554",
        "Indomethacin": "DB00328",
        "Flurbiprofen": "DB00712",
        "Ketoprofen": "DB01009",
        "Prednisone": "DB00635",
        "Dexamethasone": "DB01234",
        "Hydrocortisone": "DB00741",
        "Lornoxicam": "DB06725",
        "Aceclofenac": "DB06736",
        ## post 2025-05-15
        "Articaine": "DB09009",
        "Epinephrine": "DB00668",
        "Lidocaine": "DB00281",
        "Mepivacaine": "DB00961",
        "Bupivacaine": "DB00297",
    }

    for med, dbid in medications.items():
        text = re.sub(
            rf"(?<!\[)({med})(?![\w\s]*\])",
            rf"[{med}](https://go.drugbank.com/drugs/{dbid})",
            text,
        )

    ref_block_match = re.search(r"### References:?(.*\n){1,9}", text)
    raw_titles = []
    if ref_block_match:
        ref_block = ref_block_match.group()
        raw_titles = re.findall(r"\[\d+\]\s*(.*)\n", ref_block)
    else:
        print("no matches for References")

    ## strip elements as we might catch new lines etc
    strip_titles = list(filter(str.strip, raw_titles))

    ## Ensure only non empty elements
    non_empty_titles = list(filter(None, strip_titles))

    ## Ensure only unique elements
    titles = list(set(non_empty_titles))

    print("-----------\r\n titles: %s \r\n-----------" % titles)

    result_limit = 1
    match_found = False
    all_refs = []
    for t in titles:
        t_strip = t.strip()

        ## Prevent PubMed from getting empty terms
        if len(t_strip) == 0:
            continue

        print("--R--> getting pubmed articles for term: %s" % t_strip)

        articles = get_pubmed_articles(t_strip, count=result_limit)
        if articles:
            match_found = True
            all_refs.extend(articles)
        else:
            all_refs.append({"title": t_strip})

    ## Enrich the bulleted lines just below these headings
    refmes = ["Preoperative Precautions", "Postoperative Guidelines"]

    ## Resume the anchor order from "References" above
    anchor_index = len(all_refs)

    t_w_r_s = []

    for refme in refmes:
        ## get the bulleted list below the heading
        pattern = re.compile(rf"{re.escape(refme)}.*\n?([-\s+|\d\.\s+].*\n){{1,4}}")
        ref_block_match = pattern.search(text)

        pps = []
        if ref_block_match:
            ref_block = ref_block_match.group()

            pps = re.findall(r"[-|\d\.]\s+(.*)\n", ref_block)

            print("~~~~~~\r\n%s pps:\r\n%s\r\n~~~~~~" % (refme, pps))

            result_limit = 3

            for eye, pp_content in enumerate(pps):
                t_strip = pp_content.strip()

                ## Prevent PubMed from getting empty terms
                if len(t_strip) == 0:
                    continue

                ## 20250521 Our regex ain't perfect :/
                if t_strip == "Anesthetics":
                    continue

                ## 20250521 Our regex ain't perfect :/
                # if t_strip.startswith("**"):
                # continue

                ## Sometimes the AI return these **BLABLABLA** text in the beginning of the sentence
                t_strip = t_strip.replace("**", "")

                print("~~~~> getting pubmed articles for term: %s" % t_strip)

                # articles = get_pubmed_articles(t_strip, count=result_limit)
                articles = []

                # publication_found = True if random.randint(0, 9) < 5 else False

                print("~o~~> got %s articles for term: %s" % (len(articles), t_strip))

                ## say we found one ...
                if articles:
                    match_found = True
                    all_refs.extend(articles)

                    for article in articles:
                        anchor_index = anchor_index + 1
                        t_w_r_s.append(
                            {
                                "refme": refme,
                                "eye": eye,
                                "content": pp_content,
                                "anchor_index": anchor_index,
                            }
                        )

        else:
            print("no list matches %s" % refme)

    # extra_refs = []

    ## Append superscripts
    for eye, t_w_r in enumerate(t_w_r_s):
        ## Next one should be the publication title, not the text it was derived from
        a_c = t_w_r["content"]
        a_i = t_w_r["anchor_index"]
        a_f = "%s%s" % (a_i, anchor_random_string)
        a_c_a_f = "%s %s" % (a_c, f"[^{a_i}^](#{a_f})")
        ## append superscripts. see line 184-ish in django/novik/api/views.py
        text = text.replace(a_c, a_c_a_f)
        ## append references. see line 365-ish in django/novik/api/views.py
        # extra_refs.append("###### %s. %s {{#{%s}}}" % (a_i, a_c, a_f))

    # match_found = True if len(t_w_r_s) > 0 else match_found

    print("-----------\r\n all_refs: %s \r\n-----------" % all_refs)

    if match_found:
        ref_block = "\n\n".join(
            [
                f"###### {i + 1}. [{ref['title']}]({ref['url']}) {{#{'%s%s' % (i + 1, anchor_random_string)}}}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
                for i, ref in enumerate(all_refs)
            ]
        )
        ## Replace (and remove) references from raw AI response
        text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", "", text)

        ## Add new PubMed references on top of raw AI response
        text = "%s  \r\n\r\n%s\r\n%s  \r\n" % (
            text,
            "-------\n### PubMed References",
            ref_block,
        )
    else:
        text = "%s  \r\n\r\n%s\r\n" % (
            text,
            "**No recent PubMed references were found that specifically relate to these factors.**",
        )

    return text


def get_nouns(text):
    tokens = nltk.word_tokenize(text)
    tagged = nltk.pos_tag(tokens)
    nouns = [word for word, pos in tagged if pos.startswith("NN")]
    return nouns


def get_pubmed_articles(query, count=2):
    print("query = %s" % query)

    nouns = get_nouns(query)
    nouns_joined = " OR ".join(nouns)

    print("nouns_joined = %s" % nouns_joined)

    term = (
        "%s AND (dental[mh] OR dentistry[mh] OR endodontics[mh] OR periodontics[mh] OR maxillofacial[mh])"
        % nouns_joined
    )

    print("term = %s" % term)

    ymax = timezone.now().date().year
    ymin = ymax - 5

    # print("ymin = %s; ymax = %s" % (ymin, ymax))

    handle = Entrez.esearch(
        db="pubmed",
        term=term,
        datetype="pdat",
        mindate=ymin,
        maxdate=ymax,
        retmax=count,
        sort="relevance",
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    articles = []

    for article in records["PubmedArticle"]:
        art = article["MedlineCitation"]["Article"]
        title = art["ArticleTitle"]
        pmid = article["MedlineCitation"]["PMID"]

        # Author list (up to 6 then "et al.")
        authors = art.get("AuthorList", [])
        author_names = []
        for a in authors[:6]:
            last = a.get("LastName", "")
            initials = a.get("Initials", "")
            if last:
                author_names.append(f"{last} {initials}")
        if len(authors) > 6:
            author_names.append("et al.")
        author_str = ", ".join(author_names)

        # Journal info
        journal = art["Journal"]["Title"]
        year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "n.d.")
        volume = art["Journal"]["JournalIssue"].get("Volume", "")
        issue = art["Journal"]["JournalIssue"].get("Issue", "")
        pages = art.get("Pagination", {}).get("MedlinePgn", "")

        # DOI (if available)
        dois = [
            str(el)
            for el in art.get("ELocationID", [])
            if getattr(el, "attributes", {}).get("EIdType") == "doi"
        ]
        doi = f"doi:{dois[0]}" if dois else ""
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        articles.append(
            {
                "title": title,
                "authors": author_str,
                "journal": journal,
                "year": year,
                "volume": volume,
                "issue": issue,
                "pages": pages,
                "doi": doi,
                "url": url,
            }
        )

    return articles


@extend_schema(tags=["Banners"])
class BannerViewSet(viewsets.ModelViewSet):
    """
    ViewSet for managing banner advertisements.

    Provides CRUD operations for banners and additional actions to track views and clicks.
    """

    queryset = Banner.objects.all()
    serializer_class = BannerSerializer
    parser_classes = (MultiPartParser, FormParser)

    def get_serializer_context(self):
        # ensure .context['request'] is available to BannerSerializer
        context = super().get_serializer_context()
        context["request"] = self.request
        return context

    @extend_schema(
        description="Record a view for this banner",
        responses={
            200: OpenApiResponse(
                {"type": "object", "properties": {"status": {"type": "string"}}},
                description="View recorded successfully",
            ),
        },
    )
    @action(detail=True, methods=["post"])
    def view(self, request, pk=None):
        """
        Record a view for this banner.

        Increments the view count for this banner and records the country of origin
        based on the IP address.
        """
        banner = self.get_object()
        ip = request.META.get("REMOTE_ADDR", "")
        try:
            country = GeoIP2().country(ip).get("country_name", "Unknown")
        except Exception:
            country = "Unknown"
        stat, _ = BannerStat.objects.get_or_create(
            banner=banner, date=timezone.now().date(), country=country
        )
        stat.views = F("views") + 1
        stat.save()
        return Response({"status": "view recorded"})

    @extend_schema(
        description="Record a click for this banner",
        responses={
            200: OpenApiResponse(
                {"type": "object", "properties": {"status": {"type": "string"}}},
                description="Click recorded successfully",
            ),
        },
    )
    @action(detail=True, methods=["post"])
    def click(self, request, pk=None):
        """
        Record a click for this banner.

        Increments the click count for this banner and records the country of origin
        based on the IP address.
        """
        banner = self.get_object()
        ip = request.META.get("REMOTE_ADDR", "")
        try:
            country = GeoIP2().country(ip).get("country_name", "Unknown")
        except Exception:
            country = "Unknown"
        stat, _ = BannerStat.objects.get_or_create(
            banner=banner, date=timezone.now().date(), country=country
        )
        stat.clicks = F("clicks") + 1
        stat.save()
        return Response({"status": "click recorded"})


@extend_schema(
    tags=["Banners"],
    parameters=[
        OpenApiParameter(
            name="banner",
            location=OpenApiParameter.QUERY,
            description="Filter by banner ID",
            required=False,
            type=int,
        ),
    ],
)
class BannerStatViewSet(viewsets.ReadOnlyModelViewSet):
    """
    ViewSet for retrieving banner statistics.

    Provides read-only access to banner statistics, with optional filtering by banner ID.
    """

    queryset = BannerStat.objects.all()
    serializer_class = BannerStatSerializer

    def get_queryset(self):
        """
        Filter queryset by banner ID if provided in query parameters.
        """
        qs = super().get_queryset()
        banner_id = self.request.query_params.get("banner")
        if banner_id:
            qs = qs.filter(banner__id=banner_id)
        return qs


@extend_schema(
    tags=["Users"],
    description="Export all users to a CSV file",
    responses={
        200: OpenApiResponse(description="CSV file with user data"),
    },
)
@api_view(["GET"])
def export_users_csv(request):
    """
    Export all users to a CSV file.

    Returns a CSV file containing user data with fields: id, username, email,
    date_joined, dob, phone, occupation, country, state, city, is_staff,
    is_superuser.
    """
    qs = CustomUser.objects.all()
    fields = [
        "id",
        "username",
        "email",
        "date_joined",
        "dob",
        "phone",
        "occupation",
        "country",
        "state",
        "city",
        "is_staff",
        "is_superuser",
    ]
    resp = HttpResponse(content_type="text/csv")
    resp["Content-Disposition"] = 'attachment; filename="users.csv"'
    writer = csv.writer(resp)
    writer.writerow(fields)
    for u in qs:
        writer.writerow([getattr(u, f) for f in fields])
    return resp


@extend_schema(
    tags=["Users"],
    description="Toggle a user's active status",
    request=ToggleUserStatusRequestSerializer,
    responses={
        200: OpenApiResponse(
            response=ToggleUserStatusResponseSerializer,
            description="Status toggled successfully",
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
    },
)
@api_view(["POST"])
def toggle_user_active_status(request):
    """
    Toggle a user's active status.

    Takes a user ID and toggles the is_active field between 0 and 1.
    """
    try:
        eyed = request.data.get("id")

        uzzer = CustomUser.objects.get(pk=eyed)
        uzzer.is_active = 1 if uzzer.is_active == 0 else 0
        uzzer.save()

        return Response(
            {"status": "user status set to %s successfully" % uzzer.is_active}
        )
    except Exception:
        return Response({"status": "unable to toggle user active status"})


@extend_schema(
    tags=["Users"],
    description="Delete a user",
    request=UserDeleteRequestSerializer,
    responses={
        200: OpenApiResponse(
            response=UserDeleteResponseSerializer,
            description="User deleted successfully",
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
    },
)
@api_view(["DELETE"])
def user_delete(request):
    """
    Delete a user.

    Takes a user ID and permanently deletes the user from the database.
    """
    try:
        eyed = request.data.get("id")

        uzzer = CustomUser.objects.get(pk=eyed)
        uzzer.delete()

        return Response({"status": "user deleted successfully"})
    except Exception:
        return Response({"status": "unable to deleted user"})


@extend_schema(
    tags=["Sessions"],
    description="Clear all patient context sessions",
    responses={
        200: OpenApiResponse(
            response=ClearSessionsResponseSerializer,
            description="Sessions cleared successfully",
        ),
    },
)
@api_view(["GET"])
def clear_sessions(request):
    """
    Clear all patient context sessions.

    Removes all collected prompts from the database and returns the count
    of deleted rows.
    """
    ## Remove collected prompts
    pcs = PatientContext.objects.all()
    kount = pcs.count()
    pcs.delete()

    return Response({"status": "ok", "message": "cleared %d rows" % (kount)})
