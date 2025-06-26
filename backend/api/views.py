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
from django.contrib.auth import login
from django.contrib.gis.geoip2 import GeoIP2
from django.db.models import F
from django.http import HttpResponse
from django.utils import timezone
from dotenv import load_dotenv
from rest_framework import status, viewsets
from rest_framework.decorators import action, api_view
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.response import Response

from .models import Banner, BannerStat, CustomUser, PatientContext
from .serializers import (
    BannerSerializer,
    BannerStatSerializer,
    GoogleLoginSerializer,
    LoginSerializer,
    UserListSerializer,
    UserSerializer,
)

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")

load_dotenv()

logger = logging.getLogger(__name__)

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
ASSISTANT_ID = os.environ.get("OPENAI_ASSISTANT_ID")

Entrez.email = "stevedev0323@gmail.com"
client = openai.OpenAI(api_key=OPENAI_API_KEY)


@api_view(["POST"])
def register_view(request):
    serializer = UserSerializer(data=request.data)
    if serializer.is_valid():
        serializer.save()
        return Response(
            {"message": "Registration successful"}, status=status.HTTP_201_CREATED
        )
    logger.error(f"Registration failed: {serializer.errors}")
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["POST"])
def login_view(request):
    serializer = LoginSerializer(data=request.data)
    if serializer.is_valid():
        user = serializer.validated_data["user"]
        login(request, user)
        return Response(
            {
                "message": "Login successful",
                "user": {
                    "id": user.id,
                    "username": user.username,
                    "email": user.email,
                    "is_superuser": user.is_superuser,
                    "is_staff": user.is_staff,
                },
            },
            status=status.HTTP_200_OK,
        )
    logger.error(f"Login failed: {serializer.errors}")
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["POST"])
def google_login_view(request):
    serializer = GoogleLoginSerializer(data=request.data)
    if serializer.is_valid():
        user = serializer.validated_data["user"]
        login(request, user)
        return Response(
            {
                "message": "Google login successful",
                "user": {
                    "id": user.id,
                    "username": user.username,
                    "email": user.email,
                    "is_superuser": user.is_superuser,
                    "is_staff": user.is_staff,
                },
            },
            status=status.HTTP_200_OK,
        )
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["POST"])
def openai_response_view(request):
    try:
        content = request.data.get("message")
        if not content:
            return Response(
                {"error": "Message content is required"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        enhanced_ai_response = openai_response_shared_code(request, content)

        return Response({"message": enhanced_ai_response}, status=status.HTTP_200_OK)

    except Exception as e:
        logger.error(f"openai_response_view error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


@api_view(["POST"])
def openai_pdf_response_view(request):
    try:
        pdf_file = request.FILES.get("pdf")
        content = request.data.get("message")
        if not pdf_file:
            return Response(
                {"error": "PDF file is required"}, status=status.HTTP_400_BAD_REQUEST
            )

        # Read and decode the PDF content
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

        pdf_text_and_content = text + "\r\n" + content

        enhanced_ai_response = openai_response_shared_code(
            request, pdf_text_and_content
        )

        return Response({"message": enhanced_ai_response}, status=status.HTTP_200_OK)

    except Exception as e:
        logger.error(f"openai_pdf_response_view error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


def openai_response_shared_code(request, content):
    try:
        sessionId = request.data.get("sessionId")
        createdBy = request.data.get("createdBy", 0)

        chat = PatientContext.objects.create(
            session_id=sessionId,
            content=content,
            # created_by=request.user.id,
            created_by=createdBy,
            created_at=timezone.now().date(),
        )

        ## retrieve context i.e. previous & current chat
        chat_context = ""
        for pc in PatientContext.objects.filter(session_id=sessionId).order_by("id"):
            chat_context = chat_context + "\r\n" + pc.content

        print("chat_context:\r\n ~~~~~~~ \r\n %s \r\n ~~~~~~~" % chat_context)

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
                # f' [{index}])')
                # f' [{index+1}](#{index+1})')
                # f' <sup>[{index+1}](#{index+1})</sup>')
                # f' ^[{index+1}](#{index+1})^')
                # f' [^{index+1}^](#{index+1})')
                f" [^{index + 1}^](#{anchor_id})",
            )

            # Gather citations based on annotation attributes
            if file_citation := getattr(annotation, "file_citation", None):
                cited_file = client.files.retrieve(file_citation.file_id)
                ## "quote" has been removed from OpenAI's API
                ## See https://learn.microsoft.com/en-us/answers/questions/2123533/attributeerror-filecitation-object-has-no-attribut
                # citations.append(f'[{index}] {file_citation.quote} from {cited_file.filename}')
                citations.append(f"[{index + 1}] {cited_file.filename}")
            elif file_path := getattr(annotation, "file_path", None):
                cited_file = client.files.retrieve(file_path.file_id)
                citations.append(
                    f"[{index + 1}] Click <here> to download {cited_file.filename}"
                )
                # Note: File download functionality not implemented above for brevity

        # Add footnotes to the end of the message before displaying to user
        # ai_response += '\n' + '\n'.join(citations)
        if citations:
            # Remove existing references
            ref_index = ai_response.find("### References")
            if ref_index != -1:
                ai_response = ai_response[:ref_index]
            else:
                ref_index = ai_response.find("**References")
                if ref_index != -1:
                    ai_response = ai_response[:ref_index]

            # print("ai_response (pre custom refs):\r\n ******* \r\n %s \r\n *******" % ai_response)

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

        print(
            "ai_response (before enrich_response):\r\n ******* \r\n %s \r\n *******"
            % ai_response
        )

        updated_response = enrich_response(ai_response, anchor_random_string)

        print(
            "ai_response (after enrich_response):\r\n ******* \r\n %s \r\n *******"
            % updated_response
        )

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

    # ref_block_match = re.search(r"### References:?\s*((?:\d+\. .+\n?){1,5})", text)
    # titles = []
    # if ref_block_match:
    # ref_block = ref_block_match.group(1)
    # titles = re.findall(r'\d+\.\s*"(.+?)"', ref_block)

    ref_block_match = re.search(r"### References:?(.*\n){1,9}", text)
    raw_titles = []
    if ref_block_match:
        # print(ref_block_match.group())
        # print('---')
        # print(ref_block_match.group(1))
        # print('---')
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

    # print("-----------\r\n all_refs: %s \r\n-----------" % all_refs)

    # if match_found:
    # ref_block = "\n\n".join([
    ## f"{i+1}. [{ref['title']}]({ref['url']})\n{ref['authors']}.\n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
    ## f"###### {i+1}. [{ref['title']}]({ref['url']}) {{ #{i+1} }}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
    # f"###### {i+1}. [{ref['title']}]({ref['url']}) {{#{'%s%s' % (i+1, anchor_random_string)}}}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
    # for i, ref in enumerate(all_refs)
    # ])
    ## Replace (and remove) references from raw AI response
    ## text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", f"### References:\n{ref_block}\n", text)
    # text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", f"", text)

    ## Add new PubMed references on top of raw AI response
    # text = "%s  \r\n\r\n%s\r\n%s  \r\n" % (text, "### PubMed References:", ref_block )
    # else:
    # text = "%s  \r\n\r\n%s\r\n" % (text, "**No recent PubMed references were found that specifically relate to these factors.**")

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## Enrich the bulleted lines just below these headings
    refmes = ["Preoperative Precautions", "Postoperative Guidelines"]

    ## Resume the anchor order from "References" above
    anchor_index = len(all_refs)

    t_w_r_s = []
    publication_found = False

    for refme in refmes:
        ## get the bulleted list below the heading
        pattern = re.compile(rf"{re.escape(refme)}.*\n?([-\s+|\d\.\s+].*\n){{1,4}}")
        ref_block_match = pattern.search(text)

        pps = []
        if ref_block_match:
            # print(ref_block_match.group())
            # print('---')
            # print(ref_block_match.group(1))
            # print('---')

            ref_block = ref_block_match.group()
            # ref_block = ref_block_match.group(1)

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

                articles = get_pubmed_articles(t_strip, count=result_limit)

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
                                ## we don't really use this
                                # 'pubmed': 'XXX link for %s XXX' % anchor_index,
                                # 'pubmed': article,
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

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # match_found = True if len(t_w_r_s) > 0 else match_found

    print("-----------\r\n all_refs: %s \r\n-----------" % all_refs)

    if match_found:
        ref_block = "\n\n".join(
            [
                # f"{i+1}. [{ref['title']}]({ref['url']})\n{ref['authors']}.\n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
                # f"###### {i+1}. [{ref['title']}]({ref['url']}) {{ #{i+1} }}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
                f"###### {i + 1}. [{ref['title']}]({ref['url']}) {{#{'%s%s' % (i + 1, anchor_random_string)}}}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
                for i, ref in enumerate(all_refs)
            ]
        )
        ## Replace (and remove) references from raw AI response
        # text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", f"### References:\n{ref_block}\n", text)
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
        # term=query,
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


class BannerViewSet(viewsets.ModelViewSet):
    queryset = Banner.objects.all()
    serializer_class = BannerSerializer
    parser_classes = (MultiPartParser, FormParser)

    def get_serializer_context(self):
        # ensure .context['request'] is available to BannerSerializer
        context = super().get_serializer_context()
        context["request"] = self.request
        return context

    @action(detail=True, methods=["post"])
    def view(self, request, pk=None):
        banner = self.get_object()
        ip = request.META.get("REMOTE_ADDR", "")
        try:
            country = GeoIP2().country(ip).get("country_name", "Unknown")
        except:
            country = "Unknown"
        stat, _ = BannerStat.objects.get_or_create(
            banner=banner, date=timezone.now().date(), country=country
        )
        stat.views = F("views") + 1
        stat.save()
        return Response({"status": "view recorded"})

    ## MOD
    @action(detail=True, methods=["post"])
    def click(self, request, pk=None):
        banner = self.get_object()
        ip = request.META.get("REMOTE_ADDR", "")
        try:
            country = GeoIP2().country(ip).get("country_name", "Unknown")
        except:
            country = "Unknown"
        stat, _ = BannerStat.objects.get_or_create(
            banner=banner, date=timezone.now().date(), country=country
        )
        stat.clicks = F("clicks") + 1
        stat.save()
        return Response({"status": "click recorded"})


class BannerStatViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = BannerStat.objects.all()
    serializer_class = BannerStatSerializer

    def get_queryset(self):
        qs = super().get_queryset()
        banner_id = self.request.query_params.get("banner")
        if banner_id:
            qs = qs.filter(banner__id=banner_id)
        return qs


@api_view(["GET"])
def users_view(request):
    ## MOD : change superuser password
    # uzzer = CustomUser.objects.get(id=8)
    # uzzer.email = 'dentalnovikai@gmail.com'
    # uzzer.set_password('NovikAI1971')
    # uzzer.is_staff = True
    # uzzer.is_superuser = True
    # uzzer.save()

    # uzzer = CustomUser.objects.get(id=6)
    # uzzer.set_password('jjjound14/155C')
    # uzzer.username = 'katseye'
    # uzzer.dob = '1988-05-13'
    # uzzer.country = 'France'
    # uzzer.city = 'Nice'
    # uzzer.save()

    # Fetch all from the CustomUser table
    qs = CustomUser.objects.all()
    serializer = UserListSerializer(qs, many=True)
    return Response({"users": serializer.data}, status=status.HTTP_200_OK)


@api_view(["GET"])
def export_users_csv(request):
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


@api_view(["POST"])
def toggle_user_active_status(request):
    # print('received user id %s' % request.data.get('id'))

    try:
        eyed = request.data.get("id")

        uzzer = CustomUser.objects.get(pk=eyed)
        uzzer.is_active = 1 if uzzer.is_active == 0 else 0
        uzzer.save()

        # print('user with id %s updated' % eyed)

        return Response(
            {"status": "user status set to %s successfully" % uzzer.is_active}
        )
    except Exception:
        return Response({"status": "unable to toggle user active status"})


@api_view(["DELETE"])
def user_delete(request):
    # print('received user id %s' % request.data.get('id'))

    try:
        eyed = request.data.get("id")

        uzzer = CustomUser.objects.get(pk=eyed)
        uzzer.delete()

        # print('user with id %s updated' % eyed)

        return Response({"status": "user deleted successfully"})
    except Exception:
        return Response({"status": "unable to deleted user"})


@api_view(["GET"])
def clear_sessions(request):
    ## Remove collected prompts
    pcs = PatientContext.objects.all()
    kount = pcs.count()
    pcs.delete()

    return Response({"status": "ok", "message": "cleared %d rows" % (kount)})
