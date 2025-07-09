import logging
import os
import random
import re
import string
import time

import openai
from django.contrib.gis.geoip2 import GeoIP2
from django.db.models import F
from django.utils import timezone
from dotenv import load_dotenv
from drf_spectacular.utils import OpenApiParameter, OpenApiResponse, extend_schema
from rest_framework import status, viewsets
from rest_framework.decorators import action, api_view, permission_classes
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from .models import Banner, BannerStat, Conversation
from .pubmed import get_articles
from .schemas.openai_schemas import (
    ErrorResponseSerializer,
    OpenAIPDFRequestSerializer,
    OpenAIResponseSerializer,
)
from .serializers import (
    BannerSerializer,
    BannerStatSerializer,
    ConversationSerializer,
)

load_dotenv()

logger = logging.getLogger(__name__)

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
ASSISTANT_ID = os.environ.get("OPENAI_ASSISTANT_ID")


client = openai.OpenAI(api_key=OPENAI_API_KEY)


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


@extend_schema(tags=["Conversations"])
class ConversationViewSet(viewsets.ModelViewSet):
    """
    ViewSet for managing user conversations.

    Provides CRUD operations for conversations:
    - List all conversations for the current user
    - Create a new conversation
    - Retrieve a specific conversation
    - Update a conversation's title
    - Delete a conversation
    """

    serializer_class = ConversationSerializer
    permission_classes = [IsAuthenticated]

    def get_queryset(self):
        """
        Filter to return only conversations owned by the current user.
        """
        return Conversation.objects.filter(user=self.request.user)

    @extend_schema(
        description="Clear all messages from a conversation",
        responses={200: ConversationSerializer},
    )
    @action(detail=True, methods=["post"])
    def clear(self, request, pk=None):
        """Clear all messages from the conversation."""
        conversation = self.get_object()

        conversation.messages = []
        conversation.save()

        # Return the updated conversation
        serializer = self.get_serializer(conversation)
        return Response(serializer.data)


@extend_schema(
    tags=["Patients"],
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
        file = None
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

            if text:
                file = {"file_name": pdf_file.name, "text": text}
        else:
            if not content:
                return Response(
                    {"error": "Message content is required when no PDF is provided"},
                    status=status.HTTP_400_BAD_REQUEST,
                )

        conversation = Conversation.objects.filter(user=request.user).first()
        conversation.add_message("user", content, file=file)

        thread = client.beta.threads.create()
        for message in conversation.messages:
            content_to_process = message["content"]
            if message["file"]:
                content_to_process += f"\n\n{message['file']['text']}"
            client.beta.threads.messages.create(
                thread_id=thread.id, role=message["role"], content=content_to_process
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
        text = messages.data[0].content[0].text.value
        annotations = messages.data[0].content[0].text.annotations

        anchor_random_string = "".join(
            random.choices(string.ascii_letters + string.digits, k=5)
        )

        # Iterate over the annotations and add footnotes
        citations = []
        for index, annotation in enumerate(annotations):
            anchor_id = f"{index + 1}{anchor_random_string}"
            ## Replace the text with a footnote
            text = text.replace(
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

        if citations:
            ref_index = text.find("### References")
            if ref_index != -1:
                text = text[:ref_index]
            else:
                ref_index = text.find("**References")
                if ref_index != -1:
                    text = text[:ref_index]

            text += (
                "\r\n"
                + "### References:  "
                + "\r\n"
                + "  \r\n".join(citations)
                + "\r\n"
            )

        text = re.sub(
            r"^.*https?:\/\/(?:www\.|go\.)?drugbank\.(?:com|ca)\/drugs\/DB\d+.*\n?",
            "",
            text,
            flags=re.MULTILINE,
        )

        enriched_response = enrich_response(text, anchor_random_string)

        conversation.add_message("assistant", enriched_response)
        return Response({"message": enriched_response}, status=status.HTTP_200_OK)

    except Exception as e:
        logger.error(f"patient_assistant_view error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


def enrich_response(text: str, anchor_random_string: str) -> str:
    """
    Enrich the response text with medication links and PubMed references.

    This function:
    1. Adds drugbank.com links to medication names
    2. Extracts reference titles from the AI response
    3. Fetches PubMed articles for those references
    4. Looks for preoperative/postoperative guidelines and enhances them with references
    5. Reformats the references section with proper formatting

    Args:
        text: The AI response text to enrich
        anchor_random_string: A random string used for creating unique anchor IDs

    Returns:
        The enriched text with medication links and PubMed references
    """
    # Dictionary of medications and their DrugBank IDs
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
        "Articaine": "DB09009",
        "Epinephrine": "DB00668",
        "Lidocaine": "DB00281",
        "Mepivacaine": "DB00961",
        "Bupivacaine": "DB00297",
    }

    # Add drugbank links to medication names
    for med, dbid in medications.items():
        text = re.sub(
            rf"(?<!\[)({med})(?![\w\s]*\])",
            rf"[{med}](https://go.drugbank.com/drugs/{dbid})",
            text,
        )

    # Extract reference titles from the AI response
    ref_block_match = re.search(r"### References:?(.*\n){1,9}", text)
    raw_titles = []
    if ref_block_match:
        ref_block = ref_block_match.group()
        raw_titles = re.findall(r"\[\d+\]\s*(.*)\n", ref_block)

    # Clean up the titles
    strip_titles = list(filter(str.strip, raw_titles))
    non_empty_titles = list(filter(None, strip_titles))
    titles = list(set(non_empty_titles))

    # Get PubMed articles for the reference titles
    result_limit = 2
    match_found = False
    all_refs = []

    for title in titles:
        title = title.strip()
        if not title:
            continue

        articles = get_articles(title, count=result_limit)
        if articles:
            match_found = True
            all_refs.extend(articles)
        else:
            all_refs.append({"title": title})

    # Process preoperative/postoperative sections
    refmes = ["Preoperative Precautions", "Postoperative Guidelines"]
    anchor_index = len(all_refs)
    t_w_r_s = []

    for refme in refmes:
        # Find bulleted lists under these headings
        pattern = re.compile(rf"{re.escape(refme)}.*\n?([-\s+|\d\.\s+].*\n){{1,4}}")
        ref_block_match = pattern.search(text)

        if ref_block_match:
            ref_block = ref_block_match.group()
            # Extract bullet points
            pps = re.findall(r"[-|\d\.]\s+(.*)\n", ref_block)

            result_limit = 3

            for eye, pp_content in enumerate(pps):
                t_strip = pp_content.strip()

                # Skip empty items or known problematic patterns
                if not t_strip or t_strip == "Anesthetics":
                    continue

                # Clean up the content
                t_strip = t_strip.replace("**", "")

                # Fetch PubMed articles for this content
                articles = get_articles(t_strip, count=result_limit)

                if articles:
                    match_found = True
                    all_refs.extend(articles)

                    for article in articles:
                        anchor_index += 1
                        t_w_r_s.append(
                            {
                                "refme": refme,
                                "eye": eye,
                                "content": pp_content,
                                "anchor_index": anchor_index,
                            }
                        )

    # Add superscripts to referenced content
    for t_w_r in t_w_r_s:
        a_c = t_w_r["content"]
        a_i = t_w_r["anchor_index"]
        a_f = f"{a_i}{anchor_random_string}"
        a_c_a_f = f"{a_c} [^{a_i}^](#{a_f})"
        text = text.replace(a_c, a_c_a_f)

    # Format and add references to the output
    if match_found:
        ref_block = "\n\n".join(
            [
                f"###### {i + 1}. [{ref['title']}]({ref['url']}) {{#{i + 1}{anchor_random_string}}}  \n{ref['authors']}.  \n{ref['journal']}. {ref['year']};{ref.get('volume', '')}({ref.get('issue', '')}):{ref.get('pages', '')}. {ref.get('doi', '')}"
                for i, ref in enumerate(all_refs)
            ]
        )

        # Replace original references with our enhanced ones
        text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", "", text)

        # Add new PubMed references
        text = f"{text}  \r\n\r\n-------\n### PubMed References\r\n{ref_block}  \r\n"
    else:
        text = f"{text}  \r\n\r\n**No recent PubMed references were found that specifically relate to these factors.**\r\n"

    return text
