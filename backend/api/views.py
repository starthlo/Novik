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

        anchor = "".join(random.choices(string.ascii_letters + string.digits, k=5))

        # Iterate over the annotations and add footnotes
        citations = []
        for index, annotation in enumerate(annotations):
            ## Replace the text with a footnote
            text = text.replace(annotation.text, "", 1)
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
            for ref_marker in ["### References", "**References"]:
                ref_index = text.find(ref_marker)
                if ref_index != -1:
                    text = text[:ref_index]
                    break

            text = f"{text}\r\n\n### References:  \r\n{chr(10).join(f'  {citation}' for citation in citations)}\r\n"

        enriched_response = enrich_response(text + "\n", anchor)

        conversation.add_message("assistant", enriched_response)
        return Response({"message": enriched_response}, status=status.HTTP_200_OK)

    except Exception as e:
        logger.error(f"patient_assistant_view error: {str(e)}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


def enrich_response(text: str, anchor: str) -> str:
    """
    Enrich the response text with medication links and PubMed references.

    This function:
    1. Adds drugbank.com links to medication names
    2. Extracts PubMed search terms from the AI response
    3. Fetches PubMed articles for those search terms
    4. Adds numbered superscript references to medical statements
    5. Reformats the references section with proper formatting

    Args:
        text: The AI response text to enrich
        anchor: A random string used for creating unique anchor IDs

    Returns:
        The enriched text with medication links and PubMed references
    """
    # Dictionary of medications and their DrugBank IDs
    medications = {}
    file_path = os.path.join(os.path.dirname(__file__), "drugs/molecules.txt")

    if os.path.exists(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line == ",":
                    continue
                match = re.match(r'"(.+?)"\s*:\s*"(.+?)",?', line)
                if match:
                    name, dbid = match.groups()
                    medications[name] = dbid
    else:
        # Optionally log a warning or raise an exception
        medications = {}

    # Add drugbank links to medication names
    for med, dbid in medications.items():
        escaped_med = re.escape(med)
        text = re.sub(
            rf"(?<!\[)({escaped_med})(?![\w\s]*\])",
            rf"[{med}](https://go.drugbank.com/drugs/{dbid})",
            text,
        )

    # Initialize variables
    all_refs = []
    match_found = False

    # Extract PubMed search terms from the text
    pubmed_search_pattern = (
        r"(\*\*)?PubMed Search Terms:(\*\*)?\s*\n((?:\d+\.\s+.*?\n)*)"
    )
    pubmed_matches = re.findall(pubmed_search_pattern, text, re.MULTILINE | re.DOTALL)
    # Extract only the search terms part (third capture group)
    pubmed_matches = [match[2] for match in pubmed_matches]

    anchor_index = 0

    for pubmed_block in pubmed_matches:
        # Extract individual search terms with their ordering numbers
        search_items = re.findall(
            r"(\d+)\.\s+(.+?)(?=\n\d+\.|\Z)", pubmed_block, re.DOTALL
        )

        for order_num, search_term in search_items:
            # Clean up the search term
            search_term = search_term.strip().replace("\n", " ")
            search_term = re.sub(r"\s+", " ", search_term)

            # Fetch PubMed articles for this search term
            articles = get_articles(search_term, count=3)
            time.sleep(1)

            if articles:
                match_found = True
                # Add ordering number to each article for reference
                for article in articles:
                    anchor_index += 1
                    text = re.sub(
                        rf"\[\^{order_num}\]",
                        f"[^{anchor_index}^](#{anchor_index}{anchor}) [^{order_num}]",
                        text,
                    )
                all_refs.extend(articles)

            superscript_pattern = rf"\[\^{order_num}\]"
            text = re.sub(superscript_pattern, "", text)

    # Format and add references to the output
    if match_found and all_refs:
        # Remove existing references section
        text = re.sub(r"### References:[\s\S]*?(?=\n###|\Z)", "", text)
        text = re.sub(
            r"(\*\*)?PubMed Search Terms:(\*\*)?\s*\n(?:\d+\.\s+\([^)]+\)\s*\n)*(?:\n(?=\*\*[A-Z]\.)|(?=\n###)|$)",
            "",
            text,
            flags=re.MULTILINE,
        )

        # Create new references section
        ref_block = "\n\n".join(
            [format_reference(i, ref, anchor) for i, ref in enumerate(all_refs)]
        )

        # Add new PubMed references
        text = f"{text.rstrip()}\n\n\n### PubMed References\n\n{ref_block}\n"
    else:
        # Clean up PubMed search terms sections even if no articles found
        text = re.sub(
            r"(\*\*)?PubMed Search Terms:(\*\*)?\s*\n(?:\d+\.\s+\([^)]+\)\s*\n)*(?:\n(?=\*\*[A-Z]\.)|(?=\n###)|$)",
            "",
            text,
            flags=re.MULTILINE,
        )
        # text = f"{text.rstrip()}\n\n**No recent PubMed references were found that specifically relate to these factors.**\n"

    return text


def format_reference(i, ref, anchor):
    """Format a single reference entry"""
    # Basic info
    title_line = (
        f"###### {i + 1}. [{ref['title']}]({ref['url']}) {{#{i + 1}{anchor}}}  \n"
    )
    authors_line = f"{ref['authors']}.  \n"

    # Journal and year
    journal_year = f"{ref['journal']}. {ref['year']}"

    # Volume, issue, pages
    citation_parts = []
    if ref.get("volume"):
        citation_parts.append(f";{ref['volume']}")
    if ref.get("issue"):
        citation_parts.append(f"({ref['issue']})")
    if ref.get("pages"):
        citation_parts.append(f":{ref['pages']}")

    citation_info = "".join(citation_parts) + ". "

    # DOI
    doi_info = f"{ref['doi']}" if ref.get("doi") else ""

    return title_line + authors_line + journal_year + citation_info + doi_info
