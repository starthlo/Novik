import requests
from django.conf import settings
from django.contrib.auth import authenticate, get_user_model
from django.shortcuts import redirect
from drf_spectacular.utils import OpenApiResponse, extend_schema
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework_simplejwt.tokens import RefreshToken
from social_core.exceptions import AuthAlreadyAssociated, MissingBackend
from social_django.utils import load_backend, load_strategy

from .schemas.auth_schemas import (
    GoogleAuthTokenRequestSerializer,
    LoginRequestSerializer,
    LoginResponseSerializer,
    SocialAuthErrorResponseSerializer,
    SocialAuthSuccessResponseSerializer,
)
from .schemas.openai_schemas import ErrorResponseSerializer
from .serializers import RegisterSerializer, UserSerializer

User = get_user_model()


@extend_schema(
    tags=["Authentication"],
    description="Register a new user",
    request=RegisterSerializer,
    responses={
        201: OpenApiResponse(
            {
                "type": "object",
                "properties": {
                    "user": {"type": "object"},
                    "refresh": {"type": "string"},
                    "access": {"type": "string"},
                    "message": {"type": "string"},
                },
            },
            description="User registered successfully",
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
    },
)
@api_view(["POST"])
@permission_classes([AllowAny])
def register_view(request):
    """
    Register a new user and return access and refresh tokens.

    This endpoint creates a new user account and returns JWT tokens for authentication.
    All required fields must be provided in the request body.
    """
    serializer = RegisterSerializer(data=request.data)
    if serializer.is_valid():
        user = serializer.save()
        refresh = RefreshToken.for_user(user)

        return Response(
            {
                "user": UserSerializer(user).data,
                "refreshToken": str(refresh),
                "accessToken": str(refresh.access_token),
                "message": "User registered successfully",
            },
            status=status.HTTP_201_CREATED,
        )

    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@extend_schema(
    tags=["Authentication"],
    description="Authenticate user and get tokens",
    request=LoginRequestSerializer,
    responses={
        200: OpenApiResponse(
            response=LoginResponseSerializer, description="Login successful"
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
        401: OpenApiResponse(
            response=ErrorResponseSerializer, description="Invalid credentials"
        ),
    },
)
@api_view(["POST"])
@permission_classes([AllowAny])
def login_view(request):
    """
    Authenticate user and return access and refresh tokens.

    This endpoint authenticates a user with email and password,
    and returns JWT tokens for subsequent API calls.
    """
    email = request.data.get("email")
    password = request.data.get("password")

    if not email or not password:
        return Response(
            {"error": "Please provide both email and password"},
            status=status.HTTP_400_BAD_REQUEST,
        )

    user = authenticate(request, username=email, password=password)

    if user is None:
        return Response(
            {
                "error": "Invalid credentials. Please check your email/username and password."
            },
            status=status.HTTP_401_UNAUTHORIZED,
        )

    refresh = RefreshToken.for_user(user)

    return Response(
        {
            "user": UserSerializer(user).data,
            "refreshToken": str(refresh),
            "accessToken": str(refresh.access_token),
        }
    )


# Google OAuth endpoints
@extend_schema(
    tags=["Authentication", "OAuth"],
    description="Redirect to Google for authentication",
    responses={
        302: OpenApiResponse(description="Redirect to Google OAuth"),
    },
)
@api_view(["GET"])
@permission_classes([AllowAny])
def google_login(request):
    """
    Redirect to Google for authentication.

    This endpoint initiates the OAuth flow by redirecting the user
    to Google's authentication page.
    """
    redirect_uri = f"{request.scheme}://{request.get_host()}/api/auth/google/callback/"
    return redirect(
        f"https://accounts.google.com/o/oauth2/auth?client_id={settings.SOCIAL_AUTH_GOOGLE_OAUTH2_KEY}&redirect_uri={redirect_uri}&response_type=code&scope=email%20profile"
    )


@extend_schema(
    tags=["Authentication", "OAuth"],
    description="Process Google callback and authenticate user",
    parameters=[
        {"name": "code", "in": "query", "required": True, "schema": {"type": "string"}}
    ],
    responses={
        302: OpenApiResponse(description="Redirect with tokens"),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
        500: OpenApiResponse(
            response=ErrorResponseSerializer, description="Server error"
        ),
    },
)
@api_view(["GET"])
@permission_classes([AllowAny])
def google_callback(request):
    """
    Process Google callback and authenticate user.

    This endpoint handles the callback from Google OAuth,
    authenticates the user, and generates JWT tokens.
    """
    code = request.GET.get("code", "")
    if not code:
        return Response(
            {"error": "Code not provided"}, status=status.HTTP_400_BAD_REQUEST
        )

    try:
        # Load the Django strategy and Google backend
        strategy = load_strategy(request)
        backend = load_backend(
            strategy=strategy, name="google-oauth2", redirect_uri=None
        )

        # Complete the authentication process
        auth_info = backend.auth_complete(request)

        # Get or create the user
        user = backend.do_auth(access_token=auth_info["access_token"])

        # Generate JWT tokens
        refresh = RefreshToken.for_user(user)
        access_token = str(refresh.access_token)
        refresh_token = str(refresh)

        # Return the tokens in the response or create a redirect with tokens
        return redirect(
            f"{settings.SOCIAL_AUTH_LOGIN_REDIRECT_URL}?accessToken={access_token}&refreshToken={refresh_token}"
        )

    except (MissingBackend, AuthAlreadyAssociated) as e:
        return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)
    except Exception as e:
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


@extend_schema(
    tags=["Authentication", "OAuth"],
    description="Endpoint for successful social authentication",
    parameters=[
        {
            "name": "accessToken",
            "in": "query",
            "required": True,
            "schema": {"type": "string"},
        },
        {
            "name": "refreshToken",
            "in": "query",
            "required": True,
            "schema": {"type": "string"},
        },
    ],
    responses={
        200: OpenApiResponse(
            response=SocialAuthSuccessResponseSerializer,
            description="Authentication successful",
        ),
        400: OpenApiResponse(
            response=ErrorResponseSerializer, description="Bad request"
        ),
    },
)
@api_view(["GET"])
@permission_classes([AllowAny])
def social_auth_success(request):
    """
    Endpoint for successful authentication - provides JWT tokens.

    This endpoint is called after successful social authentication,
    and returns the JWT tokens for API access.
    """
    accessToken = request.GET.get("accessToken", "")
    refreshToken = request.GET.get("refreshToken", "")

    if not accessToken or not refreshToken:
        return Response(
            {"error": "Authentication failed"}, status=status.HTTP_400_BAD_REQUEST
        )

    return Response(
        {
            "accessToken": accessToken,
            "refreshToken": refreshToken,
            "message": "Successfully authenticated with Google",
        }
    )


@extend_schema(
    tags=["Authentication", "OAuth"],
    description="Endpoint for authentication errors",
    parameters=[
        {
            "name": "error",
            "in": "query",
            "required": False,
            "schema": {"type": "string"},
        },
    ],
    responses={
        400: OpenApiResponse(
            response=SocialAuthErrorResponseSerializer,
            description="Authentication error",
        ),
    },
)
@api_view(["GET"])
@permission_classes([AllowAny])
def social_auth_error(request):
    """
    Endpoint for authentication errors.

    This endpoint is called when social authentication fails,
    and returns error details.
    """
    return Response(
        {
            "error": "Authentication failed",
            "details": request.GET.get("error", "Unknown error"),
        },
        status=status.HTTP_400_BAD_REQUEST,
    )


@extend_schema(
    tags=["Authentication", "OAuth"],
    description="Exchange Google OAuth token for JWT tokens",
    request=GoogleAuthTokenRequestSerializer,
    responses={
        200: OpenApiResponse(
            response=LoginResponseSerializer, description="Token exchange successful"
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
@permission_classes([AllowAny])
def google_auth_token(request):
    """
    Exchange Google OAuth token for JWT tokens.

    This endpoint takes a Google OAuth token, verifies it with Google,
    and returns JWT tokens for API access.
    """
    google_token = request.data.get("token")
    if not google_token:
        return Response(
            {"error": "Token not provided"}, status=status.HTTP_400_BAD_REQUEST
        )

    try:
        # Verify the Google token and get user info
        url = f"https://www.googleapis.com/oauth2/v3/tokeninfo?id_token={google_token}"
        response = requests.get(url)
        if response.status_code != 200:
            return Response(
                {"error": "Invalid token"}, status=status.HTTP_400_BAD_REQUEST
            )

        user_data = response.json()
        email = user_data.get("email")

        if not email:
            return Response(
                {"error": "Email not found in token"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # Get or create a user with this email
        try:
            user = User.objects.get(email=email)
        except User.DoesNotExist:
            # Create a new user
            username = email.split("@")[0]
            # Ensure username is unique
            base_username = username
            i = 1
            while User.objects.filter(username=username).exists():
                username = f"{base_username}{i}"
                i += 1

            user = User.objects.create(
                username=username,
                email=email,
                first_name=user_data.get("given_name", ""),
                last_name=user_data.get("family_name", ""),
                # Set agree_to_terms to True as they're coming from Google OAuth
                agree_to_terms=True,
                # Profile is NOT complete for Google OAuth users - they need to fill additional info
                profile_completed=False,
            )
            # Set an unusable password as they'll login via Google
            user.set_unusable_password()
            user.save()

        # Generate JWT tokens
        refresh = RefreshToken.for_user(user)

        return Response(
            {
                "user": UserSerializer(user).data,
                "accessToken": str(refresh.access_token),
                "refreshToken": str(refresh),
            }
        )

    except Exception as e:
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
