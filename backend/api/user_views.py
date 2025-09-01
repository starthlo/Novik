"""
User profile management views.

Provides endpoints for:
- Getting current user profile
- Updating user profile
- Changing password
"""

from django.contrib.auth import update_session_auth_hash
from django.contrib.auth.password_validation import validate_password
from drf_spectacular.utils import extend_schema
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

from .serializers import UserSerializer


@extend_schema(
    tags=["User Profile"],
    description="Get current user profile",
    responses={
        200: UserSerializer,
        401: {"description": "Unauthorized"},
    },
)
@api_view(["GET"])
@permission_classes([IsAuthenticated])
def get_profile_view(request):
    """
    Get the current user's profile information.
    """
    serializer = UserSerializer(request.user)
    return Response(serializer.data)


@extend_schema(
    tags=["User Profile"],
    description="Update current user profile",
    request=UserSerializer,
    responses={
        200: UserSerializer,
        400: {"description": "Validation error"},
        401: {"description": "Unauthorized"},
    },
)
@api_view(["PUT", "PATCH"])
@permission_classes([IsAuthenticated])
def update_profile_view(request):
    """
    Update the current user's profile information.
    """
    serializer = UserSerializer(
        request.user, data=request.data, partial=request.method == "PATCH"
    )

    if serializer.is_valid():
        # Don't allow users to change their username via this endpoint
        if "username" in serializer.validated_data:
            del serializer.validated_data["username"]

        # Check if required fields are now filled to mark profile as complete
        user = request.user
        # required_fields = ['occupation', 'dob', 'phone', 'country', 'state', 'city']

        # # Merge existing data with new data
        # updated_data = {**serializer.validated_data}

        # # Check if all required fields will be filled after update
        # all_required_filled = all(
        #     updated_data.get(field) or getattr(user, field)
        #     for field in required_fields
        # )

        # # Also check if agree_to_terms is True (already set or being set)
        # agree_to_terms = updated_data.get('agree_to_terms', user.agree_to_terms)

        if not user.profile_completed:
            serializer.validated_data["profile_completed"] = True

        serializer.save()
        return Response(serializer.data)

    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@extend_schema(
    tags=["User Profile"],
    description="Change current user password",
    request={
        "type": "object",
        "properties": {
            "current_password": {"type": "string"},
            "new_password": {"type": "string"},
            "confirm_password": {"type": "string"},
        },
        "required": ["current_password", "new_password", "confirm_password"],
    },
    responses={
        200: {"description": "Password changed successfully"},
        400: {"description": "Validation error"},
        401: {"description": "Unauthorized"},
    },
)
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def change_password_view(request):
    """
    Change the current user's password.
    """
    current_password = request.data.get("current_password")
    new_password = request.data.get("new_password")
    confirm_password = request.data.get("confirm_password")

    # Validate that all fields are provided
    if not all([current_password, new_password, confirm_password]):
        return Response(
            {"error": "All password fields are required"},
            status=status.HTTP_400_BAD_REQUEST,
        )

    # Check if the current password is correct
    if not request.user.check_password(current_password):
        return Response(
            {"error": "Current password is incorrect"},
            status=status.HTTP_400_BAD_REQUEST,
        )

    # Check if new passwords match
    if new_password != confirm_password:
        return Response(
            {"error": "New passwords do not match"}, status=status.HTTP_400_BAD_REQUEST
        )

    # Validate the new password
    try:
        validate_password(new_password, request.user)
    except Exception as e:
        return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)

    # Set the new password
    request.user.set_password(new_password)
    request.user.save()

    # Update session to prevent logout
    update_session_auth_hash(request, request.user)

    return Response(
        {"message": "Password changed successfully"}, status=status.HTTP_200_OK
    )
