"""
Password reset views for the Novik application.
Handles forgot password and password reset functionality.
"""

import logging

from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.auth.tokens import PasswordResetTokenGenerator
from django.core.mail import EmailMultiAlternatives
from django.utils.encoding import force_bytes, force_str
from django.utils.http import urlsafe_base64_decode, urlsafe_base64_encode
from drf_spectacular.utils import extend_schema
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.response import Response

from .email_templates import (
    add_mailgun_campaign,
    add_mailgun_tags,
    add_mailgun_tracking,
    add_mailgun_variables,
    get_password_reset_confirmation_html,
    get_password_reset_confirmation_text,
    get_password_reset_email_html,
    get_password_reset_email_text,
)
from .schemas.password_reset_schemas import (
    ForgotPasswordRequestSerializer,
    ForgotPasswordResponseSerializer,
    ResetPasswordRequestSerializer,
    ResetPasswordResponseSerializer,
    ValidateResetTokenRequestSerializer,
    ValidateResetTokenResponseSerializer,
)

logger = logging.getLogger(__name__)
User = get_user_model()

# Token generator for password reset
password_reset_token_generator = PasswordResetTokenGenerator()


@extend_schema(
    request=ForgotPasswordRequestSerializer,
    responses={
        200: ForgotPasswordResponseSerializer,
        400: ForgotPasswordResponseSerializer,
        404: ForgotPasswordResponseSerializer,
    },
    tags=["Authentication"],
    operation_id="forgot_password",
    description="Request a password reset link for the given email address.",
)
@api_view(["POST"])
@permission_classes([AllowAny])
def forgot_password_view(request):
    """
    Handle forgot password requests.
    Sends a password reset email to the user if the email exists.
    """
    serializer = ForgotPasswordRequestSerializer(data=request.data)

    if not serializer.is_valid():
        return Response(
            {"success": False, "message": "Invalid email address provided."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    email = serializer.validated_data["email"]

    try:
        user = User.objects.get(email=email)
    except User.DoesNotExist:
        # Return success even if user doesn't exist (security best practice)
        # This prevents email enumeration attacks
        return Response(
            {
                "success": True,
                "message": "If an account exists with this email, you will receive a password reset link.",
            },
            status=status.HTTP_200_OK,
        )

    # Generate password reset token
    uid = urlsafe_base64_encode(force_bytes(user.pk))
    token = password_reset_token_generator.make_token(user)

    # Create reset link
    frontend_url = (
        settings.FRONTEND_URL
        if hasattr(settings, "FRONTEND_URL")
        else "http://localhost:5173"
    )
    reset_link = f"{frontend_url}/reset-password?uid={uid}&token={token}"

    # Prepare email content using templates
    email_subject = "🔐 Password Reset Request - Novik"
    user_name = user.first_name or user.username

    # Get email templates
    html_message = get_password_reset_email_html(user_name, reset_link, 24)
    plain_message = get_password_reset_email_text(user_name, reset_link, 24)

    try:
        # Create email message with both HTML and plain text
        from_email = (
            settings.MAILGUN_FROM_EMAIL
            if hasattr(settings, "MAILGUN_FROM_EMAIL")
            else "Novik <noreply@novik.ai>"
        )

        # Use EmailMultiAlternatives for better email client compatibility
        email_message = EmailMultiAlternatives(
            subject=email_subject,
            body=plain_message,
            from_email=from_email,
            to=[email],
        )

        # Attach HTML version
        email_message.attach_alternative(html_message, "text/html")

        # Add Mailgun-specific features for better tracking and analytics
        if hasattr(settings, 'EMAIL_SERVICE') and settings.EMAIL_SERVICE == 'mailgun':
            # Add tags for categorization (max 3 tags)
            add_mailgun_tags(email_message, ['password-reset', 'security', 'transactional'])
            
            # Enable click and open tracking
            add_mailgun_tracking(email_message, track_opens=True, track_clicks=True)
            
            # Add template variables for analytics
            add_mailgun_variables(email_message, {
                'user_name': user_name,
                'user_email': email,
                'action': 'password_reset_request'
            })
            
            # Add to campaign for reporting
            add_mailgun_campaign(email_message, 'password-reset-q4-2024')

        # Send the email
        email_message.send(fail_silently=False)

        logger.info(f"Password reset email sent to {email}")

        return Response(
            {
                "success": True,
                "message": "If an account exists with this email, you will receive a password reset link.",
            },
            status=status.HTTP_200_OK,
        )

    except Exception as e:
        logger.error(f"Failed to send password reset email to {email}: {str(e)}")
        return Response(
            {
                "success": False,
                "message": "Failed to send password reset email. Please try again later.",
            },
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(
    request=ValidateResetTokenRequestSerializer,
    responses={
        200: ValidateResetTokenResponseSerializer,
        400: ValidateResetTokenResponseSerializer,
    },
    tags=["Authentication"],
    operation_id="validate_reset_token",
    description="Validate a password reset token.",
)
@api_view(["POST"])
@permission_classes([AllowAny])
def validate_reset_token_view(request):
    """
    Validate a password reset token.
    This is called when the user clicks the reset link to verify it's still valid.
    """
    serializer = ValidateResetTokenRequestSerializer(data=request.data)

    if not serializer.is_valid():
        return Response(
            {"valid": False, "message": "Invalid request data."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    uid = serializer.validated_data["uid"]
    token = serializer.validated_data["token"]

    try:
        # Decode the user ID
        user_id = force_str(urlsafe_base64_decode(uid))
        user = User.objects.get(pk=user_id)
    except (TypeError, ValueError, OverflowError, User.DoesNotExist):
        return Response(
            {"valid": False, "message": "Invalid or expired reset link."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    # Check if token is valid
    if not password_reset_token_generator.check_token(user, token):
        return Response(
            {"valid": False, "message": "Invalid or expired reset link."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    return Response(
        {
            "valid": True,
            "message": "Token is valid.",
            "email": user.email,  # Return email to show in the reset form
        },
        status=status.HTTP_200_OK,
    )


@extend_schema(
    request=ResetPasswordRequestSerializer,
    responses={
        200: ResetPasswordResponseSerializer,
        400: ResetPasswordResponseSerializer,
    },
    tags=["Authentication"],
    operation_id="reset_password",
    description="Reset user password with a valid token.",
)
@api_view(["POST"])
@permission_classes([AllowAny])
def reset_password_view(request):
    """
    Reset the user's password with a valid token.
    """
    serializer = ResetPasswordRequestSerializer(data=request.data)

    if not serializer.is_valid():
        return Response(
            {
                "success": False,
                "message": "Invalid data provided.",
                "errors": serializer.errors,
            },
            status=status.HTTP_400_BAD_REQUEST,
        )

    uid = serializer.validated_data["uid"]
    token = serializer.validated_data["token"]
    password = serializer.validated_data["password"]

    try:
        # Decode the user ID
        user_id = force_str(urlsafe_base64_decode(uid))
        user = User.objects.get(pk=user_id)
    except (TypeError, ValueError, OverflowError, User.DoesNotExist):
        return Response(
            {"success": False, "message": "Invalid or expired reset link."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    # Check if token is valid
    if not password_reset_token_generator.check_token(user, token):
        return Response(
            {"success": False, "message": "Invalid or expired reset link."},
            status=status.HTTP_400_BAD_REQUEST,
        )

    # Set the new password
    user.set_password(password)
    user.save()

    logger.info(f"Password successfully reset for user {user.email}")

    # Send confirmation email
    try:
        confirmation_subject = "Password Successfully Reset - Novik"
        user_name = user.first_name or user.username

        # Get confirmation email templates
        html_message = get_password_reset_confirmation_html(user_name)
        plain_message = get_password_reset_confirmation_text(user_name)

        from_email = (
            settings.MAILGUN_FROM_EMAIL
            if hasattr(settings, "MAILGUN_FROM_EMAIL")
            else "Novik <noreply@novik.ai>"
        )

        # Create email message
        email_message = EmailMultiAlternatives(
            subject=confirmation_subject,
            body=plain_message,
            from_email=from_email,
            to=[user.email],
        )

        # Attach HTML version
        email_message.attach_alternative(html_message, "text/html")

        # Add Mailgun-specific features for confirmation email
        if hasattr(settings, 'EMAIL_SERVICE') and settings.EMAIL_SERVICE == 'mailgun':
            # Add tags for categorization
            add_mailgun_tags(email_message, ['password-reset-confirm', 'security', 'notification'])
            
            # Enable tracking (less important for confirmations)
            add_mailgun_tracking(email_message, track_opens=True, track_clicks=False)
            
            # Add template variables
            add_mailgun_variables(email_message, {
                'user_name': user_name,
                'user_email': user.email,
                'action': 'password_reset_completed'
            })
            
            # Add to campaign
            add_mailgun_campaign(email_message, 'security-notifications-2024')

        email_message.send(fail_silently=True)
    except Exception as e:
        logger.error(f"Failed to send confirmation email: {str(e)}")

    return Response(
        {"success": True, "message": "Password has been successfully reset."},
        status=status.HTTP_200_OK,
    )
