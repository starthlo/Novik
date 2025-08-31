"""
Schema serializers for password reset functionality.
"""

from django.contrib.auth.password_validation import validate_password
from rest_framework import serializers


class ForgotPasswordRequestSerializer(serializers.Serializer):
    """Request serializer for forgot password."""

    email = serializers.EmailField(
        required=True, help_text="Email address of the account to reset password for"
    )


class ForgotPasswordResponseSerializer(serializers.Serializer):
    """Response serializer for forgot password."""

    success = serializers.BooleanField()
    message = serializers.CharField()


class ValidateResetTokenRequestSerializer(serializers.Serializer):
    """Request serializer for validating reset token."""

    uid = serializers.CharField(required=True, help_text="User ID encoded in base64")
    token = serializers.CharField(required=True, help_text="Password reset token")


class ValidateResetTokenResponseSerializer(serializers.Serializer):
    """Response serializer for validating reset token."""

    valid = serializers.BooleanField()
    message = serializers.CharField()
    email = serializers.EmailField(required=False)


class ResetPasswordRequestSerializer(serializers.Serializer):
    """Request serializer for resetting password."""

    uid = serializers.CharField(required=True, help_text="User ID encoded in base64")
    token = serializers.CharField(required=True, help_text="Password reset token")
    password = serializers.CharField(
        required=True,
        write_only=True,
        style={"input_type": "password"},
        help_text="New password",
    )
    confirm_password = serializers.CharField(
        required=True,
        write_only=True,
        style={"input_type": "password"},
        help_text="Confirm new password",
    )

    def validate(self, attrs):
        """Validate that passwords match and meet requirements."""
        password = attrs.get("password")
        confirm_password = attrs.get("confirm_password")

        if password != confirm_password:
            raise serializers.ValidationError(
                {"confirm_password": "Passwords do not match."}
            )

        # Validate password strength
        validate_password(password)

        return attrs


class ResetPasswordResponseSerializer(serializers.Serializer):
    """Response serializer for resetting password."""

    success = serializers.BooleanField()
    message = serializers.CharField()
    errors = serializers.DictField(required=False)
