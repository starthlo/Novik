from rest_framework import serializers

from api.serializers import UserSerializer


class LoginRequestSerializer(serializers.Serializer):
    """
    Request serializer for login endpoint.
    """

    email = serializers.CharField(required=True, help_text="Email")
    password = serializers.CharField(
        required=True, help_text="Password", style={"input_type": "password"}
    )


class LoginResponseSerializer(serializers.Serializer):
    """
    Response serializer for login endpoint.
    """

    user = UserSerializer()
    refreshToken = serializers.CharField(help_text="JWT refresh token")
    accessToken = serializers.CharField(help_text="JWT access token")


class GoogleAuthTokenRequestSerializer(serializers.Serializer):
    """
    Request serializer for Google Auth token endpoint.
    """

    token = serializers.CharField(required=True, help_text="Google OAuth token")


class SocialAuthSuccessResponseSerializer(serializers.Serializer):
    """
    Response serializer for successful social authentication.
    """

    accessToken = serializers.CharField(help_text="JWT access token")
    refreshToken = serializers.CharField(help_text="JWT refresh token")
    message = serializers.CharField(help_text="Success message")


class SocialAuthErrorResponseSerializer(serializers.Serializer):
    """
    Response serializer for social authentication error.
    """

    error = serializers.CharField(help_text="Error message")
    details = serializers.CharField(help_text="Error details")
