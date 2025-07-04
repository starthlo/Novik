from django.contrib.auth import get_user_model
from django.contrib.auth.backends import ModelBackend


class EmailBackend(ModelBackend):
    """
    Custom authentication backend to authenticate with email and password.
    """

    def authenticate(self, request, username=None, password=None, **kwargs):
        """
        Authenticate a user by email and password.

        This method overrides the default behavior to use email for authentication.
        It can still accept a username parameter for compatibility, but treats it as an email.
        """
        UserModel = get_user_model()

        # If the username field contains an email, try to authenticate by email
        try:
            if "@" in username:
                user = UserModel.objects.get(email=username)
                # Authenticate with the password
                if user.check_password(password):
                    return user
            else:
                return super().authenticate(
                    request, username=username, password=password, **kwargs
                )
        except UserModel.DoesNotExist:
            # No user with this email or username exists
            return None

        return None
