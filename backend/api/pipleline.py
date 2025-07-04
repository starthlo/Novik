from django.contrib.auth import get_user_model

User = get_user_model()


def create_user(strategy, details, backend, user=None, *args, **kwargs):
    """
    Custom pipeline function to ensure CustomUser fields are properly set.
    """
    if user:
        return {"is_new": False}

    if not details.get("email"):
        return None

    # Get data from social provider
    email = details.get("email")
    first_name = details.get("first_name", "")
    last_name = details.get("last_name", "")

    # Create a new user with the required fields for CustomUser
    user = User.objects.create(
        username=kwargs.get("username"),
        email=email,
        first_name=first_name,
        last_name=last_name,
        agree_to_terms=True,  # Set to True for social logins
    )

    return {"is_new": True, "user": user}
