from rest_framework import serializers


class ToggleUserStatusRequestSerializer(serializers.Serializer):
    """
    Request serializer for toggling user active status.
    """

    id = serializers.IntegerField(
        required=True, help_text="User ID to toggle active status"
    )


class ToggleUserStatusResponseSerializer(serializers.Serializer):
    """
    Response serializer for toggle user status endpoint.
    """

    status = serializers.CharField(help_text="Status message")


class UserDeleteRequestSerializer(serializers.Serializer):
    """
    Request serializer for deleting a user.
    """

    id = serializers.IntegerField(required=True, help_text="User ID to delete")


class UserDeleteResponseSerializer(serializers.Serializer):
    """
    Response serializer for user delete endpoint.
    """

    status = serializers.CharField(help_text="Status message")


class ClearSessionsResponseSerializer(serializers.Serializer):
    """
    Response serializer for clearing sessions.
    """

    status = serializers.CharField(help_text="Status of the operation")
    message = serializers.CharField(help_text="Details about the operation")
