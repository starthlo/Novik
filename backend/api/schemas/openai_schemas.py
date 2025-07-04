from rest_framework import serializers


class OpenAIResponseRequestSerializer(serializers.Serializer):
    """
    Request serializer for OpenAI response endpoint.
    """

    message = serializers.CharField(
        required=True, help_text="The message to send to OpenAI"
    )
    sessionId = serializers.CharField(
        required=True, help_text="Session ID for context tracking"
    )
    createdBy = serializers.IntegerField(
        required=False, default=0, help_text="User ID of the creator"
    )


class OpenAIResponseSerializer(serializers.Serializer):
    """
    Response serializer for OpenAI response.
    """

    message = serializers.CharField(help_text="Response from OpenAI")


class OpenAIPDFRequestSerializer(serializers.Serializer):
    """
    Request serializer for OpenAI PDF response endpoint.
    """

    pdf = serializers.FileField(required=True, help_text="PDF file to analyze")
    message = serializers.CharField(
        required=True, help_text="Additional message/query about the PDF"
    )
    sessionId = serializers.CharField(
        required=True, help_text="Session ID for context tracking"
    )
    createdBy = serializers.IntegerField(
        required=False, default=0, help_text="User ID of the creator"
    )


class ErrorResponseSerializer(serializers.Serializer):
    """
    Error response serializer.
    """

    error = serializers.CharField(help_text="Error message")
