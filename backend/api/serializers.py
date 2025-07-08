from django.contrib.auth import get_user_model
from django.contrib.auth.password_validation import validate_password
from rest_framework import serializers

from .models import Banner, BannerStat, Conversation

User = get_user_model()


class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = (
            "id",
            "username",
            "email",
            "first_name",
            "last_name",
            "dob",
            "phone",
            "occupation",
            "country",
            "state",
            "city",
        )
        read_only_fields = ("id",)


class RegisterSerializer(serializers.ModelSerializer):
    password = serializers.CharField(
        write_only=True, required=True, validators=[validate_password]
    )
    password2 = serializers.CharField(write_only=True, required=True)

    class Meta:
        model = User
        fields = (
            "username",
            "password",
            "password2",
            "email",
            "first_name",
            "last_name",
            "dob",
            "phone",
            "occupation",
            "country",
            "state",
            "city",
            "agree_to_terms",
            "receive_info",
        )
        extra_kwargs = {
            "email": {"required": True},
            "first_name": {"required": False},
            "last_name": {"required": False},
            "agree_to_terms": {"required": True},
        }

    def validate(self, attrs):
        if attrs["password"] != attrs["password2"]:
            raise serializers.ValidationError(
                {"password": "Password fields didn't match."}
            )

        if not attrs.get("agree_to_terms", False):
            raise serializers.ValidationError(
                {"agree_to_terms": "You must agree to the terms and conditions."}
            )

        return attrs

    def create(self, validated_data):
        validated_data.pop("password2")

        user = User.objects.create(
            username=validated_data["username"],
            email=validated_data["email"],
            first_name=validated_data.get("first_name", ""),
            last_name=validated_data.get("last_name", ""),
            dob=validated_data.get("dob", None),
            phone=validated_data.get("phone", ""),
            occupation=validated_data.get("occupation", ""),
            country=validated_data.get("country", ""),
            state=validated_data.get("state", ""),
            city=validated_data.get("city", ""),
            agree_to_terms=validated_data.get("agree_to_terms", False),
            receive_info=validated_data.get("receive_info", False),
        )

        user.set_password(validated_data["password"])
        user.save()

        return user


class BannerSerializer(serializers.ModelSerializer):
    image_url = serializers.SerializerMethodField()

    class Meta:
        model = Banner
        fields = [
            "id",
            "title",
            "image",
            "image_url",
            "link",
            "code",
            "is_active",
            "created_at",
            "updated_at",
        ]
        read_only_fields = ["id", "image_url", "created_at", "updated_at"]

    def get_image_url(self, obj):
        """
        Return a fully‐qualified URL for the uploaded image.
        """
        request = self.context.get("request")
        if obj.image and hasattr(obj.image, "url"):
            # e.g. http://localhost:8000/media/banners/yourfile.jpg
            return request.build_absolute_uri(obj.image.url)
        return None

    def validate(self, data):
        # Ensure at least one of image or code is provided
        if not data.get("image") and not data.get("code"):
            raise serializers.ValidationError(
                "Either an image or HTML code must be provided"
            )
        return data


class BannerStatSerializer(serializers.ModelSerializer):
    class Meta:
        model = BannerStat
        fields = [
            "id",
            "banner",
            "date",
            "country",
            "views",
            "clicks",
        ]
        read_only_fields = fields


class ConversationSerializer(serializers.ModelSerializer):
    """
    Serializer for Conversation model.

    Provides fields for creating, retrieving, and updating conversations.
    Messages field is read-only as it should be modified using the model's add_message method.
    """

    user = serializers.PrimaryKeyRelatedField(
        read_only=True, default=serializers.CurrentUserDefault()
    )

    class Meta:
        model = Conversation
        fields = [
            "id",
            "user",
            "title",
            "messages",
            "created_at",
            "updated_at",
        ]
        read_only_fields = ["id", "messages", "created_at", "updated_at"]

    def create(self, validated_data):
        # Set the user to the current authenticated user
        validated_data["user"] = self.context["request"].user
        return super().create(validated_data)
