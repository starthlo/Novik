from rest_framework import serializers
from django.contrib.auth.password_validation import validate_password
from django.contrib.auth import authenticate
from google.auth.transport import requests
from google.oauth2.id_token import verify_oauth2_token
from .models import CustomUser, Banner, BannerStat

class UserSerializer(serializers.ModelSerializer):
    confirm_password = serializers.CharField(write_only=True)

    class Meta:
        model = CustomUser
        fields = ['username', 'email', 'password', 'confirm_password', 'dob', 'phone', 'occupation', 'country', 'state', 'city', 'agree_to_terms', 'receive_info']

    def validate_email(self, value):
        if CustomUser.objects.filter(email=value).exists():
            raise serializers.ValidationError("This email is already registered.")
        return value
    
    def validate(self, data):
        if data['password'] != data['confirm_password']:
            raise serializers.ValidationError("Passwords do not match")
        validate_password(data['password'])
        return data

    def create(self, validated_data):
        validated_data.pop('confirm_password')
        user = CustomUser.objects.create_user(**validated_data)
        return user

class LoginSerializer(serializers.Serializer):
    email = serializers.EmailField()
    password = serializers.CharField(write_only=True)

    def validate(self, data):
        email = data.get('email')
        password = data.get('password')

        user = authenticate(username=email, password=password)
        if not user:
            raise serializers.ValidationError('Invalid email or password')

        if not user.is_active:
            raise serializers.ValidationError('Account is inactive')

        data['user'] = user
        return data

class GoogleLoginSerializer(serializers.Serializer):
    token = serializers.CharField()

    def validate(self, data):
        token = data.get('token')

        try:
            idinfo = verify_oauth2_token(token, requests.Request())
            if idinfo['iss'] not in ['accounts.google.com', 'https://accounts.google.com']:
                raise serializers.ValidationError('Invalid issuer')

            google_id = idinfo['sub']
            email = idinfo.get('email')
            name = idinfo.get('name')

            # Create or update user based on Google info
            user, created = CustomUser.objects.get_or_create(
                email=email,
                defaults={'username': name}
            )
            data['user'] = user

            return data

        except Exception as e:
            raise serializers.ValidationError('Invalid Google token')
        
class UserListSerializer(serializers.ModelSerializer):
    """
    Serializer for listing users with non‑sensitive fields.
    """
    class Meta:
        model = CustomUser
        fields = (
            'id',
            'username',
            'email',
            'dob',
            'phone',
            'occupation',
            'country',
            'state',
            'city',
            'is_staff',
            'is_superuser',
            'is_active',
            'date_joined',
        )
        read_only_fields = fields

from .models import Banner, BannerStat

class BannerSerializer(serializers.ModelSerializer):
    image_url = serializers.SerializerMethodField()

    class Meta:
        model = Banner
        fields = [
            'id', 'title', 'image', 'image_url',
            'link', 'code', 'is_active',
            'created_at', 'updated_at'
        ]
        read_only_fields = ['id', 'image_url', 'created_at', 'updated_at']

    def get_image_url(self, obj):
        """
        Return a fully‐qualified URL for the uploaded image.
        """
        request = self.context.get('request')
        if obj.image and hasattr(obj.image, 'url'):
            # e.g. http://localhost:8000/media/banners/yourfile.jpg
            return request.build_absolute_uri(obj.image.url)
        return None

    def validate(self, data):
        # Ensure at least one of image or code is provided
        if not data.get('image') and not data.get('code'):
            raise serializers.ValidationError(
                "Either an image or HTML code must be provided"
            )
        return data

class BannerStatSerializer(serializers.ModelSerializer):
    class Meta:
        model = BannerStat
        fields = [
            'id',
            'banner',
            'date',
            'country',
            'views',
            'clicks',
        ]
        read_only_fields = fields