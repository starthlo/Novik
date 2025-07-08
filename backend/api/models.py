import uuid
from typing import Literal

from django.contrib.auth.models import AbstractUser
from django.db import models
from django.utils import timezone


class CustomUser(AbstractUser):
    dob = models.DateField(null=True, blank=True)
    phone = models.CharField(max_length=15, blank=True, null=True)
    occupation = models.CharField(max_length=100, blank=True, null=True)
    country = models.CharField(max_length=100, blank=True, null=True)
    state = models.CharField(max_length=100, blank=True, null=True)
    city = models.CharField(max_length=100, blank=True, null=True)
    agree_to_terms = models.BooleanField(default=False)
    receive_info = models.BooleanField(default=False)

    groups = models.ManyToManyField(
        "auth.Group",
        related_name="custom_user_groups",
        blank=True,
        help_text="The groups this user belongs to.",
        verbose_name="groups",
    )
    user_permissions = models.ManyToManyField(
        "auth.Permission",
        related_name="custom_user_permissions",
        blank=True,
        help_text="Specific permissions for this user.",
        verbose_name="user permissions",
    )

    def __str__(self):
        return self.username


class Banner(models.Model):
    title = models.CharField(max_length=255)
    image = models.ImageField(upload_to="banners/", blank=True, null=True)
    link = models.URLField(blank=True, null=True)
    code = models.TextField(blank=True, null=True)
    is_active = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.title


class BannerStat(models.Model):
    banner = models.ForeignKey(Banner, related_name="stats", on_delete=models.CASCADE)
    date = models.DateField(default=timezone.now)
    country = models.CharField(max_length=100, blank=True, null=True)
    views = models.IntegerField(default=0)
    clicks = models.IntegerField(default=0)

    class Meta:
        unique_together = ("banner", "date", "country")

    def __str__(self):
        return f"{self.banner.title} - {self.country or 'Unknown'} - {self.date}"


class Conversation(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(CustomUser, on_delete=models.CASCADE)
    title = models.CharField(max_length=255, blank=True)
    messages = models.JSONField(default=list)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ["-updated_at"]

    def __str__(self):
        return f"{self.title} - {self.user.username} ({self.created_at.strftime('%Y-%m-%d')})"

    def add_message(self, role: Literal["user", "assistant"], content: str):
        """Add a new message to the conversation."""
        message = {"role": role, "content": content}

        if not isinstance(self.messages, list):
            self.messages = []

        self.messages.append(message)
        self.save()

        return message
