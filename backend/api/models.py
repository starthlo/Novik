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
        return f"{self.banner.title} – {self.country or 'Unknown'} – {self.date}"


class PatientContext(models.Model):
    session_id = models.CharField(max_length=36)
    content = models.TextField(blank=True, null=True)
    created_by = models.IntegerField(default=0)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.session_id} @ {self.created_at} – {self.content}"
