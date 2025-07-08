from django.urls import path
from rest_framework.routers import DefaultRouter
from rest_framework_simplejwt.views import TokenRefreshView

from . import auth_views, views

router = DefaultRouter()
router.register("banners", views.BannerViewSet, basename="banner")
router.register("banner-stats", views.BannerStatViewSet, basename="bannerstat")
router.register("conversations", views.ConversationViewSet, basename="conversation")

urlpatterns = [
    # Authentication URLs
    path("auth/register/", auth_views.register_view, name="register"),
    path("auth/login/", auth_views.login_view, name="login"),
    path("auth/token/refresh/", TokenRefreshView.as_view(), name="token-refresh"),
    # Google OAuth URLs
    path("auth/google/", auth_views.google_login, name="google-login"),
    path("auth/google/callback/", auth_views.google_callback, name="google-callback"),
    path("auth/social/success/", auth_views.social_auth_success, name="social-success"),
    path("auth/social/error/", auth_views.social_auth_error, name="social-error"),
    path("auth/google/token/", auth_views.google_auth_token, name="google-auth-token"),
    # Novik Backend URLs
    path("patient/assistant/", views.patient_assistant_view),
]

urlpatterns += router.urls
