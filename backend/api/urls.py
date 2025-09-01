from django.urls import path
from rest_framework.routers import DefaultRouter
from rest_framework_simplejwt.views import TokenRefreshView

from . import admin_views, auth_views, password_reset_views, user_views, views

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
    # Password Reset URLs
    path(
        "auth/forgot-password/",
        password_reset_views.forgot_password_view,
        name="forgot-password",
    ),
    path(
        "auth/validate-reset-token/",
        password_reset_views.validate_reset_token_view,
        name="validate-reset-token",
    ),
    path(
        "auth/reset-password/",
        password_reset_views.reset_password_view,
        name="reset-password",
    ),
    # User Profile URLs
    path("user/profile/", user_views.get_profile_view, name="get-profile"),
    path("user/profile/update/", user_views.update_profile_view, name="update-profile"),
    path(
        "user/change-password/", user_views.change_password_view, name="change-password"
    ),
    # Novik Backend URLs
    path("patient/assistant/", views.patient_assistant_view),
    # Admin URLs
    path("admin/users/", admin_views.get_users_view, name="admin-users"),
    path(
        "admin/users/<int:user_id>/",
        admin_views.get_user_detail_view,
        name="admin-user-detail",
    ),
    path(
        "admin/users/<int:user_id>/update/",
        admin_views.update_user_view,
        name="admin-user-update",
    ),
    path(
        "admin/users/<int:user_id>/staff/",
        admin_views.make_staff_view,
        name="admin-make-staff",
    ),
    path(
        "admin/users/toggle-status/",
        admin_views.toggle_user_status_view,
        name="admin-toggle-status",
    ),
    path("admin/users/delete/", admin_views.delete_user_view, name="admin-delete-user"),
    path(
        "admin/users/export/",
        admin_views.export_users_csv_view,
        name="admin-export-users",
    ),
    path(
        "admin/dashboard/stats/",
        admin_views.get_dashboard_stats_view,
        name="admin-dashboard-stats",
    ),
]

urlpatterns += router.urls
