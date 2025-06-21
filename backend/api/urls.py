from django.urls import path
from rest_framework.routers import DefaultRouter
from . import views

router = DefaultRouter()
router.register('banners', views.BannerViewSet,      basename='banner')
router.register('banner-stats', views.BannerStatViewSet, basename='bannerstat')

urlpatterns = [
    path('register/',       views.register_view),
    path('login/',          views.login_view),
    path('google-login/',   views.google_login_view),
    path('dashboard/',      views.openai_response_view),
    path('dashboard/pdf/',  views.openai_pdf_response_view),
    path('users/',          views.users_view),
    path('users/export/',   views.export_users_csv),
    # path('user/toggle-active-status/<int:id>',   views.toggle_user_active_status),
    path('user/toggle-active-status/',   views.toggle_user_active_status),
    path('user/trash/',   views.user_delete),
    path('clear-sessions/',   views.clear_sessions),
]

urlpatterns += router.urls
