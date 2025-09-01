"""
Admin views for user management and system administration.
"""

import csv
from datetime import datetime, timedelta

from django.db.models import Count, Q
from django.http import HttpResponse
from django.utils import timezone
from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework import status
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAdminUser
from rest_framework.response import Response

from .models import Conversation, CustomUser
from .serializers import UserSerializer


@extend_schema(
    tags=["Admin"],
    description="Get all users with filtering and pagination",
    parameters=[
        OpenApiParameter(name="page", type=int, description="Page number"),
        OpenApiParameter(name="page_size", type=int, description="Items per page"),
        OpenApiParameter(name="search", type=str, description="Search term"),
        OpenApiParameter(
            name="is_active", type=bool, description="Filter by active status"
        ),
        OpenApiParameter(
            name="is_staff", type=bool, description="Filter by staff status"
        ),
        OpenApiParameter(name="country", type=str, description="Filter by country"),
    ],
)
@api_view(["GET"])
@permission_classes([IsAdminUser])
def get_users_view(request):
    """Get all users with filtering and pagination for admin."""
    try:
        # Get query parameters
        page = int(request.GET.get("page", 1))
        page_size = int(request.GET.get("page_size", 20))
        search = request.GET.get("search", "")
        is_active = request.GET.get("is_active")
        is_staff = request.GET.get("is_staff")
        country = request.GET.get("country")

        # Build query
        users = CustomUser.objects.all()

        # Apply filters
        if search:
            users = users.filter(
                Q(username__icontains=search)
                | Q(email__icontains=search)
                | Q(first_name__icontains=search)
                | Q(last_name__icontains=search)
            )

        if is_active is not None:
            users = users.filter(is_active=is_active == "true")

        if is_staff is not None:
            users = users.filter(is_staff=is_staff == "true")

        if country:
            users = users.filter(country=country)

        # Get total count before pagination
        total_count = users.count()

        # Order by date joined (newest first)
        users = users.order_by("-date_joined")

        # Apply pagination
        start = (page - 1) * page_size
        end = start + page_size
        users_page = users[start:end]

        # Serialize users
        users_data = []
        for user in users_page:
            conversation_count = Conversation.objects.filter(user=user).count()
            user_data = UserSerializer(user).data
            user_data["conversation_count"] = conversation_count
            users_data.append(user_data)

        return Response(
            {
                "users": users_data,
                "total": total_count,
                "page": page,
                "page_size": page_size,
                "total_pages": (total_count + page_size - 1) // page_size,
            }
        )
    except Exception as e:
        return Response(
            {"error": f"Failed to fetch users: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Get a single user's detailed information")
@api_view(["GET"])
@permission_classes([IsAdminUser])
def get_user_detail_view(request, user_id):
    """Get detailed information about a specific user."""
    try:
        user = CustomUser.objects.get(id=user_id)
        conversations = Conversation.objects.filter(user=user).order_by("-updated_at")[
            :5
        ]

        user_data = UserSerializer(user).data
        user_data["recent_conversations"] = [
            {
                "id": str(conv.id),
                "title": conv.title,
                "message_count": len(conv.messages),
                "updated_at": conv.updated_at,
            }
            for conv in conversations
        ]
        user_data["total_conversations"] = Conversation.objects.filter(
            user=user
        ).count()

        return Response(user_data)
    except CustomUser.DoesNotExist:
        return Response({"error": "User not found"}, status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        return Response(
            {"error": f"Failed to fetch user details: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Toggle user's active status")
@api_view(["POST"])
@permission_classes([IsAdminUser])
def toggle_user_status_view(request):
    """Toggle a user's active status."""
    try:
        user_id = request.data.get("id")
        if not user_id:
            return Response(
                {"error": "User ID is required"}, status=status.HTTP_400_BAD_REQUEST
            )

        user = CustomUser.objects.get(id=user_id)

        # Prevent deactivating superusers
        if user.is_superuser and user.is_active:
            return Response(
                {"error": "Cannot deactivate superuser accounts"},
                status=status.HTTP_403_FORBIDDEN,
            )

        user.is_active = not user.is_active
        user.save()

        return Response(
            {
                "message": f"User {'activated' if user.is_active else 'deactivated'} successfully",
                "is_active": user.is_active,
            }
        )
    except CustomUser.DoesNotExist:
        return Response({"error": "User not found"}, status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        return Response(
            {"error": f"Failed to toggle user status: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Update user information")
@api_view(["PATCH"])
@permission_classes([IsAdminUser])
def update_user_view(request, user_id):
    """Update user information as admin."""
    try:
        user = CustomUser.objects.get(id=user_id)

        # Fields that can be updated by admin
        allowed_fields = [
            "email",
            "first_name",
            "last_name",
            "dob",
            "phone",
            "occupation",
            "country",
            "state",
            "city",
            "is_staff",
        ]

        for field in allowed_fields:
            if field in request.data:
                setattr(user, field, request.data[field])

        user.save()

        return Response(
            {"message": "User updated successfully", "user": UserSerializer(user).data}
        )
    except CustomUser.DoesNotExist:
        return Response({"error": "User not found"}, status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        return Response(
            {"error": f"Failed to update user: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Delete a user account")
@api_view(["DELETE"])
@permission_classes([IsAdminUser])
def delete_user_view(request):
    """Delete a user account permanently."""
    try:
        user_id = request.data.get("id")
        if not user_id:
            return Response(
                {"error": "User ID is required"}, status=status.HTTP_400_BAD_REQUEST
            )

        user = CustomUser.objects.get(id=user_id)

        # Prevent deleting superusers
        if user.is_superuser:
            return Response(
                {"error": "Cannot delete superuser accounts"},
                status=status.HTTP_403_FORBIDDEN,
            )

        username = user.username
        user.delete()

        return Response({"message": f"User {username} deleted successfully"})
    except CustomUser.DoesNotExist:
        return Response({"error": "User not found"}, status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        return Response(
            {"error": f"Failed to delete user: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Export users to CSV")
@api_view(["GET"])
@permission_classes([IsAdminUser])
def export_users_csv_view(request):
    """Export all users to CSV file."""
    try:
        # Create the HttpResponse object with CSV header
        response = HttpResponse(content_type="text/csv")
        response["Content-Disposition"] = (
            f'attachment; filename="users_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv"'
        )

        writer = csv.writer(response)
        writer.writerow(
            [
                "ID",
                "Username",
                "Email",
                "First Name",
                "Last Name",
                "Date of Birth",
                "Phone",
                "Occupation",
                "Country",
                "State",
                "City",
                "Is Active",
                "Is Staff",
                "Is Superuser",
                "Date Joined",
                "Last Login",
                "Profile Completed",
                "Conversations Count",
            ]
        )

        users = CustomUser.objects.all().order_by("id")
        for user in users:
            conv_count = Conversation.objects.filter(user=user).count()
            writer.writerow(
                [
                    user.id,
                    user.username,
                    user.email,
                    user.first_name,
                    user.last_name,
                    user.dob,
                    user.phone,
                    user.occupation,
                    user.country,
                    user.state,
                    user.city,
                    user.is_active,
                    user.is_staff,
                    user.is_superuser,
                    user.date_joined.strftime("%Y-%m-%d %H:%M:%S")
                    if user.date_joined
                    else "",
                    user.last_login.strftime("%Y-%m-%d %H:%M:%S")
                    if user.last_login
                    else "",
                    user.profile_completed,
                    conv_count,
                ]
            )

        return response
    except Exception as e:
        return Response(
            {"error": f"Failed to export users: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Get dashboard statistics")
@api_view(["GET"])
@permission_classes([IsAdminUser])
def get_dashboard_stats_view(request):
    """Get statistics for admin dashboard."""
    try:
        now = timezone.now()
        thirty_days_ago = now - timedelta(days=30)
        seven_days_ago = now - timedelta(days=7)

        # User statistics
        total_users = CustomUser.objects.count()
        active_users = CustomUser.objects.filter(is_active=True).count()
        new_users_month = CustomUser.objects.filter(
            date_joined__gte=thirty_days_ago
        ).count()
        new_users_week = CustomUser.objects.filter(
            date_joined__gte=seven_days_ago
        ).count()

        # Users by country
        users_by_country = (
            CustomUser.objects.values("country")
            .annotate(count=Count("id"))
            .order_by("-count")[:10]
        )

        # Recent registrations
        recent_users = CustomUser.objects.order_by("-date_joined")[:5]

        # Conversation statistics
        total_conversations = Conversation.objects.count()
        conversations_month = Conversation.objects.filter(
            created_at__gte=thirty_days_ago
        ).count()

        # Users with most conversations
        top_users = CustomUser.objects.annotate(
            conv_count=Count("conversation")
        ).order_by("-conv_count")[:5]

        return Response(
            {
                "user_stats": {
                    "total": total_users,
                    "active": active_users,
                    "inactive": total_users - active_users,
                    "new_this_month": new_users_month,
                    "new_this_week": new_users_week,
                },
                "conversation_stats": {
                    "total": total_conversations,
                    "this_month": conversations_month,
                    "avg_per_user": round(total_conversations / total_users, 2)
                    if total_users > 0
                    else 0,
                },
                "users_by_country": [
                    {"country": item["country"] or "Unknown", "count": item["count"]}
                    for item in users_by_country
                ],
                "recent_registrations": [
                    {
                        "id": user.id,
                        "username": user.username,
                        "email": user.email,
                        "date_joined": user.date_joined,
                    }
                    for user in recent_users
                ],
                "top_users": [
                    {
                        "id": user.id,
                        "username": user.username,
                        "email": user.email,
                        "conversation_count": user.conv_count,
                    }
                    for user in top_users
                ],
            }
        )
    except Exception as e:
        return Response(
            {"error": f"Failed to fetch dashboard stats: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )


@extend_schema(tags=["Admin"], description="Make a user a staff member")
@api_view(["POST"])
@permission_classes([IsAdminUser])
def make_staff_view(request, user_id):
    """Grant or revoke staff privileges for a user."""
    try:
        user = CustomUser.objects.get(id=user_id)
        is_staff = request.data.get("is_staff", True)

        user.is_staff = is_staff
        user.save()

        return Response(
            {
                "message": f"User {'granted' if is_staff else 'revoked'} staff privileges successfully",
                "is_staff": user.is_staff,
            }
        )
    except CustomUser.DoesNotExist:
        return Response({"error": "User not found"}, status=status.HTTP_404_NOT_FOUND)
    except Exception as e:
        return Response(
            {"error": f"Failed to update staff status: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
