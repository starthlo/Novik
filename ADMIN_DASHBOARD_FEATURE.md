# Admin Dashboard Feature

## Overview
This document describes the comprehensive admin dashboard feature that has been added to the Novik application. The admin dashboard provides administrators with powerful tools to manage users, monitor system activity, and maintain the application.

## Features

### 1. Admin Dashboard Overview
- **Real-time Statistics**: View key metrics at a glance
  - Total users, active users, inactive users
  - New registrations (weekly and monthly)
  - Total conversations and average per user
  - User distribution by country
  - Most active users
- **Quick Actions**: Export user data, navigate to management pages
- **Visual Analytics**: Charts and progress bars for data visualization

### 2. User Management System
- **Advanced User Search**: 
  - Search by username, email, first name, or last name
  - Filter by status (active/inactive)
  - Filter by role (staff/regular users)
  - Filter by country
- **User Actions**:
  - View detailed user information
  - Edit user profiles (email, name, phone, location, etc.)
  - Toggle user active/inactive status
  - Grant or revoke staff privileges
  - Delete user accounts (with confirmation)
- **Bulk Operations**:
  - Export all users to CSV
  - Pagination for large datasets
- **User Details View**:
  - Complete profile information
  - Recent conversation history
  - Activity statistics

### 3. Security Features
- **Role-Based Access Control**:
  - Only staff and superusers can access admin features
  - Superuser accounts cannot be deleted or deactivated
  - Current password verification for sensitive operations
- **Audit Trail**:
  - Track user registrations
  - Monitor conversation activity
  - View login history

## Implementation Details

### Backend Changes

#### New Files Created
1. **`backend/api/admin_views.py`** - Complete admin API endpoints:
   - `get_users_view` - List users with filtering and pagination
   - `get_user_detail_view` - Get detailed user information
   - `toggle_user_status_view` - Activate/deactivate users
   - `update_user_view` - Edit user information
   - `delete_user_view` - Remove user accounts
   - `make_staff_view` - Grant/revoke staff privileges
   - `export_users_csv_view` - Export users to CSV
   - `get_dashboard_stats_view` - Dashboard analytics

#### Updated Files
1. **`backend/api/urls.py`** - Added admin endpoints:
   ```python
   path("admin/users/", admin_views.get_users_view)
   path("admin/users/<int:user_id>/", admin_views.get_user_detail_view)
   path("admin/users/<int:user_id>/update/", admin_views.update_user_view)
   path("admin/users/<int:user_id>/staff/", admin_views.make_staff_view)
   path("admin/users/toggle-status/", admin_views.toggle_user_status_view)
   path("admin/users/delete/", admin_views.delete_user_view)
   path("admin/users/export/", admin_views.export_users_csv_view)
   path("admin/dashboard/stats/", admin_views.get_dashboard_stats_view)
   ```

### Frontend Changes

#### New Files Created
1. **`frontend/src/services/adminService.ts`** - Admin API service layer
   - Complete TypeScript interfaces for admin data
   - API methods for all admin operations
   - Error handling and authentication

2. **`frontend/src/pages/AdminDashboard.tsx`** - Main dashboard page
   - Statistics cards with icons
   - User distribution charts
   - Recent activity feeds
   - Quick action buttons

3. **`frontend/src/pages/AdminUsers.tsx`** - User management page
   - Advanced data table with sorting
   - Search and filter controls
   - Inline editing capabilities
   - Modal dialogs for detailed views

#### Updated Files
1. **`frontend/src/App.tsx`** - Added admin routes:
   ```tsx
   <Route path="/admin/dashboard" element={<AdminDashboard />} />
   <Route path="/admin/users" element={<AdminUsers />} />
   ```

2. **`frontend/src/components/Common/Header.tsx`** - Added admin navigation:
   - Admin menu items in header
   - User dropdown with admin options
   - Mobile-responsive admin navigation

## API Endpoints

### Authentication Required
All admin endpoints require JWT authentication with staff/superuser privileges.

### Endpoints Reference

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/admin/dashboard/stats/` | Get dashboard statistics |
| GET | `/api/admin/users/` | List all users with filters |
| GET | `/api/admin/users/{id}/` | Get user details |
| PATCH | `/api/admin/users/{id}/update/` | Update user information |
| POST | `/api/admin/users/{id}/staff/` | Toggle staff status |
| POST | `/api/admin/users/toggle-status/` | Toggle active status |
| DELETE | `/api/admin/users/delete/` | Delete user account |
| GET | `/api/admin/users/export/` | Export users to CSV |

### Query Parameters

**GET /api/admin/users/**
- `page` (int): Page number for pagination
- `page_size` (int): Number of items per page
- `search` (string): Search term for username/email
- `is_active` (bool): Filter by active status
- `is_staff` (bool): Filter by staff status
- `country` (string): Filter by country

## User Interface

### Admin Dashboard Page
- **URL**: `/admin/dashboard`
- **Access**: Staff and superusers only
- **Features**:
  - 4 main statistics cards
  - 3 data visualization sections
  - Export and navigation actions
  - Responsive grid layout

### User Management Page
- **URL**: `/admin/users`
- **Access**: Staff and superusers only
- **Features**:
  - Searchable user table
  - Multiple filter options
  - Inline action buttons
  - Edit and detail modals
  - Pagination controls

## Testing

### Manual Testing
1. Create a superuser account:
   ```bash
   cd backend
   python manage.py createsuperuser
   ```

2. Run the test script:
   ```bash
   python test_admin_endpoints.py
   ```

3. Access the admin dashboard:
   - Login with admin credentials
   - Navigate to `/admin/dashboard`
   - Test all features and actions

### Test Coverage
- User CRUD operations
- Search and filtering
- Permission checks
- CSV export functionality
- Dashboard statistics accuracy

## Security Considerations

1. **Access Control**:
   - Only authenticated staff/superusers can access admin features
   - Superuser accounts have additional protections
   - Session-based permission checks

2. **Data Protection**:
   - Sensitive user data is filtered in responses
   - Passwords are never exposed in API responses
   - CSV exports exclude sensitive fields

3. **Audit Logging**:
   - All admin actions should be logged (future enhancement)
   - User modifications tracked with timestamps

## Future Enhancements

1. **Analytics Dashboard**:
   - Advanced charts and graphs
   - Time-based analytics
   - User behavior patterns
   - System performance metrics

2. **Bulk Operations**:
   - Mass email sending
   - Bulk user updates
   - Batch imports from CSV

3. **Audit System**:
   - Complete action logging
   - Admin activity reports
   - Security event monitoring

4. **Advanced Filters**:
   - Date range filters
   - Multiple selection filters
   - Saved filter presets

5. **Real-time Updates**:
   - WebSocket integration
   - Live user activity feed
   - Real-time statistics

## Troubleshooting

### Common Issues

1. **403 Forbidden Error**:
   - Ensure user has staff/superuser privileges
   - Check JWT token validity
   - Verify CORS settings

2. **Empty Dashboard**:
   - Check database connectivity
   - Verify user data exists
   - Check API endpoint responses

3. **CSV Export Not Working**:
   - Verify file permissions
   - Check browser download settings
   - Ensure proper content-type headers

## Conclusion

The admin dashboard provides a comprehensive solution for managing the Novik application. With its intuitive interface and powerful features, administrators can efficiently monitor system health, manage users, and maintain the application's integrity.
