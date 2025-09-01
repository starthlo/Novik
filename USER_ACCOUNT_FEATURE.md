# User Account Management Feature

## Overview
This document describes the user account management feature that has been added to the Novik dental assistant application. The feature allows registered users to view and edit their profile information and change their password through a dedicated account settings page.

## Implementation Details

### Backend Changes

#### 1. New API Endpoints (`backend/api/user_views.py`)
Three new endpoints have been created:

- **GET `/api/user/profile/`** - Retrieve current user's profile
- **PUT/PATCH `/api/user/profile/update/`** - Update user profile information
- **POST `/api/user/change-password/`** - Change user password

All endpoints require authentication (JWT token in Authorization header).

#### 2. User Model Fields
The existing `CustomUser` model includes the following editable fields:
- `email` - User's email address
- `first_name` - First name
- `last_name` - Last name
- `dob` - Date of birth
- `phone` - Phone number
- `occupation` - Professional occupation
- `country` - Country
- `state` - State/Province
- `city` - City

Note: `username` cannot be changed after registration for security reasons.

### Frontend Changes

#### 1. New Components/Pages

**`frontend/src/pages/Account.tsx`**
- Complete account management page with three sections:
  1. Profile Information - View/edit personal details
  2. Change Password - Update account password
  3. Account Actions - Sign out option

**`frontend/src/services/userService.ts`**
- Service layer for API communication
- Handles profile fetching, updating, and password changes

#### 2. Header Component Updates (`frontend/src/components/Common/Header.tsx`)

**Desktop View:**
- Added user account icon in top-right corner (circular person icon)
- Clicking the icon opens a dropdown menu with:
  - User's name and email
  - "My Account" option (navigates to `/account`)
  - "Sign Out" option

**Mobile View:**
- Added user section in mobile drawer
- Shows user name and email
- "My Account" menu item
- Updated "Sign Out" button styling

#### 3. Routing
- Added `/account` route as a protected route in `App.tsx`
- Only accessible to authenticated users

## User Interface Features

### Account Icon
- Located in the top-right corner of the header
- Shows a circular person icon
- Green color matching the Novik brand theme
- Opens dropdown menu on click

### Account Page Features

1. **Profile Editing**
   - Toggle between view and edit modes
   - Edit icon to enable editing
   - Save/Cancel buttons when in edit mode
   - Real-time validation

2. **Password Change**
   - Separate section for security
   - Current password verification
   - New password confirmation
   - Password visibility toggle
   - Minimum 8 character requirement

3. **Visual Feedback**
   - Success/error alerts for all actions
   - Loading states during API calls
   - Form validation messages
   - Responsive design for all screen sizes

## Security Considerations

1. **Authentication Required**
   - All profile endpoints require valid JWT token
   - Automatic redirect to login if not authenticated

2. **Password Security**
   - Current password must be verified before changes
   - Password strength validation
   - Session maintained after password change

3. **Data Protection**
   - Username cannot be changed (prevents account takeover)
   - Email changes require valid format
   - All data transmitted over HTTPS in production

## Testing

To test the implementation:

1. **Backend Testing:**
   ```bash
   cd backend
   source .venv/bin/activate
   python manage.py runserver
   ```

2. **Frontend Testing:**
   ```bash
   cd frontend
   npm run dev
   ```

3. **Manual Testing Steps:**
   - Register a new account or login with existing credentials
   - Click the account icon in the top-right corner
   - Select "My Account" from the dropdown
   - Test profile editing functionality
   - Test password change functionality
   - Verify all error states and validations

## File Structure

```
backend/
├── api/
│   ├── user_views.py        # New user profile endpoints
│   └── urls.py              # Updated with new routes

frontend/
├── src/
│   ├── pages/
│   │   └── Account.tsx      # Account management page
│   ├── services/
│   │   └── userService.ts   # User API service
│   ├── components/
│   │   └── Common/
│   │       └── Header.tsx   # Updated with user menu
│   └── App.tsx              # Updated routing
```

## Future Enhancements

Potential improvements for future iterations:

1. **Profile Picture Upload** - Allow users to upload avatar images
2. **Email Verification** - Require email confirmation for changes
3. **Two-Factor Authentication** - Add 2FA for enhanced security
4. **Account Deletion** - Allow users to delete their accounts
5. **Activity Log** - Show recent account activity
6. **Email Notifications** - Notify users of account changes
7. **Social Account Linking** - Link/unlink Google accounts
8. **Session Management** - View and manage active sessions

## Conclusion

The user account management feature provides a complete solution for users to manage their profile information and account security. The implementation follows best practices for both UX design and security, with a clean, intuitive interface that matches the Novik brand aesthetic.
