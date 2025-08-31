"""
Clean and simple email templates for Novik.
Minimalist design focused on readability and clarity.
"""


def get_email_base_html(content: str, preheader: str = "") -> str:
    """
    Base HTML email template with clean, simple design.

    Args:
        content: The main content HTML
        preheader: Preview text shown in email clients

    Returns:
        Complete HTML email with minimalist design
    """

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <title>Novik</title>
    
    <style>
        /* Reset */
        body, table, td, a {{ -webkit-text-size-adjust: 100%; -ms-text-size-adjust: 100%; }}
        table, td {{ mso-table-lspace: 0pt; mso-table-rspace: 0pt; }}
        img {{ -ms-interpolation-mode: bicubic; border: 0; outline: none; text-decoration: none; }}
        
        /* Base Styles */
        body {{
            margin: 0 !important;
            padding: 0 !important;
            background-color: #f5f5f5 !important;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif !important;
        }}
        
        /* Container */
        .email-container {{
            max-width: 600px;
            margin: 40px auto;
            background-color: #ffffff;
            border-radius: 8px;
            overflow: hidden;
        }}
        
        /* Logo */
        .logo {{
            padding: 30px;
            text-align: center;
            border-bottom: 1px solid #e5e5e5;
        }}
        
        .logo-text {{
            color: #333333 !important;
            font-size: 24px;
            font-weight: 600;
            text-decoration: none;
        }}
        
        /* Content */
        .content {{
            padding: 40px 30px;
        }}
        
        h1 {{
            color: #333333;
            font-size: 24px;
            font-weight: 600;
            margin: 0 0 20px 0;
            line-height: 1.4;
        }}
        
        p {{
            color: #555555;
            font-size: 16px;
            line-height: 1.6;
            margin: 0 0 20px 0;
        }}
        
        /* Button */
        .button {{
            display: inline-block;
            padding: 12px 30px;
            background-color: #5c6ac4;
            color: #ffffff !important;
            font-size: 16px;
            font-weight: 500;
            text-decoration: none !important;
            border-radius: 6px;
            margin: 20px 0;
        }}
        
        /* Link */
        .link-box {{
            background-color: #f8f8f8;
            border: 1px solid #e5e5e5;
            border-radius: 6px;
            padding: 15px;
            margin: 20px 0;
            word-break: break-all;
            font-family: monospace;
            font-size: 14px;
            color: #666666;
        }}
        
        /* Alert */
        .alert {{
            background-color: #fff8e1;
            border-left: 3px solid #ffc107;
            padding: 15px;
            margin: 25px 0;
        }}
        
        .alert p {{
            margin: 0;
            font-size: 14px;
        }}
        
        /* Footer */
        .footer {{
            padding: 30px;
            text-align: center;
            border-top: 1px solid #e5e5e5;
        }}
        
        .footer p {{
            color: #999999;
            font-size: 14px;
            margin: 5px 0;
        }}
        
        .footer a {{
            color: #5c6ac4 !important;
            text-decoration: none;
        }}
        
        /* Responsive */
        @media screen and (max-width: 600px) {{
            .email-container {{
                width: 100% !important;
                margin: 0 !important;
                border-radius: 0 !important;
            }}
            
            .content {{
                padding: 30px 20px !important;
            }}
            
            .button {{
                display: block !important;
                width: 100% !important;
                text-align: center !important;
            }}
        }}
    </style>
</head>
<body>
    <!-- Preview Text -->
    <div style="display: none; max-height: 0; overflow: hidden;">
        {preheader}
    </div>
    
    <!-- Email Container -->
    <div class="email-container">
        {content}
    </div>
</body>
</html>"""


def get_password_reset_email_html(
    user_name: str, reset_link: str, valid_hours: int = 24
) -> str:
    """
    Clean and simple password reset email.

    Args:
        user_name: User's name or email
        reset_link: Password reset URL
        valid_hours: Link validity period

    Returns:
        HTML email content
    """
    content = f"""
        <!-- Logo -->
        <div class="logo">
            <a href="https://novik.ai" class="logo-text">Novik</a>
        </div>
        
        <!-- Content -->
        <div class="content">
            <h1>Reset Your Password</h1>
            
            <p>Hi {user_name or "there"},</p>
            
            <p>We received a request to reset your password. Click the button below to create a new password:</p>
            
            <div style="text-align: center;">
                <a href="{reset_link}" class="button">Reset Password</a>
            </div>
            
            <p style="font-size: 14px; color: #888888;">Or copy this link:</p>
            
            <div class="link-box">
                {reset_link}
            </div>
            
            <div class="alert">
                <p><strong>⏱ This link expires in {valid_hours} hours</strong></p>
                <p>For security reasons, you'll need to request a new link after it expires.</p>
            </div>
            
            <p style="font-size: 14px; color: #888888;">
                If you didn't request this, you can safely ignore this email. Your password won't be changed.
            </p>
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>© 2024 Novik AI</p>
            <p>
                <a href="https://novik.ai/help">Help</a> · 
                <a href="https://novik.ai/privacy">Privacy</a> · 
                <a href="mailto:support@novik.ai">Support</a>
            </p>
        </div>"""

    preheader = f"Reset your password - Link expires in {valid_hours} hours"
    return get_email_base_html(content, preheader)


def get_password_reset_email_text(
    user_name: str, reset_link: str, valid_hours: int = 24
) -> str:
    """
    Plain text version of password reset email.

    Args:
        user_name: User's name or email
        reset_link: Password reset URL
        valid_hours: Link validity period

    Returns:
        Plain text email content
    """
    return f"""Reset Your Password

Hi {user_name or "there"},

We received a request to reset your password.

To reset your password, click this link:
{reset_link}

This link will expire in {valid_hours} hours.

If you didn't request this, you can safely ignore this email. Your password won't be changed.

---
© 2024 Novik AI
Support: support@novik.ai
Website: https://novik.ai"""


def get_password_reset_confirmation_html(user_name: str) -> str:
    """
    Clean and simple password reset confirmation email.

    Args:
        user_name: User's name or email

    Returns:
        HTML email content
    """
    content = f"""
        <!-- Logo -->
        <div class="logo">
            <a href="https://novik.ai" class="logo-text">Novik</a>
        </div>
        
        <!-- Content -->
        <div class="content">
            <h1>Password Successfully Reset</h1>
            
            <p>Hi {user_name or "there"},</p>
            
            <p>Your password has been successfully changed. You can now log in with your new password.</p>
            
            <div style="text-align: center;">
                <a href="https://novik.ai/login" class="button">Log In to Novik</a>
            </div>
            
            <div style="margin-top: 40px; padding-top: 20px; border-top: 1px solid #e5e5e5;">
                <p style="font-size: 14px; color: #888888;">
                    <strong>Security Tips:</strong><br>
                    • Use a strong, unique password<br>
                    • Never share your password<br>
                    • Enable two-factor authentication when available
                </p>
                
                <p style="font-size: 14px; color: #888888;">
                    If you didn't make this change, please contact us immediately at 
                    <a href="mailto:security@novik.ai" style="color: #5c6ac4;">security@novik.ai</a>
                </p>
            </div>
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>© 2024 Novik AI</p>
            <p>
                <a href="https://novik.ai/help">Help</a> · 
                <a href="https://novik.ai/privacy">Privacy</a> · 
                <a href="mailto:security@novik.ai">Security</a>
            </p>
        </div>"""

    preheader = "Your password has been successfully reset"
    return get_email_base_html(content, preheader)


def get_password_reset_confirmation_text(user_name: str) -> str:
    """
    Plain text version of password reset confirmation.

    Args:
        user_name: User's name or email

    Returns:
        Plain text email content
    """
    return f"""Password Successfully Reset

Hi {user_name or "there"},

Your password has been successfully changed. You can now log in with your new password.

Log in at: https://novik.ai/login

Security Tips:
• Use a strong, unique password
• Never share your password
• Enable two-factor authentication when available

If you didn't make this change, please contact us immediately at security@novik.ai

---
© 2024 Novik AI
Security: security@novik.ai
Website: https://novik.ai"""


# Mailgun helper functions (simplified)


def add_mailgun_tags(email_message, tags: list):
    """
    Add Mailgun tags for categorization.

    Args:
        email_message: Django EmailMessage object
        tags: List of tags (max 3)
    """
    if tags:
        email_message.extra_headers = email_message.extra_headers or {}
        email_message.extra_headers["X-Mailgun-Tag"] = tags[:3]


def add_mailgun_tracking(
    email_message, track_opens: bool = True, track_clicks: bool = True
):
    """
    Configure Mailgun tracking.

    Args:
        email_message: Django EmailMessage object
        track_opens: Enable open tracking
        track_clicks: Enable click tracking
    """
    email_message.extra_headers = email_message.extra_headers or {}
    email_message.extra_headers["X-Mailgun-Track"] = (
        "yes" if (track_opens or track_clicks) else "no"
    )
    email_message.extra_headers["X-Mailgun-Track-Opens"] = (
        "yes" if track_opens else "no"
    )
    email_message.extra_headers["X-Mailgun-Track-Clicks"] = (
        "yes" if track_clicks else "no"
    )
