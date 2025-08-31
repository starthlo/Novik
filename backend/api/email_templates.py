"""
Professional email templates for Novik - Optimized for Mailgun.
Modern, responsive, and trackable design.
"""

def get_email_base_html(content: str, preheader: str = "", footer_type: str = "default") -> str:
    """
    Base HTML email template with modern design optimized for Mailgun.
    
    Args:
        content: The main content HTML
        preheader: Preview text shown in email clients
        footer_type: Type of footer ('default', 'simple', 'security')
    
    Returns:
        Complete HTML email optimized for Mailgun
    """
    
    footer_html = get_footer_html(footer_type)
    
    return f"""<!DOCTYPE html>
<html lang="en" xmlns="http://www.w3.org/1999/xhtml" xmlns:v="urn:schemas-microsoft-com:vml" xmlns:o="urn:schemas-microsoft-com:office:office">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="x-apple-disable-message-reformatting">
    <meta name="format-detection" content="telephone=no,address=no,email=no,date=no,url=no">
    <meta name="color-scheme" content="light">
    <meta name="supported-color-schemes" content="light">
    <title>Novik - AI Medical Assistant</title>
    
    <!--[if mso]>
    <xml>
        <o:OfficeDocumentSettings>
            <o:AllowPNG/>
            <o:PixelsPerInch>96</o:PixelsPerInch>
        </o:OfficeDocumentSettings>
    </xml>
    <![endif]-->
    
    <style>
        :root {{
            color-scheme: light;
            supported-color-schemes: light;
        }}
        
        * {{
            margin: 0;
            padding: 0;
            -webkit-text-size-adjust: 100%;
            -ms-text-size-adjust: 100%;
        }}
        
        body {{
            margin: 0 !important;
            padding: 0 !important;
            min-width: 100% !important;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif !important;
            font-size: 16px;
            line-height: 1.5;
            color: #2d3748;
            background-color: #f7fafc;
        }}
        
        table {{
            border-spacing: 0 !important;
            border-collapse: collapse !important;
            table-layout: fixed !important;
        }}
        
        /* Main Styles */
        .email-container {{
            max-width: 600px;
            margin: 0 auto;
            background-color: #ffffff;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
        }}
        
        .email-header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 32px 40px;
            text-align: center;
        }}
        
        .logo-text {{
            color: #ffffff !important;
            font-size: 28px;
            font-weight: bold;
            text-decoration: none;
            letter-spacing: -0.5px;
        }}
        
        .email-body {{
            padding: 40px;
            background-color: #ffffff;
        }}
        
        h1 {{
            color: #2d3748;
            font-size: 24px;
            font-weight: 700;
            margin: 0 0 16px 0;
            line-height: 1.4;
        }}
        
        h2 {{
            color: #2d3748;
            font-size: 20px;
            font-weight: 600;
            margin: 24px 0 12px 0;
        }}
        
        p {{
            color: #4a5568;
            font-size: 16px;
            line-height: 24px;
            margin: 0 0 16px 0;
        }}
        
        .text-muted {{
            color: #718096 !important;
            font-size: 14px;
        }}
        
        /* Button Styles */
        .btn-table {{
            margin: 32px 0;
            width: 100%;
        }}
        
        .btn-primary {{
            display: inline-block;
            padding: 14px 32px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #ffffff !important;
            font-size: 16px;
            font-weight: 600;
            text-decoration: none !important;
            border-radius: 6px;
            transition: all 0.3s ease;
        }}
        
        .btn-secondary {{
            display: inline-block;
            padding: 12px 24px;
            background-color: #edf2f7;
            color: #4a5568 !important;
            font-size: 16px;
            font-weight: 600;
            text-decoration: none !important;
            border-radius: 6px;
        }}
        
        /* Alert Styles */
        .alert {{
            padding: 16px;
            border-radius: 6px;
            margin: 24px 0;
        }}
        
        .alert-warning {{
            background-color: #fef5e7;
            border-left: 4px solid #f39c12;
        }}
        
        .alert-info {{
            background-color: #e8f4fd;
            border-left: 4px solid #3498db;
        }}
        
        .alert-success {{
            background-color: #e8f8f5;
            border-left: 4px solid #27ae60;
        }}
        
        .alert-icon {{
            display: inline-block;
            margin-right: 8px;
            font-size: 18px;
            vertical-align: middle;
        }}
        
        /* Code Block */
        .code-block {{
            background-color: #f7fafc;
            border: 1px solid #e2e8f0;
            border-radius: 6px;
            padding: 16px;
            margin: 16px 0;
            font-family: 'Courier New', Courier, monospace;
            font-size: 14px;
            word-break: break-all;
            color: #2d3748;
        }}
        
        /* Divider */
        .divider {{
            height: 1px;
            background-color: #e2e8f0;
            margin: 32px 0;
        }}
        
        /* Footer Styles */
        .email-footer {{
            background-color: #f7fafc;
            padding: 32px 40px;
            text-align: center;
            border-top: 1px solid #e2e8f0;
        }}
        
        .footer-links {{
            margin: 16px 0;
        }}
        
        .footer-links a {{
            color: #667eea !important;
            text-decoration: none;
            font-size: 14px;
            margin: 0 12px;
            font-weight: 500;
        }}
        
        .social-links {{
            margin: 20px 0;
        }}
        
        .social-links a {{
            display: inline-block;
            margin: 0 8px;
        }}
        
        .copyright {{
            color: #718096;
            font-size: 12px;
            margin-top: 16px;
        }}
        
        /* Responsive Styles */
        @media screen and (max-width: 600px) {{
            .email-container {{
                width: 100% !important;
                border-radius: 0 !important;
            }}
            
            .email-body {{
                padding: 24px !important;
            }}
            
            .email-header {{
                padding: 24px !important;
            }}
            
            .email-footer {{
                padding: 24px !important;
            }}
            
            h1 {{
                font-size: 22px !important;
            }}
            
            .btn-primary {{
                display: block !important;
                width: 100% !important;
                text-align: center !important;
            }}
            
            .footer-links a {{
                display: block !important;
                margin: 8px 0 !important;
            }}
        }}
        
        /* Dark Mode Support */
        @media (prefers-color-scheme: dark) {{
            /* Mailgun tracking pixels and links work better without dark mode */
        }}
    </style>
</head>
<body style="margin: 0; padding: 0; background-color: #f7fafc;">
    <!-- Preheader Text -->
    <div style="display: none; font-size: 1px; color: #f7fafc; line-height: 1px; max-height: 0px; max-width: 0px; opacity: 0; overflow: hidden;">
        {preheader}
    </div>
    
    <!-- Email Container -->
    <table role="presentation" cellpadding="0" cellspacing="0" border="0" align="center" width="100%" style="background-color: #f7fafc; padding: 20px 0;">
        <tr>
            <td align="center">
                <table class="email-container" role="presentation" cellpadding="0" cellspacing="0" border="0" width="600" style="background-color: #ffffff; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);">
                    {content}
                    {footer_html}
                </table>
            </td>
        </tr>
    </table>
    
    <!-- Mailgun Tracking Pixel -->
    <img src="https://url.novik.ai/o/{{{{ mailgun_tracking_id }}}}" width="1" height="1" style="display: block;" alt="">
</body>
</html>"""


def get_footer_html(footer_type: str = "default") -> str:
    """Get footer HTML based on type."""
    
    if footer_type == "simple":
        return """
        <tr>
            <td class="email-footer">
                <p class="copyright">© 2024 Novik AI - Empowering Healthcare with AI</p>
            </td>
        </tr>"""
    
    elif footer_type == "security":
        return """
        <tr>
            <td class="email-footer">
                <p style="color: #718096; font-size: 14px; margin-bottom: 16px;">
                    🔒 This is a secure, automated message from Novik
                </p>
                <div class="footer-links">
                    <a href="https://novik.ai/security">Security Center</a>
                    <a href="https://novik.ai/privacy">Privacy Policy</a>
                    <a href="mailto:security@novik.ai">Report Issue</a>
                </div>
                <p class="copyright">© 2024 Novik AI - Secure Healthcare Platform</p>
            </td>
        </tr>"""
    
    else:  # default
        return """
        <tr>
            <td class="email-footer">
                <div class="footer-links">
                    <a href="https://novik.ai">Website</a>
                    <a href="https://novik.ai/help">Help Center</a>
                    <a href="mailto:support@novik.ai">Support</a>
                    <a href="https://novik.ai/privacy">Privacy</a>
                </div>
                <div class="social-links">
                    <a href="https://twitter.com/novik_ai">
                        <img src="https://cdn.novik.ai/email/twitter.png" alt="Twitter" width="24" height="24">
                    </a>
                    <a href="https://linkedin.com/company/novik-ai">
                        <img src="https://cdn.novik.ai/email/linkedin.png" alt="LinkedIn" width="24" height="24">
                    </a>
                </div>
                <p class="copyright">
                    © 2024 Novik AI - Empowering Healthcare with AI<br>
                    <span style="font-size: 11px;">You're receiving this email because you have an account with Novik.</span>
                </p>
            </td>
        </tr>"""


def get_password_reset_email_html(user_name: str, reset_link: str, valid_hours: int = 24) -> str:
    """
    Password reset email optimized for Mailgun with tracking.
    """
    content = f"""
    <!-- Header -->
    <tr>
        <td class="email-header">
            <a href="https://novik.ai" class="logo-text">Novik</a>
        </td>
    </tr>
    
    <!-- Body -->
    <tr>
        <td class="email-body">
            <h1>Reset Your Password</h1>
            
            <p>Hi {user_name or 'there'},</p>
            
            <p>
                We received a request to reset the password for your Novik account. 
                If you made this request, click the button below:
            </p>
            
            <!-- CTA Button -->
            <table role="presentation" cellpadding="0" cellspacing="0" border="0" align="center" class="btn-table">
                <tr>
                    <td align="center">
                        <a href="{reset_link}" class="btn-primary" target="_blank">
                            Reset My Password
                        </a>
                    </td>
                </tr>
            </table>
            
            <p class="text-muted">
                Or copy and paste this link into your browser:
            </p>
            
            <div class="code-block">
                {reset_link}
            </div>
            
            <!-- Warning Alert -->
            <div class="alert alert-warning">
                <span class="alert-icon">⏰</span>
                <strong>This link expires in {valid_hours} hours</strong><br>
                <span style="font-size: 14px; color: #856404;">
                    For security reasons, you'll need to request a new link if this one expires.
                </span>
            </div>
            
            <div class="divider"></div>
            
            <p class="text-muted" style="font-size: 14px;">
                <strong>Didn't request this?</strong><br>
                If you didn't request a password reset, you can safely ignore this email. 
                Your password won't be changed and your account remains secure.
            </p>
        </td>
    </tr>"""
    
    preheader = f"Reset your Novik password - This link expires in {valid_hours} hours"
    return get_email_base_html(content, preheader, "security")


def get_password_reset_email_text(user_name: str, reset_link: str, valid_hours: int = 24) -> str:
    """Plain text version for Mailgun multipart emails."""
    return f"""Reset Your Password

Hi {user_name or 'there'},

We received a request to reset the password for your Novik account.

To reset your password, click this link:
{reset_link}

This link will expire in {valid_hours} hours.

Didn't request this?
If you didn't request a password reset, you can safely ignore this email. Your password won't be changed.

---
© 2024 Novik AI
Support: support@novik.ai
Website: https://novik.ai"""


def get_password_reset_confirmation_html(user_name: str) -> str:
    """
    Password reset confirmation email optimized for Mailgun.
    """
    content = f"""
    <!-- Header -->
    <tr>
        <td class="email-header">
            <a href="https://novik.ai" class="logo-text">Novik</a>
        </td>
    </tr>
    
    <!-- Body -->
    <tr>
        <td class="email-body">
            <!-- Success Alert -->
            <div class="alert alert-success" style="text-align: center;">
                <span style="font-size: 24px;">✅</span><br>
                <strong style="font-size: 18px;">Password Successfully Reset</strong>
            </div>
            
            <h1>Your password has been changed</h1>
            
            <p>Hi {user_name or 'there'},</p>
            
            <p>
                Your Novik account password has been successfully reset. 
                You can now log in with your new password.
            </p>
            
            <!-- CTA Button -->
            <table role="presentation" cellpadding="0" cellspacing="0" border="0" align="center" class="btn-table">
                <tr>
                    <td align="center">
                        <a href="https://novik.ai/login" class="btn-primary" target="_blank">
                            Log In to Novik
                        </a>
                    </td>
                </tr>
            </table>
            
            <div class="divider"></div>
            
            <!-- Security Notice -->
            <div class="alert alert-info">
                <span class="alert-icon">🔒</span>
                <strong>Keep your account secure</strong><br>
                <span style="font-size: 14px;">
                    • Use a strong, unique password<br>
                    • Never share your password with anyone<br>
                    • Enable two-factor authentication when available
                </span>
            </div>
            
            <p class="text-muted" style="font-size: 14px;">
                <strong>Wasn't you?</strong><br>
                If you didn't make this change, please contact our support team immediately at 
                <a href="mailto:security@novik.ai" style="color: #667eea;">security@novik.ai</a>
            </p>
        </td>
    </tr>"""
    
    preheader = "Your Novik password has been successfully reset"
    return get_email_base_html(content, preheader, "security")


def get_password_reset_confirmation_text(user_name: str) -> str:
    """Plain text version for password reset confirmation."""
    return f"""Password Successfully Reset

Hi {user_name or 'there'},

Your Novik account password has been successfully reset. You can now log in with your new password.

Log in at: https://novik.ai/login

Keep your account secure:
• Use a strong, unique password
• Never share your password with anyone
• Enable two-factor authentication when available

Wasn't you?
If you didn't make this change, please contact our support team immediately at security@novik.ai

---
© 2024 Novik AI
Security: security@novik.ai
Website: https://novik.ai"""

# Mailgun-specific helper functions

def add_mailgun_variables(email_message, variables: dict):
    """
    Add Mailgun template variables for personalization.
    
    Args:
        email_message: Django EmailMessage object
        variables: Dictionary of template variables
    """
    if variables:
        email_message.extra_headers = email_message.extra_headers or {}
        email_message.extra_headers['X-Mailgun-Variables'] = str(variables)


def add_mailgun_tags(email_message, tags: list):
    """
    Add Mailgun tags for categorization and analytics.
    
    Args:
        email_message: Django EmailMessage object
        tags: List of tags (max 3 tags)
    """
    if tags:
        email_message.extra_headers = email_message.extra_headers or {}
        # Mailgun allows max 3 tags
        email_message.extra_headers['X-Mailgun-Tag'] = tags[:3]


def add_mailgun_tracking(email_message, track_opens: bool = True, track_clicks: bool = True):
    """
    Configure Mailgun tracking options.
    
    Args:
        email_message: Django EmailMessage object
        track_opens: Enable open tracking
        track_clicks: Enable click tracking
    """
    email_message.extra_headers = email_message.extra_headers or {}
    email_message.extra_headers['X-Mailgun-Track'] = 'yes' if (track_opens or track_clicks) else 'no'
    email_message.extra_headers['X-Mailgun-Track-Opens'] = 'yes' if track_opens else 'no'
    email_message.extra_headers['X-Mailgun-Track-Clicks'] = 'yes' if track_clicks else 'no'


def add_mailgun_delivery_time(email_message, delivery_time: str):
    """
    Schedule email delivery for a specific time.
    
    Args:
        email_message: Django EmailMessage object
        delivery_time: RFC 2822 formatted date string
    """
    email_message.extra_headers = email_message.extra_headers or {}
    email_message.extra_headers['X-Mailgun-Deliver-By'] = delivery_time


def add_mailgun_campaign(email_message, campaign_id: str):
    """
    Add email to a Mailgun campaign for analytics.
    
    Args:
        email_message: Django EmailMessage object
        campaign_id: Campaign identifier
    """
    email_message.extra_headers = email_message.extra_headers or {}
    email_message.extra_headers['X-Mailgun-Campaign-Id'] = campaign_id
