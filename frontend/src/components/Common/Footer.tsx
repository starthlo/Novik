import { Box, Container, Typography, Link, styled, Grid } from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';
import { novikTheme } from '../../styles/theme';
import {
  FaFacebookF,
  FaInstagram,
  FaLinkedinIn,
} from 'react-icons/fa';

// Styled components matching HTML mockup design
const FooterContainer = styled(Box)<any>({
  backgroundColor: '#111111',
  color: '#cccccc',
  padding: '2rem 1rem',
  fontFamily: novikTheme.typography.fontFamily,
});

const FooterGrid = styled(Container)({
  maxWidth: '1200px',
  margin: '0 auto',
});

const FooterColumn = styled(Box)({
  flex: '1 1 200px',
  minWidth: '200px',
  marginBottom: '2rem',
});

const FooterHeading = styled(Typography)({
  color: '#ffffff',
  marginBottom: '1rem',
  fontSize: '1.1rem',
  fontWeight: 600,
  fontFamily: novikTheme.typography.fontFamily,
});

const FooterText = styled(Typography)({
  color: '#cccccc',
  fontSize: '0.9rem',
  lineHeight: 1.6,
  fontFamily: novikTheme.typography.fontFamily,
});

const FooterLink = styled(Link)<any>({
  color: '#cccccc',
  fontSize: '0.9rem',
  textDecoration: 'none',
  display: 'block',
  marginBottom: '0.6rem',
  fontFamily: novikTheme.typography.fontFamily,
  transition: 'color 0.2s ease',
  '&:hover': {
    color: '#ffffff',
  },
});

const FooterBottom = styled(Box)({
  textAlign: 'center',
  marginTop: '2rem',
  paddingTop: '2rem',
  borderTop: '1px solid #333333',
  fontSize: '0.8rem',
  color: '#777777',
});

const SocialIconsContainer = styled(Box)({
  display: 'flex',
  gap: '15px',
  marginTop: '8px',
});

const SocialIcon = styled('a')({
  color: '#cccccc',
  transition: 'color 0.2s ease, transform 0.2s ease',
  display: 'inline-flex',
  alignItems: 'center',
  justifyContent: 'center',
  '&:hover': {
    color: '#ffffff',
    transform: 'translateY(-2px)',
  },
  '& svg': {
    width: '26px',
    height: '26px',
  },
});

const Footer = () => {
  const navigationLinks = [
    { text: 'Home', to: '/' },
    { text: 'Why free', to: '/why-free' },
    { text: 'Partners', to: '/partners' },
    { text: 'API', to: '/api-novik' },
    { text: 'FAQs', to: '/faqs' },
    { text: 'Contact', to: '/contact' },
  ];

  const legalLinks = [
    { text: 'Legal Notice', to: '/legal#notice' },
    { text: 'Terms of Service', to: '/legal#terms' },
    { text: 'Privacy Policy', to: '/legal#privacy' },
    { text: 'Cookie Policy', to: '/legal#cookies' },
  ];

  const socialLinks = [
    {
      icon: <FaFacebookF />,
      href: 'https://www.facebook.com/profile.php?id=61567745501156',
      label: 'Facebook',
    },
    {
      icon: <FaInstagram />,
      href: 'https://www.instagram.com/dentalnovik/',
      label: 'Instagram',
    },
    {
      icon: <FaLinkedinIn />,
      href: 'https://www.linkedin.com/company/novik-ai',
      label: 'LinkedIn',
    },
  ];

  return (
    <FooterContainer component="footer">
      <FooterGrid>
        <Grid container spacing={3}>
          {/* Company Info Column */}
          <Grid size={{ xs: 12, md: 3 }}>
            <FooterColumn>
              <FooterHeading variant="h4">Novik</FooterHeading>
              <FooterText>
                AI decision support for dentistry. Making complex data simple, safe and actionable.
              </FooterText>
            </FooterColumn>
          </Grid>

          {/* Navigation Column */}
          <Grid size={{ xs: 12, sm: 6, md: 3 }}>
            <FooterColumn>
              <FooterHeading variant="h4">Navigation</FooterHeading>
              {navigationLinks.map((link, index) => (
                <FooterLink
                  key={index}
                  component={RouterLink}
                  to={link.to}
                >
                  {link.text}
                </FooterLink>
              ))}
            </FooterColumn>
          </Grid>

          {/* Legal Column */}
          <Grid size={{ xs: 12, sm: 6, md: 3 }}>
            <FooterColumn>
              <FooterHeading variant="h4">Legal</FooterHeading>
              {legalLinks.map((link, index) => (
                <FooterLink
                  key={index}
                  component={RouterLink}
                  to={link.to}
                >
                  {link.text}
                </FooterLink>
              ))}
            </FooterColumn>
          </Grid>

          {/* Social Links Column */}
          <Grid size={{ xs: 12, md: 3 }}>
            <FooterColumn>
              <FooterHeading variant="h4">Follow us</FooterHeading>
              <SocialIconsContainer>
                {socialLinks.map((social, index) => (
                  <SocialIcon
                    key={index}
                    href={social.href}
                    target="_blank"
                    rel="noopener noreferrer"
                    aria-label={social.label}
                  >
                    {social.icon}
                  </SocialIcon>
                ))}
              </SocialIconsContainer>
            </FooterColumn>
          </Grid>
        </Grid>

        {/* Footer Bottom */}
        <FooterBottom>
          <Typography variant="body2" sx={{ mb: 2, color: '#777777', fontSize: '0.8rem' }}>
            © 2024 Novik. All rights reserved.
          </Typography>
          <Typography variant="body2" sx={{ color: '#777777', fontSize: '0.75rem', lineHeight: 1.6 }}>
            Novik is an experimental technology demonstrator. Novik does not provide medical advice,
            diagnosis or treatment. User questions and other inputs on Novik are not covered by
            HIPAA. It is the responsibility of the user to ensure questions do not contain protected
            health information (PHI) or any information that violates the privacy of any person.
          </Typography>
        </FooterBottom>
      </FooterGrid>
    </FooterContainer>
  );
};

export default Footer;
