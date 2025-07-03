import { ReactNode } from 'react';
import { Box, Container, Grid, Typography, Link, IconButton } from '@mui/material';
import {
  FaEnvelope,
  FaPhone,
  FaWhatsapp,
  FaFacebookF,
  FaInstagram,
  FaLinkedinIn,
} from 'react-icons/fa';
import { Link as RouterLink } from 'react-router-dom';

const Footer = () => {
  const contactLinks: { icon: ReactNode; href: string; label: string }[] = [
    { icon: <FaEnvelope />, href: 'mailto:info@novik.ai', label: 'Email' },
    { icon: <FaPhone />, href: 'tel:+34690957910', label: 'Call' },
    { icon: <FaWhatsapp />, href: 'https://wa.me/34690957910', label: 'WhatsApp' },
  ];

  const socialLinks: { icon: ReactNode; href: string; label: string }[] = [
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
    <Box
      component="footer"
      sx={{
        backgroundColor: '#ffffff',
        color: '#6B7280',
        borderTop: '1px solid #E5E7EB',
        py: 4,
        fontSize: '0.875rem',
      }}
    >
      <Container maxWidth="lg">
        <Grid container spacing={4} alignItems="flex-start">
          <Grid size={{ xs: 12, md: 4 }}>
            <Typography
              variant="h5"
              component="div"
              gutterBottom
              sx={{ fontWeight: 600, color: '#6B7280' }}
            >
              Novik
              <Box component="span" sx={{ color: '#EA580C' }}>
                .
              </Box>
            </Typography>
            <Typography variant="body2" gutterBottom>
              © 2024 Novik. All rights reserved.
            </Typography>
            <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1, mt: 1 }}>
              <Link
                component={RouterLink}
                to="/legal#terms"
                underline="hover"
                sx={{ color: '#F97316' }}
              >
                Terms of Service
              </Link>
              <Typography component="span">|</Typography>
              <Link
                component={RouterLink}
                to="/legal#privacy"
                underline="hover"
                sx={{ color: '#F97316' }}
              >
                Privacy Policy
              </Link>
              <Typography component="span">|</Typography>
              <Link
                component={RouterLink}
                to="/legal#cookies"
                underline="hover"
                sx={{ color: '#F97316' }}
              >
                Cookie Policy
              </Link>
            </Box>
          </Grid>

          <Grid size={{ xs: 12, md: 4 }}>
            <Typography variant="subtitle1" gutterBottom sx={{ color: '#6B7280', fontWeight: 600 }}>
              Contact Us
            </Typography>
            <Box>
              {contactLinks.map(({ icon, href, label }, idx) => (
                <IconButton
                  key={idx}
                  component="a"
                  href={href}
                  target={href.startsWith('http') ? '_blank' : undefined}
                  rel={href.startsWith('http') ? 'noopener noreferrer' : undefined}
                  aria-label={label}
                  sx={{
                    color: '#F97316',
                    transition: 'transform 0.2s',
                    '&:hover': { transform: 'scale(1.1)' },
                  }}
                >
                  {icon}
                </IconButton>
              ))}
            </Box>
          </Grid>

          <Grid size={{ xs: 12, md: 4 }}>
            <Typography variant="subtitle1" gutterBottom sx={{ color: '#6B7280', fontWeight: 600 }}>
              Social Links
            </Typography>
            <Box>
              {socialLinks.map(({ icon, href, label }, idx) => (
                <IconButton
                  key={idx}
                  component="a"
                  href={href}
                  target="_blank"
                  rel="noopener noreferrer"
                  aria-label={label}
                  sx={{
                    color: '#F97316',
                    transition: 'transform 0.2s',
                    '&:hover': { transform: 'scale(1.1)' },
                  }}
                >
                  {icon}
                </IconButton>
              ))}
            </Box>
          </Grid>
        </Grid>

        <Box sx={{ mt: 4 }}>
          <Typography variant="body2" sx={{ color: '#6B7280' }}>
            Novik is an experimental technology demonstrator. Novik does not provide medical advice,
            diagnosis or treatment. User questions and other inputs on Novik are not covered by
            HIPAA. It is the responsibility of the user to ensure questions do not contain protected
            health information (PHI) or any information that violates the privacy of any person.
          </Typography>
        </Box>
      </Container>
    </Box>
  );
};

export default Footer;
