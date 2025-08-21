import { useEffect, useState } from 'react';
import { Box, Button, Typography, Link, IconButton, styled } from '@mui/material';
import CloseIcon from '@mui/icons-material/Close';
import { novikTheme } from '../../styles/theme';

const ConsentContainer = styled(Box)({
  position: 'fixed',
  bottom: '24px',
  right: '24px',
  zIndex: 9999,
  '@media (max-width: 600px)': {
    bottom: '16px',
    right: '16px',
    left: '16px',
  },
});

const ConsentCard = styled(Box)({
  position: 'relative',
  backgroundColor: novikTheme.colors.background,
  border: `1px solid ${novikTheme.colors.border}`,
  borderRadius: novikTheme.borderRadius.medium,
  boxShadow: novikTheme.shadows.large,
  maxWidth: '600px',
  width: '100%',
  padding: '32px',
  fontFamily: novikTheme.typography.fontFamily,
  '@media (max-width: 600px)': {
    padding: '24px',
  },
});

const CloseButton = styled(IconButton)({
  position: 'absolute',
  top: '8px',
  right: '8px',
  color: novikTheme.colors.textMuted,
  '&:hover': {
    color: novikTheme.colors.text,
    backgroundColor: 'rgba(0, 0, 0, 0.04)',
  },
});

const ConsentTitle = styled(Typography)({
  marginBottom: '8px',
  fontWeight: 600,
  fontFamily: novikTheme.typography.fontFamily,
  color: novikTheme.colors.text,
});

const ConsentDescription = styled(Typography)({
  marginBottom: '24px',
  color: novikTheme.colors.textMuted,
  fontFamily: novikTheme.typography.fontFamily,
  lineHeight: 1.6,
});

const AcceptButton = styled(Button)({
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  fontFamily: novikTheme.typography.fontFamily,
  textTransform: 'none',
  fontWeight: 500,
  padding: '8px 20px',
  borderRadius: novikTheme.borderRadius.small,
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const SecondaryButton = styled(Button)({
  borderColor: novikTheme.colors.border,
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
  textTransform: 'none',
  fontWeight: 500,
  padding: '8px 20px',
  borderRadius: novikTheme.borderRadius.small,
  '&:hover': {
    backgroundColor: 'rgba(136, 169, 78, 0.08)',
    borderColor: novikTheme.colors.primary,
    color: novikTheme.colors.primary,
  },
});

const PolicyLink = styled(Link)({
  color: novikTheme.colors.primary,
  textDecoration: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  fontSize: '0.875rem',
  '&:hover': {
    textDecoration: 'underline',
  },
});

const CookieConsent = () => {
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const consent = localStorage.getItem('cookieConsent');
    if (!consent) setVisible(true);
  }, []);

  const handleConsent = (value: string) => {
    localStorage.setItem('cookieConsent', value);
    setVisible(false);
  };

  if (!visible) return null;

  return (
    <ConsentContainer>
      <ConsentCard>
        <CloseButton
          onClick={() => setVisible(false)}
          size="small"
          aria-label="Close cookie consent"
        >
          <CloseIcon fontSize="small" />
        </CloseButton>

        <ConsentTitle variant="h6">Cookie Consent</ConsentTitle>

        <ConsentDescription variant="body2">
          We use cookies to enhance your browsing experience and analyze our traffic. By clicking
          "Accept All", you consent to our use of cookies. You can also choose to accept only
          essential cookies or decline optional cookies.
        </ConsentDescription>

        <Box
          sx={{
            display: 'flex',
            flexWrap: 'wrap',
            gap: 1.5,
            mb: 2,
            '@media (max-width: 600px)': {
              flexDirection: 'column',
            },
          }}
        >
          <AcceptButton variant="contained" onClick={() => handleConsent('all')}>
            Accept All
          </AcceptButton>
          <SecondaryButton variant="outlined" onClick={() => handleConsent('functional')}>
            Essential Only
          </SecondaryButton>
          <SecondaryButton variant="outlined" onClick={() => handleConsent('decline')}>
            Decline
          </SecondaryButton>
        </Box>

        <Box
          sx={{
            display: 'flex',
            justifyContent: 'center',
            gap: 3,
            mt: 2,
            flexWrap: 'wrap',
          }}
        >
          <PolicyLink href="/terms-of-use" underline="hover" target="_blank">
            Terms of Use
          </PolicyLink>
          <PolicyLink href="/privacy-policy" underline="hover" target="_blank">
            Privacy Policy
          </PolicyLink>
          <PolicyLink href="/cookie-policy" underline="hover" target="_blank">
            Cookie Policy
          </PolicyLink>
        </Box>
      </ConsentCard>
    </ConsentContainer>
  );
};

export default CookieConsent;
