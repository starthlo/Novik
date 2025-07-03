import { useEffect, useState } from 'react';
import { Box, Button, Typography, Link, IconButton } from '@mui/material';
import CloseIcon from '@mui/icons-material/Close';

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
    <Box
      sx={{
        position: 'fixed',
        bottom: theme => theme.spacing(6),
        right: theme => theme.spacing(6),
        zIndex: theme => theme.zIndex.modal,
      }}
    >
      <Box
        sx={{
          position: 'relative',
          backgroundColor: '#ffffff',
          border: '1px solid #E5E7EB',
          borderRadius: 2,
          boxShadow: 24,
          maxWidth: 600,
          width: '100%',
          p: 4,
        }}
      >
        <IconButton
          onClick={() => setVisible(false)}
          sx={{
            position: 'absolute',
            top: 8,
            right: 8,
            color: '#6B7280',
            '&:hover': { color: '#000000' },
          }}
          size="small"
          aria-label="Close"
        >
          <CloseIcon fontSize="small" />
        </IconButton>

        <Typography variant="h6" sx={{ mb: 1, fontWeight: 600 }}>
          Manage cookie consent
        </Typography>

        <Typography variant="body2" sx={{ mb: 3, color: '#374151' }}>
          We use cookies to enhance your experience. You can accept all, only essential cookies, or
          decline.
        </Typography>

        <Box
          sx={{ display: 'flex', flexWrap: 'wrap', gap: 1, justifyContent: 'space-between', mb: 2 }}
        >
          <Button
            variant="contained"
            onClick={() => handleConsent('all')}
            sx={{
              backgroundColor: '#F97316',
              color: '#ffffff',
              '&:hover': { backgroundColor: '#EA580C' },
              px: 2,
              py: 1,
              fontWeight: 500,
            }}
          >
            Accept
          </Button>
          <Button
            variant="outlined"
            onClick={() => handleConsent('functional')}
            sx={{
              borderColor: '#E5E7EB',
              color: '#1F2937',
              '&:hover': {
                backgroundColor: '#EA580C',
                color: '#ffffff',
                borderColor: '#EA580C',
              },
              px: 2,
              py: 1,
              fontWeight: 500,
            }}
          >
            Accept only functional cookies
          </Button>
          <Button
            variant="outlined"
            onClick={() => handleConsent('decline')}
            sx={{
              borderColor: '#E5E7EB',
              color: '#1F2937',
              '&:hover': {
                backgroundColor: '#EA580C',
                color: '#ffffff',
                borderColor: '#EA580C',
              },
              px: 2,
              py: 1,
              fontWeight: 500,
            }}
          >
            Decline
          </Button>
        </Box>

        <Box
          sx={{ display: 'flex', justifyContent: 'center', gap: 3, mt: 2, fontSize: '0.875rem' }}
        >
          <Link
            href="/legal#terms"
            target="_blank"
            rel="noopener noreferrer"
            underline="hover"
            sx={{ color: '#F97316' }}
          >
            Terms of service
          </Link>
          <Link
            href="/legal#privacy"
            target="_blank"
            rel="noopener noreferrer"
            underline="hover"
            sx={{ color: '#F97316' }}
          >
            Privacy Policy
          </Link>
          <Link
            href="/legal#cookies"
            target="_blank"
            rel="noopener noreferrer"
            underline="hover"
            sx={{ color: '#F97316' }}
          >
            Cookie Policy
          </Link>
        </Box>
      </Box>
    </Box>
  );
};

export default CookieConsent;
