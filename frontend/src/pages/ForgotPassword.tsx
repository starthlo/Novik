import { useState } from 'react';
import { Link as RouterLink, useNavigate } from 'react-router-dom';
import {
  Box,
  Button,
  Container,
  TextField,
  Typography,
  Alert,
  Link,
  CircularProgress,
  styled,
} from '@mui/material';
import { ArrowBack } from '@mui/icons-material';
import apiClient from '../lib/apiClient';
import NovikLogo from '../assets/novik-logo.png';
import { novikTheme } from '../styles/theme';

const PageContainer = styled(Box)({
  minHeight: '100vh',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '96px 16px 40px',
  backgroundColor: '#f7f7f8',
  fontFamily: novikTheme.typography.fontFamily,
});

const ForgotPasswordCard = styled(Box)({
  width: '100%',
  maxWidth: '420px',
  backgroundColor: '#ffffff',
  border: `1px solid ${novikTheme.colors.border}`,
  borderRadius: '16px',
  boxShadow: '0 20px 50px rgba(0,0,0,0.08)',
  padding: '24px',
});

const LogoWrapper = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  gap: '0.4rem',
  margin: '6px 0 10px',
});

const Logo = styled('img')({
  width: '80px',
  height: '80px',
  cursor: 'pointer',
});

const Title = styled(Typography)({
  textAlign: 'center',
  fontSize: '1.2rem',
  fontWeight: 600,
  margin: '0 0 0.25rem',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
});

const Subtitle = styled(Typography)({
  color: novikTheme.colors.textMuted,
  textAlign: 'center',
  marginBottom: '18px',
  fontSize: '0.9rem',
  fontFamily: novikTheme.typography.fontFamily,
});

const StyledTextField = styled(TextField)({
  marginBottom: '14px',
  '& .MuiInputBase-root': {
    borderRadius: '10px',
    fontSize: '0.95rem',
    fontFamily: novikTheme.typography.fontFamily,
    backgroundColor: '#ffffff',
  },
  '& .MuiInputLabel-root': {
    fontFamily: novikTheme.typography.fontFamily,
    fontWeight: 500,
    fontSize: '0.9rem',
  },
  '& .MuiOutlinedInput-root': {
    '& fieldset': {
      borderColor: novikTheme.colors.border,
    },
    '&:hover fieldset': {
      borderColor: novikTheme.colors.primary,
    },
    '&.Mui-focused fieldset': {
      borderColor: novikTheme.colors.primary,
    },
  },
});

const PrimaryButton = styled(Button)({
  width: '100%',
  padding: '0.9rem 1rem',
  borderRadius: '10px',
  fontWeight: 700,
  fontSize: '0.95rem',
  textTransform: 'none',
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
  '&:disabled': {
    backgroundColor: '#cccccc',
  },
});

const OutlinedButton = styled(Button)({
  width: '100%',
  padding: '0.9rem 1rem',
  borderRadius: '10px',
  fontWeight: 700,
  fontSize: '0.95rem',
  textTransform: 'none',
  border: `1px solid ${novikTheme.colors.primary}`,
  color: novikTheme.colors.primary,
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: 'rgba(30, 64, 175, 0.04)',
    border: `1px solid ${novikTheme.colors.primaryDark}`,
  },
});

const LinkText = styled(Link)<any>({
  color: novikTheme.colors.primary,
  fontSize: '0.9rem',
  textDecoration: 'none',
  fontWeight: 500,
  fontFamily: novikTheme.typography.fontFamily,
  display: 'inline-flex',
  alignItems: 'center',
  '&:hover': {
    textDecoration: 'underline',
  },
});

const SuccessCard = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  textAlign: 'center',
});

export default function ForgotPasswordPage() {
  const navigate = useNavigate();
  const [email, setEmail] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setLoading(true);

    try {
      const response = await apiClient.post('/auth/forgot-password/', {
        email,
      });

      if (response.data.success) {
        setSuccess(true);
      } else {
        setError(response.data.message || 'An error occurred. Please try again.');
      }
    } catch (err: any) {
      setError(err.response?.data?.message || 'Failed to send reset email. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  if (success) {
    return (
      <PageContainer>
        <Container maxWidth="sm">
          <ForgotPasswordCard>
            <LogoWrapper>
              <Link component={RouterLink} to="/">
                <Logo src={NovikLogo} alt="Novik Logo" />
              </Link>
            </LogoWrapper>

            <SuccessCard>
              <Title variant="h1">Check Your Email</Title>
              <Subtitle>
                We've sent a password reset link to <strong>{email}</strong>
              </Subtitle>

              <Alert severity="success" sx={{ mb: 2, width: '100%' }}>
                Please check your email and follow the instructions to reset your password.
              </Alert>

              <Typography
                variant="body2"
                sx={{
                  color: novikTheme.colors.textMuted,
                  mb: 3,
                  fontSize: '0.85rem',
                }}
              >
                The link will expire in 24 hours. If you don't see the email, please check your spam
                folder.
              </Typography>

              <PrimaryButton onClick={() => navigate('/login')}>Back to Login</PrimaryButton>

              <Box sx={{ mt: 2 }}>
                <OutlinedButton
                  onClick={() => {
                    setSuccess(false);
                    setEmail('');
                  }}
                >
                  Send Another Email
                </OutlinedButton>
              </Box>
            </SuccessCard>
          </ForgotPasswordCard>
        </Container>
      </PageContainer>
    );
  }

  return (
    <PageContainer>
      <Container maxWidth="sm">
        <ForgotPasswordCard>
          <LogoWrapper>
            <Link component={RouterLink} to="/">
              <Logo src={NovikLogo} alt="Novik Logo" />
            </Link>
          </LogoWrapper>

          <Title variant="h1">Forgot Password?</Title>
          <Subtitle>
            Enter your email address and we'll send you a link to reset your password
          </Subtitle>

          {error && (
            <Alert
              severity="error"
              sx={{
                mb: 2,
                fontSize: '0.85rem',
                borderRadius: '8px',
              }}
            >
              {error}
            </Alert>
          )}

          <Box component="form" onSubmit={handleSubmit}>
            <StyledTextField
              fullWidth
              type="email"
              label="Email Address"
              value={email}
              onChange={e => setEmail(e.target.value)}
              required
              disabled={loading}
              autoFocus
            />

            <PrimaryButton type="submit" disabled={loading || !email}>
              {loading ? (
                <>
                  <CircularProgress size={20} sx={{ mr: 1, color: 'white' }} />
                  Sending...
                </>
              ) : (
                'Send Reset Link'
              )}
            </PrimaryButton>

            <Box sx={{ mt: 2, textAlign: 'center' }}>
              <LinkText component={RouterLink} to="/login">
                <ArrowBack sx={{ mr: 0.5, fontSize: 18 }} />
                Back to Login
              </LinkText>
            </Box>
          </Box>
        </ForgotPasswordCard>
      </Container>
    </PageContainer>
  );
}
