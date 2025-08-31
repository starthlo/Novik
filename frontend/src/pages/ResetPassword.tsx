import { useState, useEffect } from 'react';
import { useNavigate, useSearchParams, Link as RouterLink } from 'react-router-dom';
import {
  Box,
  Button,
  Container,
  TextField,
  Typography,
  Alert,
  CircularProgress,
  InputAdornment,
  IconButton,
  styled,
  Link,
} from '@mui/material';
import { Visibility, VisibilityOff, CheckCircle } from '@mui/icons-material';
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

const ResetCard = styled(Box)({
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

const SuccessCard = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  textAlign: 'center',
});

const HintText = styled(Typography)({
  fontSize: '0.75rem',
  color: novikTheme.colors.textMuted,
  marginTop: '-8px',
  marginBottom: '14px',
  fontFamily: novikTheme.typography.fontFamily,
});

export default function ResetPasswordPage() {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();

  const [password, setPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [showPassword, setShowPassword] = useState(false);
  const [showConfirmPassword, setShowConfirmPassword] = useState(false);
  const [loading, setLoading] = useState(false);
  const [validating, setValidating] = useState(true);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState(false);
  const [tokenValid, setTokenValid] = useState(false);
  const [userEmail, setUserEmail] = useState('');

  const uid = searchParams.get('uid');
  const token = searchParams.get('token');

  useEffect(() => {
    // Validate the reset token when component mounts
    const validateToken = async () => {
      if (!uid || !token) {
        setError('Invalid reset link. Please request a new password reset.');
        setValidating(false);
        return;
      }

      try {
        const response = await apiClient.post('/auth/validate-reset-token/', {
          uid,
          token,
        });

        if (response.data.valid) {
          setTokenValid(true);
          setUserEmail(response.data.email || '');
        } else {
          setError(response.data.message || 'Invalid or expired reset link.');
        }
      } catch (err: any) {
        setError(
          err.response?.data?.message ||
          'Invalid or expired reset link. Please request a new password reset.'
        );
      } finally {
        setValidating(false);
      }
    };

    validateToken();
  }, [uid, token]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');

    // Validate passwords match
    if (password !== confirmPassword) {
      setError('Passwords do not match.');
      return;
    }

    // Basic password validation
    if (password.length < 8) {
      setError('Password must be at least 8 characters long.');
      return;
    }

    setLoading(true);

    try {
      const response = await apiClient.post('/auth/reset-password/', {
        uid,
        token,
        password,
        confirm_password: confirmPassword,
      });

      if (response.data.success) {
        setSuccess(true);
      } else {
        setError(response.data.message || 'Failed to reset password.');
      }
    } catch (err: any) {
      if (err.response?.data?.errors) {
        // Handle validation errors
        const errors = err.response.data.errors;
        const errorMessages = Object.values(errors).flat().join(' ');
        setError(errorMessages);
      } else {
        setError(err.response?.data?.message || 'Failed to reset password. Please try again.');
      }
    } finally {
      setLoading(false);
    }
  };

  if (validating) {
    return (
      <PageContainer>
        <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 2 }}>
          <CircularProgress />
          <Typography sx={{ color: novikTheme.colors.textMuted }}>
            Validating reset link...
          </Typography>
        </Box>
      </PageContainer>
    );
  }

  if (success) {
    return (
      <PageContainer>
        <Container maxWidth="sm">
          <ResetCard>
            <SuccessCard>
              <CheckCircle sx={{ fontSize: 60, color: 'success.main', mb: 2 }} />
              <Title variant="h1">Password Reset Successful!</Title>
              <Subtitle>
                Your password has been successfully reset. You can now log in with your new
                password.
              </Subtitle>

              <PrimaryButton onClick={() => navigate('/login')}>Go to Login</PrimaryButton>
            </SuccessCard>
          </ResetCard>
        </Container>
      </PageContainer>
    );
  }

  if (!tokenValid) {
    return (
      <PageContainer>
        <Container maxWidth="sm">
          <ResetCard>
            <LogoWrapper>
              <Link component={RouterLink} to="/">
                <Logo src={NovikLogo} alt="Novik Logo" />
              </Link>
            </LogoWrapper>

            <Title variant="h1">Invalid Reset Link</Title>

            <Alert
              severity="error"
              sx={{
                mb: 3,
                fontSize: '0.85rem',
                borderRadius: '8px',
              }}
            >
              {error}
            </Alert>

            <PrimaryButton onClick={() => navigate('/forgot-password')}>
              Request New Password Reset
            </PrimaryButton>
          </ResetCard>
        </Container>
      </PageContainer>
    );
  }

  return (
    <PageContainer>
      <Container maxWidth="sm">
        <ResetCard>
          <LogoWrapper>
            <Link component={RouterLink} to="/">
              <Logo src={NovikLogo} alt="Novik Logo" />
            </Link>
          </LogoWrapper>

          <Title variant="h1">Reset Your Password</Title>
          {userEmail && (
            <Subtitle>
              Resetting password for: <strong>{userEmail}</strong>
            </Subtitle>
          )}

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
              type={showPassword ? 'text' : 'password'}
              label="New Password"
              value={password}
              onChange={e => setPassword(e.target.value)}
              required
              disabled={loading}
              slotProps={{
                input: {
                  endAdornment: (
                    <InputAdornment position="end">
                      <IconButton
                        onClick={() => setShowPassword(!showPassword)}
                        edge="end"
                        size="small"
                      >
                        {showPassword ? <VisibilityOff /> : <Visibility />}
                      </IconButton>
                    </InputAdornment>
                  ),
                },
              }}
            />

            <StyledTextField
              fullWidth
              type={showConfirmPassword ? 'text' : 'password'}
              label="Confirm New Password"
              value={confirmPassword}
              onChange={e => setConfirmPassword(e.target.value)}
              required
              disabled={loading}
              slotProps={{
                input: {
                  endAdornment: (
                    <InputAdornment position="end">
                      <IconButton
                        onClick={() => setShowConfirmPassword(!showConfirmPassword)}
                        edge="end"
                        size="small"
                      >
                        {showConfirmPassword ? <VisibilityOff /> : <Visibility />}
                      </IconButton>
                    </InputAdornment>
                  ),
                },
              }}
            />

            <HintText>Password must be at least 8 characters long</HintText>

            <PrimaryButton type="submit" disabled={loading || !password || !confirmPassword}>
              {loading ? (
                <>
                  <CircularProgress size={20} sx={{ mr: 1, color: 'white' }} />
                  Resetting Password...
                </>
              ) : (
                'Reset Password'
              )}
            </PrimaryButton>
          </Box>
        </ResetCard>
      </Container>
    </PageContainer>
  );
}
