import React, { useState } from 'react';
import {
  Container,
  Box,
  Typography,
  TextField,
  Button,
  Link,
  CircularProgress,
  styled,
} from '@mui/material';
import { GoogleOAuthProvider, GoogleLogin } from '@react-oauth/google';
import { Link as RouterLink, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../stores/auth';
import { authService } from '../services';
import { novikTheme } from '../styles/theme';
import NovikLogo from '../assets/novik-logo.png';

const PageContainer = styled(Box)({
  minHeight: '100vh',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '96px 16px 40px',
  backgroundColor: '#f7f7f8',
  fontFamily: novikTheme.typography.fontFamily,
});

const LoginCard = styled(Box)({
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

const DividerText = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  gap: '12px',
  color: novikTheme.colors.textMuted,
  fontSize: '0.85rem',
  margin: '14px 0',
  fontFamily: novikTheme.typography.fontFamily,
  '&::before, &::after': {
    content: '""',
    flex: 1,
    height: '1px',
    backgroundColor: novikTheme.colors.border,
  },
});

const LinkText = styled(Link)<any>({
  color: novikTheme.colors.primary,
  fontSize: '0.9rem',
  textDecoration: 'none',
  fontWeight: 500,
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    textDecoration: 'underline',
  },
});

function Login() {
  const navigate = useNavigate();
  const [formData, setFormData] = useState({ email: '', password: '' });
  const [isLogging, setIsLogging] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isLogin, setIsLogin] = useState(true); // Toggle between login and initial registration screen

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setFormData(prev => ({ ...prev, [e.target.name]: e.target.value }));
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setIsLogging(true);
    try {
      const { accessToken, user } = await authService.login(formData.email, formData.password);
      useAuthStore.setState({ accessToken, isAuthorized: true, user });
      navigate('/dashboard');
    } catch (err: any) {
      setError(err.response?.data?.error || 'Login failed. Please check your credentials.');
    } finally {
      setIsLogging(false);
    }
  };

  const handleEmailContinue = () => {
    if (formData.email) {
      setIsLogin(true);
    }
  };

  const handleGoogleSuccess = async (res: any) => {
    setError(null);
    setIsLogging(true);
    try {
      const { accessToken, user } = await authService.loginWithGoogle(res.credential);
      useAuthStore.setState({ accessToken, isAuthorized: true, user });
      navigate('/dashboard');
    } catch (err: any) {
      setError(err.response?.data?.error || 'Google login failed');
    } finally {
      setIsLogging(false);
    }
  };

  const handleGoogleError = () => {
    setError('Google login failed');
  };

  return (
    <PageContainer>
      <Container maxWidth="sm">
        <LoginCard>
          <LogoWrapper>
            <Link component={RouterLink} to="/">
              <Logo src={NovikLogo} alt="Novik Logo" />
            </Link>
          </LogoWrapper>

          {!isLogin ? (
            <>
              <Title variant="h1">Welcome to Novik</Title>
              <Subtitle>Use your professional email to get started</Subtitle>

              <StyledTextField
                fullWidth
                name="email"
                label="Username or Email Address"
                value={formData.email}
                onChange={handleChange}
                required
                autoFocus
              />

              <PrimaryButton onClick={handleEmailContinue} disabled={!formData.email}>
                Continue
              </PrimaryButton>

              <DividerText>OR</DividerText>

              <GoogleOAuthProvider clientId="415749549321-2g2mhh6ugbk8fhjfdcd4jo7sk00dfa8v.apps.googleusercontent.com">
                <Box sx={{ width: '100%' }}>
                  <GoogleLogin
                    onSuccess={handleGoogleSuccess}
                    onError={handleGoogleError}
                    width="100%"
                    theme="outline"
                    size="large"
                    text="continue_with"
                  />
                </Box>
              </GoogleOAuthProvider>

              <Box sx={{ mt: 2, textAlign: 'center' }}>
                <Typography variant="body2" sx={{ color: novikTheme.colors.textMuted }}>
                  Already have an account?{' '}
                  <LinkText component="button" onClick={() => setIsLogin(true)}>
                    Sign in
                  </LinkText>
                </Typography>
              </Box>
            </>
          ) : (
            <>
              <Title variant="h1">Sign in to Novik</Title>
              <Subtitle>Enter your credentials to access the dental assistant</Subtitle>

              <Box component="form" onSubmit={handleSubmit}>
                <StyledTextField
                  fullWidth
                  name="email"
                  label="Username or Email Address"
                  value={formData.email}
                  onChange={handleChange}
                  required
                />

                <StyledTextField
                  fullWidth
                  name="password"
                  type="password"
                  label="Password"
                  value={formData.password}
                  onChange={handleChange}
                  required
                  autoFocus={!!formData.email}
                />

                {error && (
                  <Typography
                    variant="body2"
                    sx={{
                      color: novikTheme.colors.danger || '#d14343',
                      mb: 2,
                      fontSize: '0.85rem',
                      textAlign: 'center',
                    }}
                  >
                    {error}
                  </Typography>
                )}

                <PrimaryButton type="submit" disabled={isLogging}>
                  {isLogging ? <CircularProgress size={24} sx={{ color: 'white' }} /> : 'Sign In'}
                </PrimaryButton>
              </Box>

              <Box sx={{ mt: 2, display: 'flex', justifyContent: 'space-between' }}>
                <LinkText component={RouterLink} to="/register">
                  Create account
                </LinkText>
                <LinkText component={RouterLink} to="/forgot-password">
                  Forgot password?
                </LinkText>
              </Box>

              <DividerText>OR</DividerText>

              <GoogleOAuthProvider clientId="415749549321-2g2mhh6ugbk8fhjfdcd4jo7sk00dfa8v.apps.googleusercontent.com">
                <Box sx={{ width: '100%' }}>
                  <GoogleLogin
                    onSuccess={handleGoogleSuccess}
                    onError={handleGoogleError}
                    width="100%"
                    theme="outline"
                    size="large"
                    text="signin_with"
                  />
                </Box>
              </GoogleOAuthProvider>
            </>
          )}
        </LoginCard>
      </Container>
    </PageContainer>
  );
}

export default Login;
