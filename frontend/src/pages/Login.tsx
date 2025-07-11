import React, { useState } from 'react';
import {
  Container,
  Box,
  Typography,
  TextField,
  FormControlLabel,
  Checkbox,
  Button,
  Link,
  CircularProgress,
} from '@mui/material';
import { Person as PersonIcon, Lock as LockIcon } from '@mui/icons-material';
import { InputAdornment } from '@mui/material';
import { GoogleOAuthProvider, GoogleLogin } from '@react-oauth/google';
import Header from '../components/Common/Header';
import { useAuthStore } from '../stores/auth';
import { authService } from '../services';
import FrontImage from '../assets/Front Image.png';

function Login() {
  const [formData, setFormData] = useState({ email: '', password: '' });
  const [isLogging, setIsLogging] = useState(false);
  const [error, setError] = useState<string | null>(null);

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
    } catch (err: any) {
      setError(err.response?.data?.error || 'Login failed');
    } finally {
      setIsLogging(false);
    }
  };

  const handleSuccess = async (res: any) => {
    setError(null);
    setIsLogging(true);
    try {
      const { accessToken, user } = await authService.loginWithGoogle(res.credential);
      useAuthStore.setState({ accessToken, isAuthorized: true, user });
    } catch (err: any) {
      setError(err.response?.data?.error || 'Google login failed');
    } finally {
      setIsLogging(false);
    }
  };

  const handleError = () => {
    setError('Google login failed');
  };

  return (
    <Box
      sx={{
        minHeight: '100vh',
        backgroundImage: `url(${FrontImage})`,
        backgroundSize: 'cover',
        backgroundPosition: 'center',
        backgroundRepeat: 'no-repeat',
      }}
    >
      <Header />
      <Container maxWidth="sm" sx={{ pt: 8 }}>
        <Box
          component="form"
          onSubmit={handleSubmit}
          sx={{ mt: 4, p: 3, bgcolor: 'background.paper', borderRadius: 2, boxShadow: 3 }}
        >
          <Typography variant="h6" gutterBottom>
            Login
          </Typography>
          <Typography variant="body2" color="textSecondary" gutterBottom>
            Login to use the NOVIK Dental Assistant
          </Typography>

          <TextField
            fullWidth
            margin="normal"
            name="email"
            label="Username or Email Address"
            value={formData.email}
            onChange={handleChange}
            required
            slotProps={{
              input: {
                startAdornment: (
                  <InputAdornment position="start">
                    <PersonIcon />
                  </InputAdornment>
                ),
              },
            }}
          />

          <TextField
            fullWidth
            margin="normal"
            name="password"
            type="password"
            label="Password"
            value={formData.password}
            onChange={handleChange}
            required
            slotProps={{
              input: {
                startAdornment: (
                  <InputAdornment position="start">
                    <LockIcon />
                  </InputAdornment>
                ),
              },
            }}
          />

          <FormControlLabel control={<Checkbox />} label="Remember Me" sx={{ mt: 1 }} />

          {error && (
            <Typography variant="body2" color="error" sx={{ mt: 1 }}>
              {error}
            </Typography>
          )}

          <Button
            type="submit"
            fullWidth
            variant="contained"
            disabled={isLogging}
            sx={{ mt: 2, position: 'relative' }}
            style={{ backgroundColor: '#F97316' }}
          >
            {isLogging ? <CircularProgress size={24} sx={{ color: 'white' }} /> : 'Log In'}
          </Button>

          <Box sx={{ mt: 3, textAlign: 'center' }}>
            <GoogleOAuthProvider clientId="415749549321-2g2mhh6ugbk8fhjfdcd4jo7sk00dfa8v.apps.googleusercontent.com">
              <GoogleLogin onSuccess={handleSuccess} onError={handleError} />
            </GoogleOAuthProvider>
          </Box>

          <Box sx={{ mt: 2, display: 'flex', justifyContent: 'space-between' }}>
            <Link href="/register" variant="body2" color="warning">
              Sign Up
            </Link>
            <Link href="/forgot-password" variant="body2" color="warning">
              Forgot Your Password?
            </Link>
          </Box>
        </Box>
      </Container>
    </Box>
  );
}

export default Login;
