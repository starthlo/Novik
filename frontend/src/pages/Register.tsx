import React, { useState, useEffect } from 'react';
import {
  Container,
  Box,
  Typography,
  TextField,
  Button,
  Link,
  CircularProgress,
  Checkbox,
  FormControlLabel,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  styled,
  List,
  ListItem,
  Autocomplete,
} from '@mui/material';
import { Check, Close } from '@mui/icons-material';
import { GoogleOAuthProvider, GoogleLogin } from '@react-oauth/google';
import { Link as RouterLink, useNavigate } from 'react-router-dom';
import { Country, State, City } from 'country-state-city';
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

const RegisterCard = styled(Box)({
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
  transition: 'transform 0.2s ease',
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

const StyledFormControl = styled(FormControl)({
  marginBottom: '14px',
  '& .MuiInputBase-root': {
    borderRadius: '10px',
    fontSize: '0.95rem',
    fontFamily: novikTheme.typography.fontFamily,
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

const CheckboxGroup = styled(Box)({
  display: 'flex',
  alignItems: 'flex-start',
  gap: '6px',
  marginBottom: '14px',
  '& .MuiCheckbox-root': {
    padding: '0 6px 0 0',
  },
  '& .MuiFormControlLabel-label': {
    fontWeight: 400,
    fontSize: '0.8rem',
    lineHeight: 1.2,
    fontFamily: novikTheme.typography.fontFamily,
  },
});

const PasswordChecklist = styled(List)({
  border: `1px solid ${novikTheme.colors.border}`,
  borderRadius: '12px',
  padding: '10px 12px',
  marginTop: '10px',
  marginBottom: '14px',
  fontSize: '0.9rem',
});

const ChecklistItem = styled(ListItem)<{ valid: boolean }>(({ valid }) => ({
  listStyle: 'none',
  margin: '0.35rem 0',
  padding: 0,
  color: valid ? '#179b4d' : '#d14343',
  fontSize: '0.9rem',
  display: 'flex',
  alignItems: 'center',
  gap: '8px',
  fontFamily: novikTheme.typography.fontFamily,
}));

const LinkText = styled(Link)<any>({
  color: novikTheme.colors.primary,
  textDecoration: 'none',
  fontWeight: 500,
  '&:hover': {
    textDecoration: 'underline',
  },
});

function Register() {
  const navigate = useNavigate();
  const [step, setStep] = useState(1);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const [formData, setFormData] = useState({
    email: '',
    fullName: '',
    occupation: '',
    licenseId: '',
    password: '',
    confirmPassword: '',
    attestProfessional: false,
    agreeToTerms: false,
    dob: '',
    country: '',
    state: '',
    city: '',
  });

  // Location data states
  const [countries, setCountries] = useState<any[]>([]);
  const [states, setStates] = useState<any[]>([]);
  const [cities, setCities] = useState<any[]>([]);
  const [selectedCountry, setSelectedCountry] = useState<any>(null);
  const [selectedState, setSelectedState] = useState<any>(null);
  const [selectedCity, setSelectedCity] = useState<any>(null);

  // Load countries on component mount
  useEffect(() => {
    const allCountries = Country.getAllCountries();
    setCountries(allCountries);
  }, []);

  // Load states when country changes
  useEffect(() => {
    if (selectedCountry) {
      const countryStates = State.getStatesOfCountry(selectedCountry.isoCode);
      setStates(countryStates);
      setSelectedState(null);
      setSelectedCity(null);
      setCities([]);
    }
  }, [selectedCountry]);

  // Load cities when state changes
  useEffect(() => {
    if (selectedCountry && selectedState) {
      const stateCities = City.getCitiesOfState(selectedCountry.isoCode, selectedState.isoCode);
      setCities(stateCities);
      setSelectedCity(null);
    }
  }, [selectedState, selectedCountry]);

  const [passwordValid, setPasswordValid] = useState({
    length: false,
    complexity: false,
  });

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value, type } = e.target;
    const checked = (e.target as HTMLInputElement).checked;

    setFormData(prev => ({
      ...prev,
      [name]: type === 'checkbox' ? checked : value,
    }));

    // Check password validity
    if (name === 'password') {
      checkPasswordValidity(value);
    }
  };

  const handleSelectChange = (e: any) => {
    setFormData(prev => ({ ...prev, occupation: e.target.value }));
  };

  const checkPasswordValidity = (password: string) => {
    const hasLength = password.length >= 8;
    const hasLower = /[a-z]/.test(password);
    const hasUpper = /[A-Z]/.test(password);
    const hasNumber = /[0-9]/.test(password);
    const hasSpecial = /[!@#$%^&*]/.test(password);
    const complexityCount = [hasLower, hasUpper, hasNumber, hasSpecial].filter(Boolean).length;

    setPasswordValid({
      length: hasLength,
      complexity: complexityCount >= 3,
    });
  };

  const handleStep1Continue = () => {
    if (formData.email) {
      setStep(2);
    }
  };

  const handleStep2Continue = () => {
    // Validate age (must be at least 18)
    if (formData.dob) {
      const birthDate = new Date(formData.dob);
      const today = new Date();
      const age = today.getFullYear() - birthDate.getFullYear();
      const monthDiff = today.getMonth() - birthDate.getMonth();

      if (
        age < 18 ||
        (age === 18 && monthDiff < 0) ||
        (age === 18 && monthDiff === 0 && today.getDate() < birthDate.getDate())
      ) {
        setError('You must be at least 18 years old to register');
        return;
      }
    }

    if (
      formData.fullName &&
      formData.occupation &&
      formData.licenseId &&
      formData.dob &&
      formData.country &&
      formData.state &&
      formData.city &&
      formData.attestProfessional &&
      formData.agreeToTerms
    ) {
      setStep(3);
    } else {
      setError('Please complete all required fields');
    }
  };

  const handleSubmit = async () => {
    if (!passwordValid.length || !passwordValid.complexity) {
      setError('Password does not meet requirements');
      return;
    }

    if (formData.password !== formData.confirmPassword) {
      setError('Passwords do not match');
      return;
    }

    setIsLoading(true);
    setError(null);

    try {
      // Split full name into first and last name
      const nameParts = formData.fullName.trim().split(' ');
      const firstName = nameParts[0] || '';
      const lastName = nameParts.slice(1).join(' ') || '';

      const registerData = {
        email: formData.email,
        username: formData.email.split('@')[0],
        password: formData.password,
        password2: formData.confirmPassword,
        firstName: firstName,
        lastName: lastName,
        occupation: formData.occupation,
        agreeToTerms: formData.agreeToTerms,
        receiveInfo: false,
        dob: formData.dob,
        phone: '',
        country: formData.country,
        state: formData.state,
        city: formData.city,
      };

      const { accessToken, user } = await authService.register(registerData);
      useAuthStore.setState({ accessToken, isAuthorized: true, user });
      navigate('/dashboard');
    } catch (err: any) {
      setError(err.response?.data?.message || 'Registration failed. Please try again.');
    } finally {
      setIsLoading(false);
    }
  };

  const handleGoogleSuccess = async (res: any) => {
    setError(null);
    setIsLoading(true);
    try {
      const { accessToken, user } = await authService.loginWithGoogle(res.credential);
      useAuthStore.setState({ accessToken, isAuthorized: true, user });
      navigate('/dashboard');
    } catch (err: any) {
      setError(err.response?.data?.error || 'Google registration failed');
    } finally {
      setIsLoading(false);
    }
  };

  const handleGoogleError = () => {
    setError('Google registration failed');
  };

  return (
    <PageContainer>
      <Container maxWidth="sm">
        <RegisterCard>
          <LogoWrapper>
            <Link component={RouterLink} to="/">
              <Logo src={NovikLogo} alt="Novik Logo" />
            </Link>
          </LogoWrapper>

          {/* Step 1: Email */}
          {step === 1 && (
            <>
              <Title>Please register to access Novik</Title>
              <Subtitle>Use your professional email. We will verify your credentials.</Subtitle>

              <StyledTextField
                fullWidth
                name="email"
                type="email"
                label="Email"
                value={formData.email}
                onChange={handleChange}
                required
                autoFocus
              />

              <PrimaryButton onClick={handleStep1Continue} disabled={!formData.email}>
                Continue
              </PrimaryButton>

              <DividerText>OR</DividerText>

              <GoogleOAuthProvider clientId="415749549321-2g2mhh6ugbk8fhjfdcd4jo7sk00dfa8v.apps.googleusercontent.com">
                <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
                  <GoogleLogin
                    onSuccess={handleGoogleSuccess}
                    onError={handleGoogleError}
                    theme="outline"
                    width={370}
                    size="large"
                    text="continue_with"
                  />
                </Box>
              </GoogleOAuthProvider>

              <Box sx={{ mt: 2, textAlign: 'center' }}>
                <Typography variant="body2" sx={{ color: novikTheme.colors.textMuted }}>
                  Already have an account?{' '}
                  <LinkText component={RouterLink} to="/login">
                    Sign in
                  </LinkText>
                </Typography>
              </Box>
            </>
          )}

          {/* Step 2: Professional Information */}
          {step === 2 && (
            <>
              <Title>Complete your registration</Title>
              <Subtitle>We verify HCP status to provide safe clinical features.</Subtitle>

              <StyledTextField
                fullWidth
                name="fullName"
                type="text"
                label="Full name"
                value={formData.fullName}
                onChange={handleChange}
                required
                autoFocus
              />

              <StyledTextField
                fullWidth
                name="dob"
                type="date"
                label="Date of Birth"
                value={formData.dob}
                onChange={handleChange}
                slotProps={{
                  inputLabel: { shrink: true },
                  htmlInput: {
                    max: new Date().toISOString().split('T')[0], // Prevent future dates
                    min: new Date(new Date().setFullYear(new Date().getFullYear() - 100))
                      .toISOString()
                      .split('T')[0], // Max age 100
                  },
                }}
                required
              />

              <StyledFormControl fullWidth>
                <InputLabel>Occupation*</InputLabel>
                <Select
                  value={formData.occupation}
                  onChange={handleSelectChange}
                  label="Occupation*"
                  required
                >
                  <MenuItem value="">Select your role</MenuItem>
                  <MenuItem value="Dentist">Dentist</MenuItem>
                  <MenuItem value="Dental Hygienist">Dental Hygienist</MenuItem>
                  <MenuItem value="Dental Assistant">Dental Assistant</MenuItem>
                  <MenuItem value="Orthodontist">Orthodontist</MenuItem>
                  <MenuItem value="Oral Surgeon">Oral Surgeon</MenuItem>
                  <MenuItem value="Student">Student</MenuItem>
                  <MenuItem value="Other">Other Healthcare Professional</MenuItem>
                </Select>
              </StyledFormControl>

              <StyledTextField
                fullWidth
                name="licenseId"
                type="text"
                label="Professional ID"
                value={formData.licenseId}
                onChange={handleChange}
                required
              />

              <Autocomplete
                value={selectedCountry}
                onChange={(_, newValue) => {
                  setSelectedCountry(newValue);
                  setFormData(prev => ({ ...prev, country: newValue?.name || '' }));
                }}
                options={countries}
                getOptionLabel={option => option.name}
                renderInput={params => <StyledTextField {...params} label="Country" required />}
                sx={{ mb: '14px' }}
              />

              <Autocomplete
                value={selectedState}
                onChange={(_, newValue) => {
                  setSelectedState(newValue);
                  setFormData(prev => ({ ...prev, state: newValue?.name || '' }));
                }}
                options={states}
                getOptionLabel={option => option.name}
                renderInput={params => (
                  <StyledTextField {...params} label="State/Province" required />
                )}
                disabled={!selectedCountry}
                sx={{ mb: '14px' }}
              />

              <Autocomplete
                value={selectedCity}
                onChange={(_, newValue) => {
                  setSelectedCity(newValue);
                  setFormData(prev => ({ ...prev, city: newValue?.name || '' }));
                }}
                options={cities}
                getOptionLabel={option => option.name}
                renderInput={params => <StyledTextField {...params} label="City" required />}
                disabled={!selectedState}
                sx={{ mb: '14px' }}
              />

              <CheckboxGroup>
                <FormControlLabel
                  control={
                    <Checkbox
                      name="attestProfessional"
                      checked={formData.attestProfessional}
                      onChange={handleChange}
                      size="small"
                    />
                  }
                  label="I confirm I am a licensed healthcare professional."
                />
              </CheckboxGroup>

              <CheckboxGroup>
                <FormControlLabel
                  control={
                    <Checkbox
                      name="agreeToTerms"
                      checked={formData.agreeToTerms}
                      onChange={handleChange}
                      size="small"
                    />
                  }
                  label={
                    <span>
                      I agree to the{' '}
                      <LinkText component={RouterLink} to="/legal#terms">
                        Terms
                      </LinkText>{' '}
                      and{' '}
                      <LinkText component={RouterLink} to="/legal#privacy">
                        Privacy Policy
                      </LinkText>
                      .
                    </span>
                  }
                />
              </CheckboxGroup>

              {error && (
                <Typography
                  variant="body2"
                  sx={{
                    color: '#d14343',
                    mb: 2,
                    fontSize: '0.85rem',
                    textAlign: 'center',
                  }}
                >
                  {error}
                </Typography>
              )}

              <PrimaryButton onClick={handleStep2Continue}>Continue</PrimaryButton>

              <Box sx={{ mt: 2, textAlign: 'center' }}>
                <LinkText component="button" onClick={() => setStep(1)}>
                  ← Back
                </LinkText>
              </Box>
            </>
          )}

          {/* Step 3: Password */}
          {step === 3 && (
            <>
              <Title>Set your password</Title>
              <Subtitle>Create a secure password to finish.</Subtitle>

              <StyledTextField
                fullWidth
                name="password"
                type="password"
                label="Password*"
                value={formData.password}
                onChange={handleChange}
                required
                autoFocus
              />

              <PasswordChecklist>
                <ChecklistItem valid={passwordValid.length}>
                  {passwordValid.length ? <Check fontSize="small" /> : <Close fontSize="small" />}
                  At least 8 characters
                </ChecklistItem>
                <ChecklistItem valid={passwordValid.complexity}>
                  {passwordValid.complexity ? (
                    <Check fontSize="small" />
                  ) : (
                    <Close fontSize="small" />
                  )}
                  At least 3 of: lowercase, uppercase, numbers, special chars
                </ChecklistItem>
              </PasswordChecklist>

              <StyledTextField
                fullWidth
                name="confirmPassword"
                type="password"
                label="Confirm Password*"
                value={formData.confirmPassword}
                onChange={handleChange}
                required
              />

              {error && (
                <Typography
                  variant="body2"
                  sx={{
                    color: '#d14343',
                    mb: 2,
                    fontSize: '0.85rem',
                    textAlign: 'center',
                  }}
                >
                  {error}
                </Typography>
              )}

              <PrimaryButton
                onClick={handleSubmit}
                disabled={isLoading || !passwordValid.length || !passwordValid.complexity}
              >
                {isLoading ? <CircularProgress size={24} sx={{ color: 'white' }} /> : 'Finish'}
              </PrimaryButton>

              <Box sx={{ mt: 2, textAlign: 'center' }}>
                <LinkText component="button" onClick={() => setStep(2)}>
                  ← Back
                </LinkText>
              </Box>
            </>
          )}
        </RegisterCard>
      </Container>
    </PageContainer>
  );
}

export default Register;
