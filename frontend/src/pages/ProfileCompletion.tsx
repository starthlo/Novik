import { useState, useEffect } from 'react';
import {
  Box,
  Container,
  Typography,
  TextField,
  Button,
  Alert,
  Paper,
  Grid,
  Checkbox,
  FormControlLabel,
  FormHelperText,
  styled,
  CircularProgress,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  Autocomplete,
} from '@mui/material';
import { useNavigate } from 'react-router-dom';
import { Country, State, City } from 'country-state-city';
import { userService } from '../services/userService';
import { useAuthStore } from '../stores/auth';
import { novikTheme } from '../styles/theme';
import NovikLogo from '../assets/novik-logo.png';

const PageContainer = styled(Box)({
  minHeight: '100vh',
  backgroundColor: '#f7f7f8',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '2rem',
});

const CompletionCard = styled(Paper)({
  maxWidth: '800px',
  width: '100%',
  padding: '3rem',
  borderRadius: '16px',
  boxShadow: '0 20px 50px rgba(0,0,0,0.08)',
});

const LogoWrapper = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  marginBottom: '2rem',
});

const Logo = styled('img')({
  width: '100px',
  height: '100px',
  marginBottom: '1rem',
});

const Title = styled(Typography)({
  fontSize: '1.8rem',
  fontWeight: 700,
  marginBottom: '0.5rem',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
  textAlign: 'center',
});

const Subtitle = styled(Typography)({
  fontSize: '1rem',
  color: novikTheme.colors.textMuted,
  fontFamily: novikTheme.typography.fontFamily,
  textAlign: 'center',
  marginBottom: '2rem',
});

const StyledTextField = styled(TextField)({
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

const SubmitButton = styled(Button)({
  padding: '0.75rem 2rem',
  borderRadius: '10px',
  fontWeight: 600,
  fontSize: '1rem',
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

const ProfileCompletion = () => {
  const navigate = useNavigate();
  const { user } = useAuthStore();
  const [loading, setLoading] = useState(false);
  const [errors, setErrors] = useState<Record<string, string>>({});

  const [formData, setFormData] = useState({
    occupation: '',
    dob: '',
    phone: '',
    country: '',
    state: '',
    city: '',
    agreeToTerms: false,
    receiveInfo: false,
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

  useEffect(() => {
    // If user is not logged in, redirect to login
    if (!user) {
      navigate('/login');
      return;
    }

    // If profile is already completed, redirect to dashboard
    if (user.profileCompleted) {
      navigate('/dashboard');
      return;
    }

    // Pre-fill any existing data
    setFormData({
      occupation: user.occupation || '',
      dob: user.dob || '',
      phone: user.phone || '',
      country: user.country || '',
      state: user.state || '',
      city: user.city || '',
      agreeToTerms: user.agreeToTerms || false,
      receiveInfo: user.receiveInfo || false,
    });

    // Pre-select location if available
    if (user.country) {
      const country = countries.find(c => c.name === user.country);
      if (country) {
        setSelectedCountry(country);
      }
    }
  }, [user, navigate, countries]);

  const handleChange = (field: string, value: any) => {
    setFormData(prev => ({ ...prev, [field]: value }));
    // Clear error for this field when user starts typing
    if (errors[field]) {
      setErrors(prev => ({ ...prev, [field]: '' }));
    }
  };

  const validateForm = () => {
    const newErrors: Record<string, string> = {};

    if (!formData.occupation) newErrors.occupation = 'Occupation is required';
    if (!formData.dob) newErrors.dob = 'Date of birth is required';
    if (!formData.phone) newErrors.phone = 'Phone number is required';
    if (!formData.country) newErrors.country = 'Country is required';
    if (!formData.state) newErrors.state = 'State/Province is required';
    if (!formData.city) newErrors.city = 'City is required';
    if (!formData.agreeToTerms)
      newErrors.agreeToTerms = 'You must agree to the terms and conditions';

    if (formData.phone && !/^[\d\s\-\+\(\)]+$/.test(formData.phone)) {
      newErrors.phone = 'Please enter a valid phone number';
    }

    // Validate date of birth (must be in the past and reasonable age)
    if (formData.dob) {
      const dob = new Date(formData.dob);
      const today = new Date();
      const age = today.getFullYear() - dob.getFullYear();

      if (dob > today) {
        newErrors.dob = 'Date of birth must be in the past';
      } else if (age < 18) {
        newErrors.dob = 'You must be at least 18 years old';
      } else if (age > 120) {
        newErrors.dob = 'Please enter a valid date of birth';
      }
    }

    setErrors(newErrors);
    return Object.keys(newErrors).length === 0;
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    if (!validateForm()) {
      return;
    }

    try {
      setLoading(true);
      const updatedProfile = await userService.updateProfile(formData);

      // Update the user in auth store
      const authStore = useAuthStore.getState();
      authStore.user = { ...authStore.user, ...updatedProfile };

      // Redirect to dashboard
      navigate('/dashboard');
    } catch (error: any) {
      console.error('Error updating profile:', error);
      setErrors({
        submit: error.response?.data?.error || 'Failed to update profile. Please try again.',
      });
    } finally {
      setLoading(false);
    }
  };

  return (
    <PageContainer>
      <Container maxWidth="md">
        <CompletionCard>
          <LogoWrapper>
            <Logo src={NovikLogo} alt="Novik Logo" />
            <Title>Complete Your Profile</Title>
            <Subtitle>
              Welcome to Novik! Please provide the following information to complete your
              registration. This information helps us provide you with a better experience.
            </Subtitle>
          </LogoWrapper>

          <Box component="form" onSubmit={handleSubmit}>
            <Grid container spacing={3}>
              <Grid size={{ xs: 12 }}>
                <StyledFormControl fullWidth required error={!!errors.occupation}>
                  <InputLabel>Occupation</InputLabel>
                  <Select
                    value={formData.occupation}
                    onChange={e => handleChange('occupation', e.target.value)}
                    label="Occupation"
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
                  {errors.occupation && <FormHelperText>{errors.occupation}</FormHelperText>}
                </StyledFormControl>
              </Grid>

              <Grid size={{ xs: 12, md: 6 }}>
                <StyledTextField
                  fullWidth
                  required
                  label="Date of Birth"
                  type="date"
                  value={formData.dob}
                  onChange={e => handleChange('dob', e.target.value)}
                  error={!!errors.dob}
                  helperText={errors.dob}
                  InputLabelProps={{ shrink: true }}
                />
              </Grid>

              <Grid size={{ xs: 12, md: 6 }}>
                <StyledTextField
                  fullWidth
                  required
                  label="Phone Number"
                  placeholder="+1 (555) 123-4567"
                  value={formData.phone}
                  onChange={e => handleChange('phone', e.target.value)}
                  error={!!errors.phone}
                  helperText={errors.phone}
                />
              </Grid>

              <Grid size={{ xs: 12, md: 4 }}>
                <Autocomplete
                  value={selectedCountry}
                  onChange={(_, newValue) => {
                    setSelectedCountry(newValue);
                    handleChange('country', newValue?.name || '');
                  }}
                  options={countries}
                  getOptionLabel={option => option.name}
                  renderInput={params => (
                    <StyledTextField
                      {...params}
                      label="Country"
                      required
                      error={!!errors.country}
                      helperText={errors.country}
                    />
                  )}
                />
              </Grid>

              <Grid size={{ xs: 12, md: 4 }}>
                <Autocomplete
                  value={selectedState}
                  onChange={(_, newValue) => {
                    setSelectedState(newValue);
                    handleChange('state', newValue?.name || '');
                  }}
                  options={states}
                  getOptionLabel={option => option.name}
                  renderInput={params => (
                    <StyledTextField
                      {...params}
                      label="State/Province"
                      required
                      error={!!errors.state}
                      helperText={errors.state}
                    />
                  )}
                  disabled={!selectedCountry}
                />
              </Grid>

              <Grid size={{ xs: 12, md: 4 }}>
                <Autocomplete
                  value={selectedCity}
                  onChange={(_, newValue) => {
                    setSelectedCity(newValue);
                    handleChange('city', newValue?.name || '');
                  }}
                  options={cities}
                  getOptionLabel={option => option.name}
                  renderInput={params => (
                    <StyledTextField
                      {...params}
                      label="City"
                      required
                      error={!!errors.city}
                      helperText={errors.city}
                    />
                  )}
                  disabled={!selectedState}
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={formData.agreeToTerms}
                      onChange={e => handleChange('agreeToTerms', e.target.checked)}
                      color="primary"
                    />
                  }
                  label={
                    <Typography
                      sx={{ fontSize: '0.9rem', fontFamily: novikTheme.typography.fontFamily }}
                    >
                      I agree to the Terms of Service and Privacy Policy *
                    </Typography>
                  }
                />
                {errors.agreeToTerms && (
                  <FormHelperText error>{errors.agreeToTerms}</FormHelperText>
                )}
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={formData.receiveInfo}
                      onChange={e => handleChange('receiveInfo', e.target.checked)}
                      color="primary"
                    />
                  }
                  label={
                    <Typography
                      sx={{ fontSize: '0.9rem', fontFamily: novikTheme.typography.fontFamily }}
                    >
                      I would like to receive updates and information about Novik
                    </Typography>
                  }
                />
              </Grid>
            </Grid>

            {errors.submit && (
              <Alert severity="error" sx={{ mt: 2 }}>
                {errors.submit}
              </Alert>
            )}

            <Box sx={{ mt: 4, display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
              <SubmitButton type="submit" disabled={loading}>
                {loading ? (
                  <CircularProgress size={24} sx={{ color: 'white' }} />
                ) : (
                  'Complete Profile'
                )}
              </SubmitButton>
            </Box>
          </Box>
        </CompletionCard>
      </Container>
    </PageContainer>
  );
};

export default ProfileCompletion;
