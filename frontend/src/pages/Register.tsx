import { useEffect, useState } from 'react';
import { Country, State, City, ICountry, IState, ICity } from 'country-state-city';
import Header from '../components/Common/Header';
import {
  Box,
  Checkbox,
  Container,
  FormControl,
  FormControlLabel,
  Grid,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  SelectChangeEvent,
  Snackbar,
  TextField,
  Typography,
} from '@mui/material';
import { LoadingButton } from '@mui/lab';
import { authService } from '../services';
import { useAuthStore } from '../stores/auth';
import FrontImage from '../assets/Front Image.png';

export default function Register() {
  const [formData, setFormData] = useState({
    username: '',
    email: '',
    password: '',
    confirm_password: '',
    dob: '',
    phone: '',
    occupation: '',
    country: '',
    state: '',
    city: '',
    agreeToTerms: false,
    receiveInfo: false,
  });
  const [countries] = useState<ICountry[]>(Country.getAllCountries());
  const [states, setStates] = useState<IState[]>([]);
  const [cities, setCities] = useState<ICity[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [snackbar, setSnackbar] = useState<{
    open: boolean;
    message: string;
    severity: 'success' | 'error';
  }>({
    open: false,
    message: '',
    severity: 'success',
  });

  useEffect(() => {
    if (formData.country) {
      const sel = countries.find(c => c.name === formData.country);
      setStates(sel ? State.getStatesOfCountry(sel.isoCode) : []);
      if (!sel) setFormData(f => ({ ...f, state: '', city: '' }));
    }
  }, [formData.country]);

  useEffect(() => {
    const countryObj = countries.find(c => c.name === formData.country);
    const stateObj = states.find(s => s.name === formData.state);
    setCities(
      countryObj && stateObj ? City.getCitiesOfState(countryObj.isoCode, stateObj.isoCode) : []
    );
    if (!stateObj) setFormData(f => ({ ...f, city: '' }));
  }, [formData.state, countries, states]);

  const handleChange = (
    e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement> | SelectChangeEvent<string>
  ) => {
    const { name, value, type, checked } = e.target as HTMLInputElement & {
      name: string;
      value: any;
    };
    setFormData(f => ({
      ...f,
      [name]: type === 'checkbox' ? checked : value,
    }));
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsLoading(true);
    try {
      const { data } = await authService.register(formData);
      useAuthStore.setState({ accessToken: data.accessToken, isAuthorized: true, user: data.user });
      setSnackbar({ open: true, message: 'Registration successful!', severity: 'success' });
    } catch (err: any) {
      console.error(err);
      setSnackbar({
        open: true,
        message: err.response?.data?.message || 'Registration failed',
        severity: 'error',
      });
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Box
      sx={{
        backgroundImage: `url(${FrontImage})`,
        backgroundSize: 'cover',
        backgroundPosition: 'center',
        backgroundRepeat: 'no-repeat',
        minHeight: '100vh',
      }}
    >
      <Header />
      <Container maxWidth="sm" sx={{ py: 4 }}>
        <Paper sx={{ p: 4 }}>
          <Typography variant="h4" align="center" gutterBottom>
            User Registration
          </Typography>
          <Typography color="text.secondary" align="center" paragraph>
            To use the NOVIK Dental Assistant you need to register. Once you have registered you
            will be able to use Novik and discover its full potential.
          </Typography>

          <Box component="form" onSubmit={handleSubmit} noValidate>
            <Grid container spacing={2}>
              <Grid size={{ xs: 8 }}>
                <TextField
                  label="Name"
                  name="username"
                  value={formData.username}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>
              <Grid size={{ xs: 4 }}>
                <TextField
                  label="Date of Birth"
                  name="dob"
                  type="date"
                  InputLabelProps={{ shrink: true }}
                  value={formData.dob}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>

              <Grid size={{ xs: 6 }}>
                <TextField
                  label="Email"
                  name="email"
                  type="email"
                  value={formData.email}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>
              <Grid size={{ xs: 6 }}>
                <TextField
                  label="Phone"
                  name="phone"
                  value={formData.phone}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>

              <Grid size={{ xs: 6 }}>
                <TextField
                  label="Password"
                  name="password"
                  type="password"
                  value={formData.password}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>
              <Grid size={{ xs: 6 }}>
                <TextField
                  label="Confirm Password"
                  name="confirm_password"
                  type="password"
                  value={formData.confirm_password}
                  onChange={handleChange}
                  fullWidth
                  disabled={isLoading}
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControl fullWidth disabled={isLoading}>
                  <InputLabel>Occupation</InputLabel>
                  <Select
                    label="Occupation"
                    name="occupation"
                    value={formData.occupation}
                    onChange={handleChange}
                  >
                    <MenuItem value="">
                      <em>None</em>
                    </MenuItem>
                    <MenuItem value="General Dentistry">General Dentistry</MenuItem>
                    <MenuItem value="Endodontics">Endodontics</MenuItem>
                    <MenuItem value="Orthodontics">Orthodontics</MenuItem>
                    <MenuItem value="Oral and Maxillofacial Surgery">
                      Oral & Maxillofacial Surgery
                    </MenuItem>
                    <MenuItem value="Pediatric Dentistry">Pediatric Dentistry</MenuItem>
                    <MenuItem value="Prosthodontics / Oral Rehabilitation">
                      Prosthodontics / Oral Rehab
                    </MenuItem>
                    <MenuItem value="Cosmetic Dentistry">Cosmetic Dentistry</MenuItem>
                    <MenuItem value="Preventive and Community Dentistry">
                      Preventive & Community Dentistry
                    </MenuItem>
                    <MenuItem value="Digital Dentistry / CAD-CAM">
                      Digital Dentistry / CAD-CAM
                    </MenuItem>
                    <MenuItem value="Student">Student</MenuItem>
                  </Select>
                </FormControl>
              </Grid>

              <Grid size={{ xs: 4 }}>
                <FormControl fullWidth disabled={isLoading}>
                  <InputLabel>Country</InputLabel>
                  <Select
                    label="Country"
                    name="country"
                    value={formData.country}
                    onChange={handleChange}
                  >
                    <MenuItem value="">
                      <em>None</em>
                    </MenuItem>
                    {countries.map(c => (
                      <MenuItem key={c.isoCode} value={c.name}>
                        {c.name}
                      </MenuItem>
                    ))}
                  </Select>
                </FormControl>
              </Grid>
              <Grid size={{ xs: 4 }}>
                <FormControl fullWidth disabled={isLoading || !states.length}>
                  <InputLabel>State</InputLabel>
                  <Select label="State" name="state" value={formData.state} onChange={handleChange}>
                    <MenuItem value="">
                      <em>None</em>
                    </MenuItem>
                    {states.map(s => (
                      <MenuItem key={s.isoCode} value={s.name}>
                        {s.name}
                      </MenuItem>
                    ))}
                  </Select>
                </FormControl>
              </Grid>
              <Grid size={{ xs: 4 }}>
                <FormControl fullWidth disabled={isLoading || !cities.length}>
                  <InputLabel>City</InputLabel>
                  <Select label="City" name="city" value={formData.city} onChange={handleChange}>
                    <MenuItem value="">
                      <em>None</em>
                    </MenuItem>
                    {cities.map(city => (
                      <MenuItem key={city.name} value={city.name}>
                        {city.name}
                      </MenuItem>
                    ))}
                  </Select>
                </FormControl>
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControlLabel
                  control={
                    <Checkbox
                      name="agreeToTerms"
                      checked={formData.agreeToTerms}
                      onChange={handleChange}
                      disabled={isLoading}
                    />
                  }
                  label={
                    <Typography variant="body2">
                      I agree to the{' '}
                      <a href="/terms" style={{ color: '#FF9800' }}>
                        Legal Notice
                      </a>{' '}
                      and{' '}
                      <a href="/privacy" style={{ color: '#FF9800' }}>
                        Privacy Policy
                      </a>
                    </Typography>
                  }
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControlLabel
                  control={
                    <Checkbox
                      name="receiveInfo"
                      checked={formData.receiveInfo}
                      onChange={handleChange}
                      disabled={isLoading}
                    />
                  }
                  label="I would like to receive commercial information"
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <LoadingButton
                  type="submit"
                  variant="contained"
                  fullWidth
                  loading={isLoading}
                  style={{ backgroundColor: '#F97316' }}
                >
                  Submit
                </LoadingButton>
              </Grid>
            </Grid>
          </Box>
        </Paper>
      </Container>

      <Snackbar
        open={snackbar.open}
        autoHideDuration={4000}
        onClose={() => setSnackbar(s => ({ ...s, open: false }))}
        message={snackbar.message}
      />
    </Box>
  );
}
