import { useState, useEffect } from 'react';
import {
  Box,
  Container,
  Typography,
  TextField,
  Button,
  Alert,
  Snackbar,
  Paper,
  Grid,
  CircularProgress,
  styled,
  IconButton,
  InputAdornment,
} from '@mui/material';
import { Visibility, VisibilityOff, Edit } from '@mui/icons-material';
import { useNavigate } from 'react-router-dom';
import { userService, UserProfile, ChangePasswordData } from '../services/userService';
import { useAuthStore } from '../stores/auth';
import { novikTheme } from '../styles/theme';

const PageContainer = styled(Box)({
  minHeight: '100vh',
  backgroundColor: '#f7f7f8',
  paddingTop: '2rem',
  paddingBottom: '4rem',
});

const SectionPaper = styled(Paper)({
  padding: '2rem',
  marginBottom: '2rem',
  borderRadius: '12px',
  boxShadow: '0 2px 8px rgba(0,0,0,0.08)',
});

const SectionTitle = styled(Typography)({
  fontSize: '1.5rem',
  fontWeight: 600,
  marginBottom: '1.5rem',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
});

const FormButton = styled(Button)({
  padding: '0.5rem 1.5rem',
  borderRadius: '8px',
  fontWeight: 600,
  fontSize: '0.95rem',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
});

const SaveButton = styled(FormButton)({
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const CancelButton = styled(FormButton)({
  backgroundColor: '#f5f5f5',
  color: novikTheme.colors.text,
  marginRight: '1rem',
  '&:hover': {
    backgroundColor: '#e0e0e0',
  },
});

type AlertType = {
  open: boolean;
  message: string;
  severity: 'error' | 'warning' | 'info' | 'success';
};

const Account = () => {
  const navigate = useNavigate();
  const { logout } = useAuthStore();
  const [loading, setLoading] = useState(true);
  const [editing, setEditing] = useState(false);
  const [changingPassword, setChangingPassword] = useState(false);

  const [profile, setProfile] = useState<UserProfile | null>(null);
  const [editedProfile, setEditedProfile] = useState<Partial<UserProfile>>({});

  const [passwordData, setPasswordData] = useState<ChangePasswordData>({
    current_password: '',
    new_password: '',
    confirm_password: '',
  });

  const [showPasswords, setShowPasswords] = useState({
    current: false,
    new: false,
    confirm: false,
  });

  const [alert, setAlert] = useState<AlertType>({
    open: false,
    message: '',
    severity: 'info',
  });

  useEffect(() => {
    loadProfile();
  }, []);

  const loadProfile = async () => {
    try {
      setLoading(true);
      const data = await userService.getProfile();
      setProfile(data);
      setEditedProfile(data);
    } catch (error) {
      showAlert('Failed to load profile', 'error');
      console.error('Error loading profile:', error);
    } finally {
      setLoading(false);
    }
  };

  const showAlert = (message: string, severity: AlertType['severity'] = 'info') => {
    setAlert({ open: true, message, severity });
  };

  const handleAlertClose = () => {
    setAlert(prev => ({ ...prev, open: false }));
  };

  const handleEditToggle = () => {
    if (editing) {
      // Cancel editing - restore original values
      setEditedProfile(profile || {});
    }
    setEditing(!editing);
  };

  const handleProfileChange = (field: keyof UserProfile, value: string) => {
    setEditedProfile(prev => ({ ...prev, [field]: value }));
  };

  const handleProfileSave = async () => {
    try {
      setLoading(true);
      const updatedProfile = await userService.updateProfile(editedProfile);
      setProfile(updatedProfile);
      setEditing(false);
      showAlert('Profile updated successfully', 'success');

      // Update the user in auth store
      const authStore = useAuthStore.getState();
      authStore.user = { ...authStore.user, ...updatedProfile };
    } catch (error: any) {
      showAlert(error.response?.data?.error || 'Failed to update profile', 'error');
    } finally {
      setLoading(false);
    }
  };

  const handlePasswordChange = (field: keyof ChangePasswordData, value: string) => {
    setPasswordData(prev => ({ ...prev, [field]: value }));
  };

  const handlePasswordSave = async () => {
    // Validate passwords
    if (
      !passwordData.current_password ||
      !passwordData.new_password ||
      !passwordData.confirm_password
    ) {
      showAlert('All password fields are required', 'error');
      return;
    }

    if (passwordData.new_password !== passwordData.confirm_password) {
      showAlert('New passwords do not match', 'error');
      return;
    }

    if (passwordData.new_password.length < 8) {
      showAlert('New password must be at least 8 characters long', 'error');
      return;
    }

    try {
      setLoading(true);
      const response = await userService.changePassword(passwordData);
      showAlert(response.message, 'success');
      setChangingPassword(false);
      setPasswordData({
        current_password: '',
        new_password: '',
        confirm_password: '',
      });
    } catch (error: any) {
      showAlert(error.response?.data?.error || 'Failed to change password', 'error');
    } finally {
      setLoading(false);
    }
  };

  const togglePasswordVisibility = (field: 'current' | 'new' | 'confirm') => {
    setShowPasswords(prev => ({ ...prev, [field]: !prev[field] }));
  };

  if (loading && !profile) {
    return (
      <PageContainer>
        <Container maxWidth="md">
          <Box sx={{ display: 'flex', justifyContent: 'center', py: 4 }}>
            <CircularProgress />
          </Box>
        </Container>
      </PageContainer>
    );
  }

  return (
    <PageContainer>
      <Container maxWidth="md">
        <Typography
          variant="h4"
          sx={{
            fontWeight: 700,
            mb: 3,
            color: novikTheme.colors.text,
            fontFamily: novikTheme.typography.fontFamily,
          }}
        >
          Account Settings
        </Typography>

        {/* Profile Information Section */}
        <SectionPaper>
          <Box
            sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}
          >
            <SectionTitle>Profile Information</SectionTitle>
            {!editing ? (
              <IconButton onClick={handleEditToggle} color="primary">
                <Edit />
              </IconButton>
            ) : null}
          </Box>

          <Grid container spacing={3}>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="Username"
                value={profile?.username || ''}
                disabled
                helperText="Username cannot be changed"
              />
            </Grid>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="Email"
                type="email"
                value={editedProfile.email || ''}
                onChange={e => handleProfileChange('email', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="First Name"
                value={editedProfile.firstName || ''}
                onChange={e => handleProfileChange('firstName', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="Last Name"
                value={editedProfile.lastName || ''}
                onChange={e => handleProfileChange('lastName', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="Date of Birth"
                type="date"
                value={editedProfile.dob || ''}
                onChange={e => handleProfileChange('dob', e.target.value)}
                disabled={!editing}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 6 }}>
              <TextField
                fullWidth
                label="Phone"
                value={editedProfile.phone || ''}
                onChange={e => handleProfileChange('phone', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12 }}>
              <TextField
                fullWidth
                label="Occupation"
                value={editedProfile.occupation || ''}
                onChange={e => handleProfileChange('occupation', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 4 }}>
              <TextField
                fullWidth
                label="Country"
                value={editedProfile.country || ''}
                onChange={e => handleProfileChange('country', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 4 }}>
              <TextField
                fullWidth
                label="State/Province"
                value={editedProfile.state || ''}
                onChange={e => handleProfileChange('state', e.target.value)}
                disabled={!editing}
              />
            </Grid>
            <Grid size={{ xs: 12, md: 4 }}>
              <TextField
                fullWidth
                label="City"
                value={editedProfile.city || ''}
                onChange={e => handleProfileChange('city', e.target.value)}
                disabled={!editing}
              />
            </Grid>
          </Grid>

          {editing && (
            <Box sx={{ mt: 3, display: 'flex', justifyContent: 'flex-end' }}>
              <CancelButton onClick={handleEditToggle}>Cancel</CancelButton>
              <SaveButton onClick={handleProfileSave} disabled={loading}>
                Save Changes
              </SaveButton>
            </Box>
          )}
        </SectionPaper>

        {/* Change Password Section */}
        <SectionPaper>
          <SectionTitle>Change Password</SectionTitle>

          {!changingPassword ? (
            <Button
              variant="outlined"
              onClick={() => setChangingPassword(true)}
              sx={{
                borderColor: novikTheme.colors.primary,
                color: novikTheme.colors.primary,
                '&:hover': {
                  borderColor: novikTheme.colors.primaryDark,
                  backgroundColor: 'rgba(136, 169, 78, 0.05)',
                },
              }}
            >
              Change Password
            </Button>
          ) : (
            <>
              <Grid container spacing={3}>
                <Grid size={{ xs: 12 }}>
                  <TextField
                    fullWidth
                    label="Current Password"
                    type={showPasswords.current ? 'text' : 'password'}
                    value={passwordData.current_password}
                    onChange={e => handlePasswordChange('current_password', e.target.value)}
                    InputProps={{
                      endAdornment: (
                        <InputAdornment position="end">
                          <IconButton
                            onClick={() => togglePasswordVisibility('current')}
                            edge="end"
                          >
                            {showPasswords.current ? <VisibilityOff /> : <Visibility />}
                          </IconButton>
                        </InputAdornment>
                      ),
                    }}
                  />
                </Grid>
                <Grid size={{ xs: 12 }}>
                  <TextField
                    fullWidth
                    label="New Password"
                    type={showPasswords.new ? 'text' : 'password'}
                    value={passwordData.new_password}
                    onChange={e => handlePasswordChange('new_password', e.target.value)}
                    InputProps={{
                      endAdornment: (
                        <InputAdornment position="end">
                          <IconButton onClick={() => togglePasswordVisibility('new')} edge="end">
                            {showPasswords.new ? <VisibilityOff /> : <Visibility />}
                          </IconButton>
                        </InputAdornment>
                      ),
                    }}
                  />
                </Grid>
                <Grid size={{ xs: 12 }}>
                  <TextField
                    fullWidth
                    label="Confirm New Password"
                    type={showPasswords.confirm ? 'text' : 'password'}
                    value={passwordData.confirm_password}
                    onChange={e => handlePasswordChange('confirm_password', e.target.value)}
                    InputProps={{
                      endAdornment: (
                        <InputAdornment position="end">
                          <IconButton
                            onClick={() => togglePasswordVisibility('confirm')}
                            edge="end"
                          >
                            {showPasswords.confirm ? <VisibilityOff /> : <Visibility />}
                          </IconButton>
                        </InputAdornment>
                      ),
                    }}
                  />
                </Grid>
              </Grid>

              <Box sx={{ mt: 3, display: 'flex', justifyContent: 'flex-end' }}>
                <CancelButton
                  onClick={() => {
                    setChangingPassword(false);
                    setPasswordData({
                      current_password: '',
                      new_password: '',
                      confirm_password: '',
                    });
                  }}
                >
                  Cancel
                </CancelButton>
                <SaveButton onClick={handlePasswordSave} disabled={loading}>
                  Update Password
                </SaveButton>
              </Box>
            </>
          )}
        </SectionPaper>

        {/* Account Actions Section */}
        <SectionPaper>
          <SectionTitle>Account Actions</SectionTitle>
          <Typography sx={{ mb: 2, color: novikTheme.colors.textMuted }}>
            Need to sign out or delete your account?
          </Typography>
          <Box sx={{ display: 'flex', gap: 2 }}>
            <Button
              variant="outlined"
              color="error"
              onClick={() => {
                logout();
                navigate('/login');
              }}
            >
              Sign Out
            </Button>
          </Box>
        </SectionPaper>
      </Container>

      <Snackbar
        open={alert.open}
        autoHideDuration={4000}
        onClose={handleAlertClose}
        anchorOrigin={{ vertical: 'top', horizontal: 'center' }}
      >
        <Alert
          onClose={handleAlertClose}
          severity={alert.severity}
          variant="filled"
          sx={{
            width: '100%',
            ...(alert.severity === 'success' && {
              backgroundColor: novikTheme.colors.primary,
              color: '#ffffff',
              '& .MuiAlert-icon': {
                color: '#ffffff',
              },
            }),
          }}
        >
          {alert.message}
        </Alert>
      </Snackbar>
    </PageContainer>
  );
};

export default Account;
