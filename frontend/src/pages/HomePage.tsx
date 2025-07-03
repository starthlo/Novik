import { Box, Typography } from '@mui/material';
import CookieConsent from '../components/Common/CookieConsent';
import NovikLogo from '../assets/Novik.png';

const HomePage = () => {
  return (
    <>
      <CookieConsent />
      <Box
        component="main"
        sx={{
          flexGrow: 1,
          overflowY: 'auto',
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          pt: 8,
          px: 2,
        }}
      >
        <Box component="img" src={NovikLogo} alt="Novik" sx={{ height: 80, mb: 2 }} />
        <Typography variant="h5" align="center" sx={{ color: '#374151', maxWidth: 600 }}>
          Your smart AI Dental assistant for safe clinical decisions
        </Typography>
      </Box>
    </>
  );
};

export default HomePage;
