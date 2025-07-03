import { Box, Typography } from '@mui/material';
import Header from '../components/Common/Header';
import PoweredBy from '../components/Common/PoweredBy';
import Partners from '../components/Common/Partners';
import Footer from '../components/Common/Footer';
import CookieConsent from '../components/Common/CookieConsent';
import NovikLogo from '../assets/Novik.png';
import FrontImage from '../assets/Front Image.png';

const HomePage = () => {
  return (
    <Box
      sx={{
        backgroundImage: `url(${FrontImage})`,
        backgroundSize: 'cover',
        backgroundPosition: 'center',
        backgroundRepeat: 'no-repeat',
        width: '100vw',
        minHeight: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <Header />
      <CookieConsent />
      <Box sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
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
        <PoweredBy />
        <Partners />
        <Footer />
      </Box>
    </Box>
  );
};

export default HomePage;
