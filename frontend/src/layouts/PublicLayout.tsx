import { Outlet } from 'react-router-dom';
import { Box } from '@mui/material';
import Header from '../components/Common/Header';
import PoweredBy from '../components/Common/PoweredBy';
import Partners from '../components/Common/Partners';
import Footer from '../components/Common/Footer';
import FrontImage from '../assets/Front Image.png';

const PublicLayout = () => {
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
      <Box sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
        <Outlet />
      </Box>
      <PoweredBy />
      <Partners />
      <Footer />
    </Box>
  );
};

export default PublicLayout;
