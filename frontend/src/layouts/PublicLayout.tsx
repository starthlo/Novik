import { Outlet } from 'react-router-dom';
import { Box } from '@mui/material';
import Header from '../components/Common/Header';
import FrontImage from '../assets/Front Image.png';

const PublicLayout = () => {
  return (
    <Box
      sx={{
        backgroundImage: `url(${FrontImage})`,
        backgroundSize: 'cover',
        backgroundPosition: 'center',
        backgroundRepeat: 'no-repeat',
        minHeight: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <Header />
      <Box sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
        <Outlet />
      </Box>
    </Box>
  );
};

export default PublicLayout;
