import { useState } from 'react';
import { Box, Container, Typography, Link, Grid, Paper } from '@mui/material';
import osstemLogo from '../../assets/ossem_logo.png';
// import other logos here

const partners = [
  {
    name: 'Osstem',
    logo: osstemLogo,
    url: 'https://marketing.osstem.es/',
  },
];

const Partners = () => {
  const [clicks, setClicks] = useState(0);

  const handleClick = () => {
    setClicks(prev => prev + 1);
    fetch('/api/banner-click', { method: 'POST' });
  };

  return (
    <Box component="section" sx={{ py: 6 }}>
      <Container maxWidth="lg">
        <Typography
          align="center"
          gutterBottom
          sx={{ fontSize: 20, fontWeight: 500, color: '#333' }}
        >
          Our Partners
        </Typography>
        <Grid container spacing={4} justifyContent="center" alignItems="center" sx={{ mt: 2 }}>
          {partners.map(({ name, logo, url }) => (
            <Grid key={name} size={{ xs: 6, sm: 4, md: 3, lg: 2 }}>
              <Paper
                elevation={1}
                onClick={handleClick}
                component={Link}
                href={url}
                target="_blank"
                rel="noopener noreferrer"
                sx={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  p: 2,
                  cursor: 'pointer',
                  transition: 'transform 0.2s',
                  '&:hover': { transform: 'scale(1.05)' },
                }}
              >
                <Box
                  component="img"
                  src={logo}
                  alt={name}
                  sx={{ maxHeight: 64, objectFit: 'contain' }}
                />
              </Paper>
            </Grid>
          ))}
        </Grid>
        <Typography align="center" sx={{ fontSize: 12, color: '#666', mt: 3 }}>
          Clicks: {clicks}
        </Typography>
      </Container>
    </Box>
  );
};

export default Partners;
