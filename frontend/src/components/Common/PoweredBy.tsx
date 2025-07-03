import { Box, Container, Typography } from '@mui/material';
import chatgptLogo from '../../assets/chatgpt.png';
import drugbankLogo from '../../assets/drugbank.png';
import pubmedLogo from '../../assets/pubmed.png';

const PoweredBy = () => (
  <Box
    component="section"
    sx={{
      mt: 6,
      textAlign: 'center',
      px: 2,
    }}
  >
    <Container maxWidth="md">
      <Typography
        variant="h6"
        component="h2"
        sx={{
          color: 'text.secondary',
          fontWeight: 500,
          mb: 2,
        }}
      >
        Powered By
      </Typography>

      <Box
        sx={{
          display: 'flex',
          flexDirection: { xs: 'column', sm: 'row' },
          justifyContent: 'center',
          alignItems: 'center',
          gap: { xs: 3, sm: 4 },
          flexWrap: 'wrap',
        }}
      >
        {[chatgptLogo, drugbankLogo, pubmedLogo].map((src, idx) => (
          <Box
            key={idx}
            component="img"
            src={src}
            alt=""
            sx={{
              width: { xs: 96, sm: 128, md: 192 },
              objectFit: 'contain',
            }}
          />
        ))}
      </Box>
    </Container>
  </Box>
);

export default PoweredBy;
