import { useState } from 'react';
import {
  Box,
  Container,
  Typography,
  Grid,
  Paper,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Button,
  SelectChangeEvent,
} from '@mui/material';
import Header from '../components/Common/Header';
import PoweredBy from '../components/Common/PoweredBy';
import Partners from '../components/Common/Partners';
import Footer from '../components/Common/Footer';
import FrontImage from '../assets/Front Image.png';

interface FormData {
  name: string;
  email: string;
  category: string;
  phone: string;
  message: string;
}

const ContactUs = () => {
  const [formData, setFormData] = useState<FormData>({
    name: '',
    email: '',
    category: '',
    phone: '',
    message: '',
  });

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
  };

  const handleSelectChange = (e: SelectChangeEvent<string>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name as string]: value }));
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    console.log(formData);
    alert('Form submitted successfully!');
  };

  return (
    <Box
      sx={{
        backgroundImage: `url(${FrontImage})`,
        backgroundSize: 'cover',
        backgroundPosition: 'center',
        minHeight: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <Header />

      <Container maxWidth="md" sx={{ flexGrow: 1, py: 6 }}>
        <Paper
          elevation={4}
          sx={{
            p: 4,
            backgroundColor: 'rgba(255,255,255,0.9)',
            borderRadius: 2,
          }}
        >
          <Typography variant="h4" gutterBottom color="text.secondary">
            We want to hear from you!
          </Typography>

          <Typography variant="body1" paragraph color="text.secondary">
            We're thrilled you made it this far, because it only means one thing: you want to get in
            touch with us!
          </Typography>
          <Typography variant="body1" paragraph color="text.secondary">
            Whether you've got a question, a suggestion, a crazy idea, or just want to say hello
            (and hey, we love that too), no one here will judge your curiosity—we're huge fans of
            it.
          </Typography>
          <Typography variant="body1" paragraph color="text.secondary">
            If you're looking for answers, a collaboration, or even just a friendly chat, we're here
            for you.
          </Typography>
          <Typography variant="body1" paragraph color="text.secondary">
            We promise not to hide behind endless forms or leave you on "read." We love conversation
            as much as we love coffee (and that's saying a lot).
          </Typography>
          <Typography variant="body1" paragraph color="text.secondary">
            Choose how you want to reach us and send your message our way.
          </Typography>

          <Typography variant="h4" gutterBottom color="text.secondary" mt={4}>
            Contact Us
          </Typography>

          <Box component="form" onSubmit={handleSubmit} sx={{ mt: 2 }}>
            <Grid container spacing={2}>
              <Grid size={{ xs: 12, sm: 6 }}>
                <TextField
                  name="name"
                  label="Name"
                  value={formData.name}
                  onChange={handleInputChange}
                  fullWidth
                  required
                />
              </Grid>

              <Grid size={{ xs: 12, sm: 6 }}>
                <TextField
                  name="email"
                  label="Email"
                  type="email"
                  value={formData.email}
                  onChange={handleInputChange}
                  fullWidth
                  required
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <FormControl fullWidth required>
                  <InputLabel>Category</InputLabel>
                  <Select
                    name="category"
                    value={formData.category}
                    label="Category"
                    onChange={handleSelectChange}
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
                    <MenuItem value="Oral Implantology">Oral Implantology</MenuItem>
                    <MenuItem value="Pediatric Dentistry">Pediatric Dentistry</MenuItem>
                    <MenuItem value="Prosthodontics / Oral Rehabilitation">
                      Prosthodontics / Oral Rehabilitation
                    </MenuItem>
                    <MenuItem value="Cosmetic Dentistry">Cosmetic Dentistry</MenuItem>
                    <MenuItem value="Oral & Maxillofacial Pathology">
                      Oral & Maxillofacial Pathology
                    </MenuItem>
                    <MenuItem value="Oral Medicine">Oral Medicine</MenuItem>
                    <MenuItem value="Oral & Maxillofacial Radiology">
                      Oral & Maxillofacial Radiology
                    </MenuItem>
                    <MenuItem value="Preventive & Community Dentistry">
                      Preventive & Community Dentistry
                    </MenuItem>
                    <MenuItem value="Forensic Dentistry">Forensic Dentistry</MenuItem>
                    <MenuItem value="Dental Sleep Medicine">Dental Sleep Medicine</MenuItem>
                    <MenuItem value="Geriatric Dentistry">Geriatric Dentistry</MenuItem>
                    <MenuItem value="Restorative Dentistry">Restorative Dentistry</MenuItem>
                    <MenuItem value="Digital Dentistry / CAD-CAM">
                      Digital Dentistry / CAD-CAM
                    </MenuItem>
                    <MenuItem value="Minimally Invasive Dentistry">
                      Minimally Invasive Dentistry
                    </MenuItem>
                    <MenuItem value="Dental Biomaterials & Bioengineering">
                      Dental Biomaterials & Bioengineering
                    </MenuItem>
                    <MenuItem value="Dental Education & Research">
                      Dental Education & Research
                    </MenuItem>
                    <MenuItem value="Practice Management & Administration">
                      Practice Management & Administration
                    </MenuItem>
                    <MenuItem value="Student">Student</MenuItem>
                  </Select>
                </FormControl>
              </Grid>

              <Grid size={{ xs: 12, sm: 6 }}>
                <TextField
                  name="phone"
                  label="Phone"
                  value={formData.phone}
                  onChange={handleInputChange}
                  fullWidth
                  required
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <TextField
                  name="message"
                  label="Message"
                  value={formData.message}
                  onChange={handleInputChange}
                  multiline
                  rows={4}
                  fullWidth
                  required
                />
              </Grid>

              <Grid size={{ xs: 12 }}>
                <Button
                  type="submit"
                  variant="contained"
                  fullWidth
                  sx={{
                    backgroundColor: '#F97316',
                    '&:hover': { backgroundColor: '#EA580C' },
                    py: 1.5,
                  }}
                >
                  Send
                </Button>
              </Grid>
            </Grid>
          </Box>
        </Paper>
      </Container>

      <PoweredBy />
      <Partners />
      <Footer />
    </Box>
  );
};

export default ContactUs;
