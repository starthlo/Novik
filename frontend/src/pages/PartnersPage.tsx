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

interface FormData {
  name: string;
  email: string;
  category: string;
  phone: string;
  message: string;
}

const PartnersPage = () => {
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
          Hello, potential partners and allies!
        </Typography>
        <Typography variant="body1" color="text.secondary">
          We're thrilled to have you here because it means something about Novik has already caught
          your attention.
        </Typography>
        <Typography variant="body1" color="text.secondary">
          And let me tell you the excitement is mutual! We're always on the lookout for connections
          that feel more like partnerships...
        </Typography>
        <Typography variant="body1" color="text.secondary">
          So go ahead, if you think Novik could be a great place for your support, we're just a
          message away.
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
                  <MenuItem value="Practice Management and Administration">
                    Practice Management and Administration
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
  );
};

export default PartnersPage;
