import { useState } from 'react';
import {
  Container,
  Box,
  Typography,
  TextField,
  Button,
  styled,
  List,
  ListItem,
} from '@mui/material';
import { novikTheme } from '../styles/theme';

const PageContainer = styled(Box)({
  fontFamily: novikTheme.typography.fontFamily,
  padding: '7rem 1rem 4rem',
  minHeight: '100vh',
  backgroundColor: novikTheme.colors.background,
});

const ContentContainer = styled(Container)({
  maxWidth: '900px !important',
});

const PageTitle = styled(Typography)({
  fontSize: '2rem',
  fontWeight: 700,
  marginBottom: '1rem',
  textAlign: 'center',
  color: novikTheme.colors.text,
});

const IntroText = styled(Typography)({
  marginBottom: '1rem',
  lineHeight: 1.6,
  color: novikTheme.colors.textMuted,
});

const BenefitsList = styled(List)({
  listStyle: 'none',
  paddingLeft: 0,
  marginBottom: '2rem',
});

const BenefitItem = styled(ListItem)({
  display: 'flex',
  alignItems: 'flex-start',
  marginBottom: '1rem',
  paddingLeft: 0,
  color: novikTheme.colors.textMuted,
  lineHeight: 1.6,
  '&::before': {
    content: '""',
    minWidth: '10px',
    height: '10px',
    backgroundColor: novikTheme.colors.primary,
    borderRadius: '50%',
    marginRight: '15px',
    marginTop: '8px',
    flexShrink: 0,
  },
  '& .list-content': {
    flex: 1,
  },
  '& strong': {
    color: novikTheme.colors.text,
    fontWeight: 600,
  },
});

const FormContainer = styled(Box)({
  marginTop: '2rem',
  backgroundColor: novikTheme.colors.section,
  padding: '2rem',
  borderRadius: novikTheme.borderRadius.medium,
});

const StyledTextField = styled(TextField)({
  marginBottom: '1rem',
  '& .MuiInputBase-root': {
    borderRadius: '6px',
    fontFamily: novikTheme.typography.fontFamily,
    backgroundColor: '#ffffff',
  },
  '& .MuiInputLabel-root': {
    fontFamily: novikTheme.typography.fontFamily,
    fontWeight: 500,
  },
  '& .MuiOutlinedInput-root': {
    backgroundColor: '#ffffff',
    '& fieldset': {
      borderColor: '#ddd',
    },
    '&:hover fieldset': {
      borderColor: novikTheme.colors.primary,
    },
    '&.Mui-focused fieldset': {
      borderColor: novikTheme.colors.primary,
    },
  },
});

const SubmitButton = styled(Button)({
  backgroundColor: novikTheme.colors.primary,
  color: '#fff',
  padding: '0.6rem 1.4rem',
  borderRadius: '20px',
  fontWeight: 600,
  fontSize: '1rem',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  transition: 'background 0.2s ease',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const HighlightBox = styled(Box)({
  backgroundColor: '#fff',
  border: `2px solid ${novikTheme.colors.primary}`,
  borderRadius: novikTheme.borderRadius.medium,
  padding: '1.5rem',
  marginBottom: '2rem',
  marginTop: '2rem',
});

interface FormData {
  name: string;
  company: string;
  email: string;
  phone: string;
  message: string;
}

const Partners = () => {
  const [formData, setFormData] = useState<FormData>({
    name: '',
    company: '',
    email: '',
    phone: '',
    message: '',
  });

  const [isSubmitting, setIsSubmitting] = useState(false);

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);

    // Simulate form submission
    setTimeout(() => {
      console.log('Partnership inquiry submitted:', formData);
      alert('Thank you for your interest in partnering with Novik! We will contact you soon.');
      setFormData({
        name: '',
        company: '',
        email: '',
        phone: '',
        message: '',
      });
      setIsSubmitting(false);
    }, 1000);
  };

  const benefits = [
    {
      title: 'Targeted exposure',
      description:
        'Connect with thousands of practising dentists and dental students who rely on Novik for clinical decisions.',
    },
    {
      title: 'High visibility placements',
      description:
        'Enjoy banner placements, sponsored webinars and custom content within our ecosystem.',
    },
    {
      title: 'Trusted environment',
      description: 'Associate your brand with a tool that prioritises evidence and patient safety.',
    },
    {
      title: 'Data‑driven reporting',
      description: 'Receive detailed metrics on impressions, clicks and engagement.',
    },
  ];

  return (
    <PageContainer>
      <ContentContainer>
        <PageTitle variant="h1">Partner with Novik</PageTitle>

        <IntroText>
          As the first AI decision‑support platform tailored to dentistry, Novik reaches a global
          community of forward‑thinking clinicians. By sponsoring Novik, you align your brand with
          innovation, precision and safety in oral healthcare.
        </IntroText>

        <BenefitsList>
          {benefits.map((benefit, index) => (
            <BenefitItem key={index}>
              <span className="list-content">
                <strong>{benefit.title}.</strong> {benefit.description}
              </span>
            </BenefitItem>
          ))}
        </BenefitsList>

        <HighlightBox>
          <Typography
            variant="h6"
            sx={{
              mb: 1,
              fontWeight: 600,
              color: novikTheme.colors.primaryDark,
              textAlign: 'center',
            }}
          >
            Current Partners & Achievements
          </Typography>
          <Typography sx={{ color: novikTheme.colors.textMuted, textAlign: 'center', mb: 1 }}>
            🏆 Finalist for Best Innovation at Barcelona Dental Show 2025
          </Typography>
          <Typography sx={{ color: novikTheme.colors.textMuted, textAlign: 'center' }}>
            🥇 Winner of Best Technological Advance at the 28th Gaceta Dental Awards
          </Typography>
        </HighlightBox>

        <IntroText>
          If you're interested in partnering with us, please fill out the form below and our team
          will be in touch.
        </IntroText>

        <FormContainer>
          <Box component="form" onSubmit={handleSubmit}>
            <StyledTextField
              fullWidth
              label="Name"
              name="name"
              value={formData.name}
              onChange={handleInputChange}
              required
              variant="outlined"
            />

            <StyledTextField
              fullWidth
              label="Company"
              name="company"
              value={formData.company}
              onChange={handleInputChange}
              variant="outlined"
            />

            <StyledTextField
              fullWidth
              label="Email"
              name="email"
              type="email"
              value={formData.email}
              onChange={handleInputChange}
              required
              variant="outlined"
            />

            <StyledTextField
              fullWidth
              label="Phone"
              name="phone"
              type="tel"
              value={formData.phone}
              onChange={handleInputChange}
              variant="outlined"
            />

            <StyledTextField
              fullWidth
              label="Message"
              name="message"
              value={formData.message}
              onChange={handleInputChange}
              required
              multiline
              rows={5}
              variant="outlined"
            />

            <SubmitButton type="submit" disabled={isSubmitting}>
              {isSubmitting ? 'Sending...' : 'Send message'}
            </SubmitButton>
          </Box>
        </FormContainer>
      </ContentContainer>
    </PageContainer>
  );
};

export default Partners;
