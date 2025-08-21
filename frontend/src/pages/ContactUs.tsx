import { useState } from 'react';
import { Container, Box, Typography, TextField, Button, styled } from '@mui/material';
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
  fontSize: '2.2rem',
  fontWeight: 700,
  marginBottom: '1.2rem',
  textAlign: 'center',
  color: novikTheme.colors.primaryDark,
});

const IntroText = styled(Typography)({
  marginBottom: '1rem',
  lineHeight: 1.6,
  color: novikTheme.colors.textMuted,
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
  padding: '0.7rem 1.6rem',
  borderRadius: '25px',
  fontWeight: 600,
  fontSize: '1rem',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  transition: 'background 0.3s ease',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

interface FormData {
  name: string;
  email: string;
  subject: string;
  message: string;
}

const ContactUs = () => {
  const [formData, setFormData] = useState<FormData>({
    name: '',
    email: '',
    subject: '',
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
      console.log('Form submitted:', formData);
      alert('Thank you for your message! We will get back to you soon.');
      setFormData({
        name: '',
        email: '',
        subject: '',
        message: '',
      });
      setIsSubmitting(false);
    }, 1000);
  };

  return (
    <PageContainer>
      <ContentContainer>
        <PageTitle variant="h1">Contact us</PageTitle>

        <IntroText>
          We would be delighted to hear from you. At Novik, we value meaningful connections and
          believe that every message is an opportunity to build something together. Whether you want
          to share your impressions, suggest improvements, explore collaboration opportunities or
          simply ask a question, our team is here to listen attentively and respond with care.
        </IntroText>

        <IntroText>
          Your insights and feedback help us grow, refine our platform and continue pushing the
          boundaries of innovation in dentistry. We aim to cultivate long-term relationships based
          on trust, professionalism and transparency, and we welcome conversations that help us
          achieve that vision.
        </IntroText>

        <IntroText>
          Please use the form below to reach out. Every enquiry will be reviewed with the attention
          it deserves, and we will make sure you receive a thoughtful reply from our team.
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
              label="Subject"
              name="subject"
              value={formData.subject}
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

export default ContactUs;
