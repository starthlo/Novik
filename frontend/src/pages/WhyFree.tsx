import { Container, Box, Typography, styled, List, ListItem } from '@mui/material';
import { novikTheme } from '../styles/theme';
import { useEffect } from 'react';

const PageContainer = styled(Box)({
  fontFamily: novikTheme.typography.fontFamily,
  padding: '6.5rem 1rem 4rem',
  minHeight: '100vh',
  backgroundColor: novikTheme.colors.background,
});

const ContentContainer = styled(Container)({
  maxWidth: '900px !important',
});

const PageTitle = styled(Typography)({
  fontSize: '2rem',
  fontWeight: 700,
  textAlign: 'center',
  marginBottom: '1rem',
  color: novikTheme.colors.text,
});

const IntroText = styled(Typography)({
  lineHeight: 1.65,
  marginBottom: '2rem',
  color: novikTheme.colors.textMuted,
  fontSize: '1rem',
});

const StyledList = styled(List)({
  listStyle: 'none',
  paddingLeft: 0,
  margin: 0,
});

const StyledListItem = styled(ListItem)({
  display: 'flex',
  alignItems: 'flex-start',
  marginBottom: '1.5rem',
  paddingLeft: 0,
  color: novikTheme.colors.textMuted,
  lineHeight: 1.65,
  fontSize: '1rem',
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

const InfoBox = styled(Box)({
  marginTop: '3rem',
  padding: '2rem',
  backgroundColor: novikTheme.colors.section,
  borderRadius: novikTheme.borderRadius.medium,
  border: `1px solid ${novikTheme.colors.borderLight}`,
});

const WhyFree = () => {
  useEffect(() => {
    window.scrollTo({ top: 0, behavior: 'instant' });
  }, []);

  const reasons = [
    {
      title: 'Open beta',
      description:
        'Novik is in a public beta phase. We offer free access to validate our models with practising clinicians and gather feedback.',
    },
    {
      title: 'Sponsorship & grants',
      description:
        'Our work is funded by partners and grants. We do not sell personal data or run intrusive advertising.',
    },
    {
      title: 'Anonymous metrics',
      description:
        'We collect aggregated usage data to improve the system; individual patient data remains private and secure.',
    },
    {
      title: 'Transparent future',
      description:
        'When we introduce paid plans, the free tier will remain for core functionality. Early users will keep their benefits.',
    },
  ];

  return (
    <PageContainer>
      <ContentContainer>
        <PageTitle variant="h1">Why Novik is free</PageTitle>

        <IntroText>
          We believe that every dental professional deserves access to accurate, evidence‑based
          decision support. Here's why you can use Novik at no cost during our open beta:
        </IntroText>

        <StyledList>
          {reasons.map((reason, index) => (
            <StyledListItem key={index}>
              <span className="list-content">
                <strong>{reason.title}.</strong> {reason.description}
              </span>
            </StyledListItem>
          ))}
        </StyledList>

        <InfoBox>
          <Typography
            variant="h6"
            sx={{
              mb: 2,
              fontWeight: 600,
              color: novikTheme.colors.primaryDark,
            }}
          >
            Our Commitment
          </Typography>
          <Typography sx={{ color: novikTheme.colors.textMuted, lineHeight: 1.65 }}>
            We're committed to maintaining a free tier that provides essential functionality to all
            dental professionals. Our mission is to democratize access to AI-powered clinical
            decision support, ensuring that quality dental care isn't limited by financial
            constraints.
          </Typography>
          <Typography sx={{ mt: 2, color: novikTheme.colors.textMuted, lineHeight: 1.65 }}>
            As we grow, premium features may be introduced for advanced use cases, but the core
            Novik experience will always remain accessible to those who need it most.
          </Typography>
        </InfoBox>

        <Box sx={{ mt: 4, textAlign: 'center' }}>
          <Typography sx={{ color: novikTheme.colors.textMuted }}>
            Questions about our pricing model?{' '}
            <a
              href="/contact"
              style={{
                color: novikTheme.colors.primary,
                textDecoration: 'none',
                fontWeight: 600,
              }}
            >
              Contact us
            </a>
          </Typography>
        </Box>
      </ContentContainer>
    </PageContainer>
  );
};

export default WhyFree;
