import { Container, Box, Typography, styled, Chip, List, ListItem } from '@mui/material';
import { novikTheme } from '../styles/theme';

// Styled components following the HTML mockup design
const PageContainer = styled(Box)({
  fontFamily: novikTheme.typography.fontFamily,
  padding: '7rem 1rem 4rem',
  minHeight: '100vh',
  backgroundColor: novikTheme.colors.background,
});

const ContentContainer = styled(Container)({
  maxWidth: '900px !important',
  textAlign: 'center',
});

const ComingSoonBadge = styled(Chip)({
  marginBottom: '1rem',
  padding: '0.5rem 1rem',
  backgroundColor: novikTheme.colors.primaryDark,
  color: '#fff',
  fontWeight: 500,
  fontSize: '0.9rem',
});

const PageTitle = styled(Typography)({
  fontSize: '2.2rem',
  fontWeight: 700,
  marginBottom: '1rem',
  color: novikTheme.colors.text,
});

const IntroText = styled(Typography)({
  marginBottom: '2rem',
  lineHeight: 1.6,
  textAlign: 'left',
  color: novikTheme.colors.textMuted,
  '& strong': {
    color: novikTheme.colors.text,
    fontWeight: 600,
  },
});

const SectionTitle = styled(Typography)({
  fontSize: '1.5rem',
  fontWeight: 600,
  marginTop: '2rem',
  marginBottom: '1rem',
  color: novikTheme.colors.primaryDark,
});

const StyledList = styled(List)({
  listStyle: 'none',
  paddingLeft: 0,
  textAlign: 'left',
  margin: '0 0 2rem 0',
});

const StyledListItem = styled(ListItem)({
  display: 'flex',
  alignItems: 'flex-start',
  marginBottom: '1rem',
  paddingLeft: 0,
  color: novikTheme.colors.textMuted,
  lineHeight: 1.6,
  fontSize: '0.95rem',
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

const ContactBox = styled(Box)({
  marginTop: '3rem',
  padding: '2rem',
  backgroundColor: novikTheme.colors.section,
  borderRadius: novikTheme.borderRadius.medium,
  border: `1px solid ${novikTheme.colors.borderLight}`,
});

const ApiPage = () => {
  return (
    <PageContainer>
      <ContentContainer>
        <ComingSoonBadge label="Coming soon" />

        <PageTitle variant="h1">API & Developers</PageTitle>

        <IntroText>
          The Novik API is designed to seamlessly integrate with your existing digital ecosystem.
          While clinicians benefit from Novik through the platform, the API opens powerful
          possibilities for <strong>software providers and enterprise partners</strong>. If you are
          developing management platforms, EHR systems, or other health‑tech solutions, Novik's API
          allows you to embed advanced dental clinical decision support into your own product.
        </IntroText>

        <SectionTitle variant="h2">What You Can Do</SectionTitle>
        <StyledList>
          <StyledListItem>
            <span className="list-content">
              <strong>Integration with Third‑Party Software:</strong> Incorporate Novik's
              intelligence into dental management platforms, EHR systems, or enterprise solutions to
              add value to your product offering.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Decision Support Integration:</strong> Automate evidence‑based suggestions for
              medications, anesthetics, and antibiotic prophylaxis, surfaced directly in your
              software tools.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Workflow Automation:</strong> Enable bulk submission of patient histories and
              receive structured risk assessments and treatment guidelines programmatically.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Custom Analytics:</strong> Provide your clients with dashboards that compare
              clinic performance, monitor compliance, and track prescribing trends.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Research Applications:</strong> Pull curated PubMed references related to
              patient data to support academic studies and scientific publishing.
            </span>
          </StyledListItem>
        </StyledList>

        <SectionTitle variant="h2">Authentication</SectionTitle>
        <StyledList>
          <StyledListItem>
            <span className="list-content">
              Secure token‑based authentication to ensure only authorised clients can access the
              API.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              OAuth 2.0 workflows for user‑level access when needed.
            </span>
          </StyledListItem>
        </StyledList>

        <SectionTitle variant="h2">Endpoints</SectionTitle>
        <StyledList>
          <StyledListItem>
            <span className="list-content">
              <strong>Submit case:</strong> Send patient history, medications and procedure data to
              receive a personalised decision support plan.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Retrieve recommendations:</strong> Access past recommendations for a case,
              including dose adjustments and bibliographies.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              <strong>Status webhooks:</strong> Receive asynchronous updates when long‑running
              analyses are complete.
            </span>
          </StyledListItem>
        </StyledList>

        <SectionTitle variant="h2">SDKs</SectionTitle>
        <StyledList>
          <StyledListItem>
            <span className="list-content">
              Client libraries will be available in JavaScript, Python and other popular languages
              to simplify integration.
            </span>
          </StyledListItem>
        </StyledList>

        <SectionTitle variant="h2">Status & Changelog</SectionTitle>
        <StyledList>
          <StyledListItem>
            <span className="list-content">
              A public status page will provide real‑time updates on API availability and
              performance.
            </span>
          </StyledListItem>
          <StyledListItem>
            <span className="list-content">
              Versioned endpoints and detailed changelogs will ensure backwards compatibility.
            </span>
          </StyledListItem>
        </StyledList>

        <ContactBox>
          <Typography
            variant="h6"
            sx={{ mb: 2, fontWeight: 600, color: novikTheme.colors.primaryDark }}
          >
            Get Early Access
          </Typography>
          <Typography sx={{ color: novikTheme.colors.textMuted, lineHeight: 1.6 }}>
            For API access inquiries, contact us at{' '}
            <a
              href="mailto:info@novik.ai"
              style={{
                color: novikTheme.colors.primary,
                textDecoration: 'none',
                fontWeight: 600,
              }}
            >
              info@novik.ai
            </a>
          </Typography>
          <Typography sx={{ mt: 2, color: novikTheme.colors.textMuted, fontSize: '0.9rem' }}>
            We're currently onboarding select partners for our beta API program.
          </Typography>
        </ContactBox>
      </ContentContainer>
    </PageContainer>
  );
};

export default ApiPage;
