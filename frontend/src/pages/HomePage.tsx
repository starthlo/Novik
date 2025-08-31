import {
  Box,
  Container,
  Typography,
  Button,
  styled,
  Card,
  List,
  ListItem,
  Grid,
} from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';
import { novikTheme } from '../styles/theme';
import CookieConsent from '../components/Common/CookieConsent';
import Footer from '../components/Common/Footer';
import FrontImage from '../assets/Front Image.png';
import NovikLogoWhite from '../assets/novik-logo-white.png';
import PubMedWhite from '../assets/pubmed-white.png';
import DrugBankWhite from '../assets/drugbank-white.png';

const HeroSection = styled(Box)({
  minHeight: '90vh',
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  position: 'relative',
  textAlign: 'center',
  color: '#ffffff',
  backgroundSize: 'cover',
  backgroundImage: `url(${FrontImage})`,
  backgroundPosition: 'center',
  backgroundColor: '#777777',
  backgroundBlendMode: 'multiply',
  '&::before': {
    content: '""',
    position: 'absolute',
    inset: 0,
    background: 'rgba(50, 50, 50, 0.4)',
    backdropFilter: 'blur(2px)',
  },
});

const HeroContent = styled(Box)({
  position: 'relative',
  maxWidth: '800px',
  zIndex: 1,
});

const NovikLogo = styled('img')({
  height: '280px',
  width: 'auto',
  filter: 'brightness(1.1)',
});

const PartnersSection = styled(Box)({
  width: '420px',
  marginTop: '100px',
  zIndex: 1,
  position: 'relative',
});

const PartnerLogo = styled('img')({
  height: '150px',
  width: 'auto',
  filter: 'brightness(1.1)',
});

const PartnerLogoSecond = styled('img')({
  position: 'absolute',
  right: '0',
  bottom: '-53px',
  height: '240px',
  width: 'auto',
  filter: 'brightness(1.1)',
});

const HeroTitle = styled(Typography)({
  fontSize: '2.5rem',
  marginBottom: '1rem',
  fontWeight: 700,
  fontFamily: novikTheme.typography.fontFamily,
  '@media (max-width: 600px)': {
    fontSize: '1.8rem',
  },
});

const HeroSubtitle = styled(Typography)({
  fontSize: '1.2rem',
  marginBottom: '2rem',
  color: '#eaeaea',
  fontFamily: novikTheme.typography.fontFamily,
});

const CTAButton = styled(Button)<any>({
  padding: '0.75rem 1.6rem',
  borderRadius: '30px',
  fontWeight: 600,
  fontSize: '1rem',
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  transition: 'background 0.3s ease',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const Section = styled(Box)<{ variant?: 'default' | 'alt' }>(({ variant }) => ({
  padding: '4rem 1rem',
  backgroundColor: variant === 'alt' ? novikTheme.colors.section : novikTheme.colors.background,
}));

const SectionTitle = styled(Box)({
  textAlign: 'center',
  marginBottom: '2rem',
});

const SectionHeading = styled(Typography)({
  fontSize: '2rem',
  fontWeight: 700,
  marginBottom: '0.5rem',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
});

const SectionSubheading = styled(Typography)({
  maxWidth: '700px',
  margin: '0.5rem auto',
  color: novikTheme.colors.textMuted,
  fontFamily: novikTheme.typography.fontFamily,
  lineHeight: 1.6,
});

const FeatureCard = styled(Card)({
  background: '#fff',
  borderRadius: '20px',
  boxShadow: '0 2px 5px rgba(0,0,0,0.05)',
  padding: '1.5rem',
  height: '100%',
  transition: 'transform 0.2s ease, box-shadow 0.2s ease',
  '&:hover': {
    transform: 'translateY(-5px)',
    boxShadow: '0 4px 12px rgba(0,0,0,0.1)',
  },
});

const ProblemCard = styled(FeatureCard)({
  textAlign: 'center',
  display: 'flex',
  flexDirection: 'column',
  justifyContent: 'center',
  alignItems: 'center',
});

const CardNumber = styled(Box)<any>({
  display: 'inline-flex',
  alignItems: 'center',
  justifyContent: 'center',
  backgroundColor: novikTheme.colors.primary,
  color: '#fff',
  borderRadius: '50%',
  width: '36px',
  height: '36px',
  fontWeight: 700,
  marginRight: '0.5rem',
});

const StyledList = styled(List)({
  listStyle: 'none',
  maxWidth: '800px',
  margin: '0 auto',
  padding: 0,
});

const StyledListItem = styled(ListItem)({
  display: 'flex',
  alignItems: 'flex-start',
  marginBottom: '1rem',
  paddingLeft: 0,
  color: novikTheme.colors.textMuted,
  lineHeight: 1.65,
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

const FAQItem = styled('details')({
  backgroundColor: '#fff',
  borderRadius: '8px',
  padding: '1rem 1.5rem',
  marginBottom: '1rem',
  boxShadow: '0 1px 3px rgba(0,0,0,0.08)',
  '& summary': {
    cursor: 'pointer',
    fontWeight: 600,
    listStyle: 'none',
    position: 'relative',
    '&::-webkit-details-marker': {
      display: 'none',
    },
    '&::after': {
      content: '"▼"',
      position: 'absolute',
      right: 0,
      top: 0,
      transition: 'transform 0.3s',
    },
  },
  '&[open] summary::after': {
    transform: 'rotate(180deg)',
  },
  '& p': {
    marginTop: '0.5rem',
    color: novikTheme.colors.textMuted,
  },
});

const HomePage = () => {
  const problemsData = [
    {
      icon: (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" style={{ width: '60px' }}>
          <circle cx="50" cy="50" r="48" fill={novikTheme.colors.primary} />
          <text x="50%" y="55%" textAnchor="middle" fill="#ffffff" fontSize="28" fontWeight="700">
            70%
          </text>
        </svg>
      ),
      title: 'Polypharmacy prevalence',
      description:
        'More than 70% of dental patients have at least one underlying condition or take multiple medications, increasing the chance of interactions.',
    },
    {
      icon: (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" style={{ width: '60px' }}>
          <circle
            cx="50"
            cy="50"
            r="45"
            fill="none"
            stroke={novikTheme.colors.primary}
            strokeWidth="6"
          />
          <line
            x1="50"
            y1="50"
            x2="50"
            y2="25"
            stroke={novikTheme.colors.primary}
            strokeWidth="6"
            strokeLinecap="round"
          />
          <line
            x1="50"
            y1="50"
            x2="70"
            y2="60"
            stroke={novikTheme.colors.primary}
            strokeWidth="6"
            strokeLinecap="round"
          />
        </svg>
      ),
      title: 'Time constraints',
      description:
        'Dentists rarely have the bandwidth to review drug interactions or keep up with evolving protocols before each procedure.',
    },
    {
      icon: (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" style={{ width: '60px' }}>
          <polygon
            points="50,10 90,85 10,85"
            fill="none"
            stroke={novikTheme.colors.primary}
            strokeWidth="6"
          />
          <line
            x1="50"
            y1="35"
            x2="50"
            y2="60"
            stroke={novikTheme.colors.primary}
            strokeWidth="6"
            strokeLinecap="round"
          />
          <circle cx="50" cy="75" r="4" fill={novikTheme.colors.primary} />
        </svg>
      ),
      title: 'Risk of clinical errors',
      description:
        'Human oversights can lead to ethical, legal and health complications – especially when managing high-risk patients.',
    },
  ];

  const howItWorksData = [
    {
      number: 1,
      title: 'Input clinical data',
      description:
        'Provide patient history, current medications and your intended dental procedure. Upload optional PDF reports for richer context.',
    },
    {
      number: 2,
      title: 'AI analysis & risk detection',
      description:
        'Our engine screens for contraindications, calculates adjusted dosages and flags potential hazards based on evidence.',
    },
    {
      number: 3,
      title: 'Get AI support & suggestions',
      description:
        'Receive tailored suggestions, safe alternatives, adjusted dosages and curated references to support your clinical decision-making.',
    },
  ];

  const featuresData = [
    {
      number: 1,
      title: 'Risk triage',
      description:
        'Automatically classifies each case as low, moderate or high complexity, routing it through the appropriate AI model for optimal recommendations.',
    },
    {
      number: 2,
      title: 'Evidence engine',
      description:
        'Queries medical literature from the last five years to deliver 3–7 citations in Vancouver format, complete with PubMed links.',
    },
    {
      number: 3,
      title: 'Drug intelligence',
      description:
        "Suggests antibiotics, analgesics and local anaesthetics with dose limits, alternatives and interactions based on each patient's history.",
    },
    {
      number: 4,
      title: 'Data ingestion',
      description:
        'Upload PDFs of medical histories and lab results. Novik extracts key findings to enrich its recommendations.',
    },
    {
      number: 5,
      title: 'Multilingual',
      description:
        'Works in multiple languages and automatically responds in the language you use.',
    },
    {
      number: 6,
      title: 'Traceability',
      description:
        'Every recommendation comes with its justification and source, so you know exactly why a suggestion was made.',
    },
  ];

  const evidencePoints = [
    {
      title: 'Methodology',
      description:
        'We filter PubMed by condition, drug and dental procedure, prioritising systematic reviews, clinical guidelines and meta-analyses.',
    },
    {
      title: 'Scope',
      description:
        'Covers antibiotics, NSAIDs, local anaesthetics, anticoagulants, bisphosphonates, diabetes, hypertension and more.',
    },
    {
      title: 'Formatting',
      description:
        'References are delivered in Vancouver format with direct links to abstracts and DrugBank details.',
    },
    {
      title: 'Transparency',
      description:
        'Each data source and update timestamp is clearly shown, along with known limitations and biases.',
    },
  ];

  const faqItems = [
    {
      question: 'What is Novik?',
      answer:
        'Novik is an AI-powered assistant designed exclusively for dental professionals. It synthesises patient data and medical literature to support treatment planning without replacing your clinical expertise.',
    },
    {
      question: 'Is Novik really free?',
      answer:
        'Yes. During our open beta, all core features are free for licensed clinicians. We are funded by sponsorships and grants, not by selling your data.',
    },
    {
      question: 'Is my patient data secure?',
      answer:
        "Absolutely. We comply with the EU GDPR, encrypt data in transit and at rest, and only collect anonymised usage metrics. You remain in control of your patients' information.",
    },
    {
      question: "Does Novik replace a dentist's judgement?",
      answer:
        'No. Novik provides decision support to help you make informed choices, but the final judgement always resides with the treating clinician.',
    },
    {
      question: 'Which languages does Novik support?',
      answer:
        'Novik responds in the same language you use. It is multilingual and automatically replies in your language.',
    },
  ];

  const aboutPoints = [
    {
      title: 'Mission',
      description:
        'Empower dental professionals with trustworthy AI to deliver safer, evidence-based care.',
    },
    {
      title: "What we're not",
      description:
        'Novik is a decision support tool and never a substitute for clinical judgement.',
    },
    {
      title: 'Team',
      description:
        'A collective of dentists, AI researchers and software engineers united by a passion for improving oral health.',
    },
    {
      title: 'Milestones',
      description:
        'Finalist for Best Innovation at Barcelona Dental Show 2025 and winner of the Best Technological Advance at the 28th Gaceta Dental Awards in Spain.',
    },
  ];

  return (
    <>
      <CookieConsent />

      {/* Hero Section */}
      <HeroSection>
        <NovikLogo src={NovikLogoWhite} alt="Novik" />
        <HeroContent>
          <Box>
            <HeroTitle variant="h1">Smarter decisions for faster, safer dentistry</HeroTitle>
            <HeroSubtitle>
              From diagnosis to treatment, Novik guides your clinical decisions with precision.
            </HeroSubtitle>
            <CTAButton component={RouterLink} to="/contact">
              Start working with Novik for free
            </CTAButton>
          </Box>
        </HeroContent>
        <PartnersSection>
          <PartnerLogo src={PubMedWhite} alt="PubMed" />
          <PartnerLogoSecond src={DrugBankWhite} alt="DrugBank" />
        </PartnersSection>
      </HeroSection>

      {/* Problem Section */}
      <Section variant="alt">
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">The problem we solve</SectionHeading>
            <SectionSubheading>
              Complex medical histories, polypharmacy and limited time make treatment planning risky
              and stressful for dental professionals. Novik simplifies this process.
            </SectionSubheading>
          </SectionTitle>
          <Grid container spacing={3}>
            {problemsData.map((problem, index) => (
              <Grid size={{ xs: 12, md: 4 }} key={index}>
                <ProblemCard>
                  <Box sx={{ mb: 2 }}>{problem.icon}</Box>
                  <Typography variant="h6" fontWeight={600} gutterBottom>
                    {problem.title}
                  </Typography>
                  <Typography variant="body2" color="text.secondary">
                    {problem.description}
                  </Typography>
                </ProblemCard>
              </Grid>
            ))}
          </Grid>
        </Container>
      </Section>

      {/* How It Works Section */}
      <Section>
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">How it works</SectionHeading>
            <SectionSubheading>
              Novik streamlines decision-making in three simple steps.
            </SectionSubheading>
          </SectionTitle>
          <Grid container spacing={3}>
            {howItWorksData.map((step, index) => (
              <Grid size={{ xs: 12, md: 4 }} key={index}>
                <FeatureCard>
                  <Typography variant="h6" fontWeight={600} gutterBottom>
                    <CardNumber component="span">{step.number}</CardNumber>
                    {step.title}
                  </Typography>
                  <Typography variant="body2" color="text.secondary">
                    {step.description}
                  </Typography>
                </FeatureCard>
              </Grid>
            ))}
          </Grid>
        </Container>
      </Section>

      {/* Evidence Section */}
      <Section variant="alt">
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">Evidence at your fingertips</SectionHeading>
            <SectionSubheading>
              Transparent methodology ensures recommendations are grounded in high-quality research.
            </SectionSubheading>
          </SectionTitle>
          <StyledList>
            {evidencePoints.map((point, index) => (
              <StyledListItem key={index}>
                <span className="list-content">
                  <strong>{point.title}.</strong> {point.description}
                </span>
              </StyledListItem>
            ))}
          </StyledList>
        </Container>
      </Section>

      {/* Product Features Section */}
      <Section>
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">What Novik does</SectionHeading>
            <SectionSubheading>
              A digital co-pilot that turns complex data into clear, actionable treatment plans.
            </SectionSubheading>
          </SectionTitle>
          <Grid container spacing={3}>
            {featuresData.map((feature, index) => (
              <Grid size={{ xs: 12, md: 4 }} key={index}>
                <FeatureCard>
                  <Typography variant="h6" fontWeight={600} gutterBottom>
                    <CardNumber component="span">{feature.number}</CardNumber>
                    {feature.title}
                  </Typography>
                  <Typography variant="body2" color="text.secondary">
                    {feature.description}
                  </Typography>
                </FeatureCard>
              </Grid>
            ))}
          </Grid>
        </Container>
      </Section>

      {/* FAQ Section */}
      <Section variant="alt">
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">Frequently Asked Questions</SectionHeading>
            <SectionSubheading>Answers to common questions about Novik.</SectionSubheading>
          </SectionTitle>
          <Box sx={{ maxWidth: '800px', margin: '0 auto' }}>
            {faqItems.map((faq, index) => (
              <FAQItem key={index}>
                <summary>{faq.question}</summary>
                <p>{faq.answer}</p>
              </FAQItem>
            ))}
          </Box>
          <Box sx={{ textAlign: 'center', mt: 3 }}>
            <CTAButton component={RouterLink} to="/faqs">
              More FAQs
            </CTAButton>
          </Box>
        </Container>
      </Section>

      {/* About Section */}
      <Section>
        <Container maxWidth="lg">
          <SectionTitle>
            <SectionHeading variant="h2">About Novik</SectionHeading>
            <SectionSubheading>Our mission, team and values.</SectionSubheading>
          </SectionTitle>
          <StyledList>
            {aboutPoints.map((point, index) => (
              <StyledListItem key={index}>
                <span className="list-content">
                  <strong>{point.title}.</strong> {point.description}
                </span>
              </StyledListItem>
            ))}
          </StyledList>
        </Container>
      </Section>

      <Footer />
    </>
  );
};

export default HomePage;
