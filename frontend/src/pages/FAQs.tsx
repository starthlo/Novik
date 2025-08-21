import { useState } from 'react';
import { Container, Box, Typography, styled } from '@mui/material';
import { novikTheme } from '../styles/theme';

const PageContainer = styled(Box)({
  fontFamily: novikTheme.typography.fontFamily,
  paddingTop: '2rem',
  paddingBottom: '4rem',
  minHeight: '100vh',
  backgroundColor: novikTheme.colors.background,
});

const ContentContainer = styled(Container)({
  maxWidth: '960px !important',
});

const PageTitle = styled(Typography)({
  textAlign: 'center',
  fontSize: '2rem',
  fontWeight: 700,
  marginBottom: '1.25rem',
  marginTop: '1.25rem',
  padding: '2rem',
  borderRadius: '12px',
  background: `linear-gradient(135deg, ${novikTheme.colors.primary} 0%, ${novikTheme.colors.primaryDark} 100%)`,
  color: '#fff',
});

const SectionTitle = styled(Typography)({
  fontSize: '1.25rem',
  fontWeight: 600,
  margin: '2rem 0 1rem',
  color: novikTheme.colors.primaryDark,
});

const FAQItem = styled('details')({
  backgroundColor: '#fff',
  borderRadius: novikTheme.borderRadius.medium,
  padding: '1rem 1.25rem',
  marginBottom: '0.9rem',
  boxShadow: novikTheme.shadows.small,
  transition: 'box-shadow 0.3s ease',
  '&:hover': {
    boxShadow: novikTheme.shadows.medium,
  },
});

const FAQSummary = styled('summary')({
  cursor: 'pointer',
  fontWeight: 600,
  listStyle: 'none',
  position: 'relative',
  fontSize: '1rem',
  color: novikTheme.colors.text,
  '&::-webkit-details-marker': {
    display: 'none',
  },
  '&::after': {
    content: '"▾"',
    position: 'absolute',
    right: 0,
    top: 0,
    transition: 'transform 0.3s',
    color: novikTheme.colors.primary,
  },
});

const FAQAnswer = styled(Typography)({
  marginTop: '0.6rem',
  color: novikTheme.colors.textMuted,
  lineHeight: 1.65,
  fontSize: '0.95rem',
});

const ContactSection = styled(Box)({
  marginTop: '3rem',
  textAlign: 'center',
  padding: '2rem',
  backgroundColor: novikTheme.colors.section,
  borderRadius: novikTheme.borderRadius.medium,
});

interface FAQItemData {
  question: string;
  answer: string;
}

interface FAQSectionData {
  title: string;
  items: FAQItemData[];
}

const FAQs = () => {
  const [openItems, setOpenItems] = useState<Set<string>>(new Set());

  const handleToggle = (id: string) => {
    const newOpenItems = new Set(openItems);
    if (newOpenItems.has(id)) {
      newOpenItems.delete(id);
    } else {
      newOpenItems.add(id);
    }
    setOpenItems(newOpenItems);
  };

  const faqSections: FAQSectionData[] = [
    {
      title: 'General',
      items: [
        {
          question: 'What is Novik?',
          answer:
            "Novik is an AI‑powered clinical decision support tool built for dentistry. It synthesises patient data with medical literature to support safe treatment planning. It never replaces a clinician's judgement.",
        },
        {
          question: 'Who can use Novik?',
          answer:
            'Licensed dentists and dental care professionals, plus dental students in supervised settings. It is not intended for laypersons.',
        },
        {
          question: 'Where does Novik work?',
          answer: 'Novik is web‑based and runs in a modern browser. No installation is required.',
        },
        {
          question: 'Can Novik be used outside dentistry?',
          answer:
            'No. Novik is designed specifically for dentistry. However, its methodology may inspire future tools in other medical specialties.',
        },
      ],
    },
    {
      title: 'Usage & Features',
      items: [
        {
          question: 'What data can I input?',
          answer:
            'Age, weight, medical history, current medications, allergies and planned dental procedure. Optional PDF uploads let you add lab reports or previous notes.',
        },
        {
          question: 'What procedures and drugs are covered?',
          answer:
            'Common procedures (prophylaxis, fillings, extractions, endodontics, implants) and core drug classes (antibiotics, analgesics/NSAIDs, local anaesthetics). Anticoagulants, bisphosphonates and chronic conditions are also considered.',
        },
        {
          question: 'How does triage work?',
          answer:
            'Each case is classified by complexity (low / moderate / high). Higher‑risk cases trigger more rigorous checks and the most capable reasoning models.',
        },
        {
          question: 'Does Novik diagnose?',
          answer:
            'No. Novik provides evidence‑based suggestions and dose adjustments. Diagnosis and final decisions always remain with the treating clinician.',
        },
        {
          question: 'Can I upload images or radiographs?',
          answer:
            'Currently Novik supports PDF uploads. Image and radiograph support is planned for a future release.',
        },
        {
          question: 'Does Novik work offline?',
          answer:
            'No. Novik requires an internet connection to access its models and up‑to‑date medical references.',
        },
      ],
    },
    {
      title: 'Evidence & Methodology',
      items: [
        {
          question: 'Where do the references come from?',
          answer:
            'From articles indexed in PubMed and other internationally recognised scientific societies and publications. Citations are returned in Vancouver style with links to abstracts where available.',
        },
        {
          question: 'How current is the evidence?',
          answer:
            'Recent literature is prioritised (typically the last 5 years) while preserving seminal sources when relevant. Each recommendation includes an update timestamp.',
        },
        {
          question: 'Does Novik include drug interactions?',
          answer:
            'Yes. Novik checks for drug–drug and drug–condition interactions based on the information you provide and links to DrugBank profiles.',
        },
        {
          question: 'Can I trust Novik over guidelines?',
          answer:
            'Novik aligns with clinical guidelines but is not a substitute. Always prioritise official protocols in your jurisdiction.',
        },
      ],
    },
    {
      title: 'Languages & Compatibility',
      items: [
        {
          question: 'What languages does Novik support?',
          answer:
            'Novik is multilingual and replies in the same language you use (for example: English, Spanish, Portuguese, Italian, French, German — expanding over time).',
        },
        {
          question: 'Can I integrate Novik with my software?',
          answer:
            'Yes. Our API enables integration with practice management systems and custom tools. See the API page for details.',
        },
        {
          question: 'Which browsers are supported?',
          answer:
            'Latest versions of Chrome, Safari, Firefox and Edge are recommended. Mobile browsers are supported but desktop use provides the best experience.',
        },
      ],
    },
    {
      title: 'Security & Compliance',
      items: [
        {
          question: 'Is my data secure?',
          answer:
            'Data is encrypted in transit and at rest. Access is restricted to authorised users, and we follow EU GDPR requirements. We only collect anonymised usage metrics to improve the service.',
        },
        {
          question: 'Where is data stored?',
          answer: 'On secure EU‑based servers with regular backups and access logging.',
        },
        {
          question: 'Do you share data with third parties?',
          answer:
            'No personal clinical data is shared with third parties. Aggregated anonymous usage statistics may be used for service improvements.',
        },
        {
          question: 'Is Novik HIPAA compliant?',
          answer:
            'Novik complies with GDPR and is exploring HIPAA alignment for use in the United States.',
        },
      ],
    },
    {
      title: 'Plans & Pricing',
      items: [
        {
          question: 'Is Novik free?',
          answer:
            'Yes. During open beta the core features are free. Future paid tiers may add higher limits and premium capabilities while maintaining a free option.',
        },
        {
          question: 'Will my account change when pricing launches?',
          answer:
            'Early users will keep free access to core functionality. Any changes will be communicated well in advance.',
        },
        {
          question: 'Will there be an enterprise version?',
          answer:
            'Yes. We are planning enterprise tiers with analytics, team dashboards and advanced integrations for large clinics or networks.',
        },
      ],
    },
    {
      title: 'Account & Access',
      items: [
        {
          question: 'How do I sign up?',
          answer: 'Use the "Login/Register" button in the header to create your account.',
        },
        {
          question: 'I forgot my password. What now?',
          answer:
            'Use the password reset link on the login page. If issues persist, contact us at info@novik.ai.',
        },
        {
          question: 'Can I delete my account?',
          answer:
            'Yes. Contact support and request permanent deletion of your account and associated data.',
        },
        {
          question: 'Can I use Novik on multiple devices?',
          answer: 'Yes. You can log in from any device with a modern browser.',
        },
      ],
    },
    {
      title: 'Support',
      items: [
        {
          question: 'How can I send feedback or request a feature?',
          answer:
            'Use the Contact page form or email info@novik.ai. We review all clinical and usability feedback during the beta.',
        },
        {
          question: 'Do you offer training?',
          answer:
            'No. Novik does not provide training. However, the interface is intuitive and supported by clear documentation and FAQs.',
        },
        {
          question: 'How quickly do you respond to support requests?',
          answer:
            'We aim to respond to all inquiries within 2 business days during the beta phase.',
        },
      ],
    },
  ];

  return (
    <PageContainer>
      <ContentContainer>
        <PageTitle variant="h1">Frequently Asked Questions</PageTitle>

        {faqSections.map((section, sectionIndex) => (
          <Box key={sectionIndex} sx={{ mb: 3 }}>
            <SectionTitle variant="h2">{section.title}</SectionTitle>

            {section.items.map((item, itemIndex) => {
              const itemId = `${sectionIndex}-${itemIndex}`;
              return (
                <FAQItem
                  key={itemIndex}
                  open={openItems.has(itemId)}
                  onToggle={() => handleToggle(itemId)}
                  style={{
                    cursor: 'pointer',
                  }}
                >
                  <FAQSummary
                    style={{
                      transform: openItems.has(itemId) ? 'none' : 'none',
                    }}
                  >
                    {item.question}
                  </FAQSummary>
                  <FAQAnswer>{item.answer}</FAQAnswer>
                </FAQItem>
              );
            })}
          </Box>
        ))}

        <ContactSection>
          <Typography variant="h6" sx={{ mb: 2, fontWeight: 600, color: novikTheme.colors.text }}>
            Can't find what you're looking for?
          </Typography>
          <Typography sx={{ color: novikTheme.colors.textMuted }}>
            Contact us at{' '}
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
        </ContactSection>
      </ContentContainer>
    </PageContainer>
  );
};

export default FAQs;
