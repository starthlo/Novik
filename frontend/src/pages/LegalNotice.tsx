import { useState, useEffect } from 'react';
import { Box, Container, Typography, Button, styled, Divider } from '@mui/material';
import { novikTheme } from '../styles/theme';

const PageContainer = styled(Container)({
  paddingTop: '2rem',
  paddingBottom: '4rem',
  maxWidth: '900px',
  '& h1': {
    fontSize: '2.5rem',
    fontWeight: 700,
    marginBottom: '1rem',
    color: novikTheme.colors.text,
    fontFamily: novikTheme.typography.fontFamily,
  },
  '& h2': {
    fontSize: '1.8rem',
    fontWeight: 600,
    marginTop: '2.5rem',
    marginBottom: '1rem',
    color: novikTheme.colors.text,
    fontFamily: novikTheme.typography.fontFamily,
  },
  '& h3': {
    fontSize: '1.3rem',
    fontWeight: 600,
    marginTop: '1.5rem',
    marginBottom: '0.8rem',
    color: novikTheme.colors.text,
    fontFamily: novikTheme.typography.fontFamily,
  },
  '& p': {
    fontSize: '1rem',
    lineHeight: 1.7,
    marginBottom: '1rem',
    color: novikTheme.colors.textMuted,
    fontFamily: novikTheme.typography.fontFamily,
  },
  '& strong': {
    fontWeight: 600,
    color: novikTheme.colors.text,
  },
});

const LanguageToggle = styled(Box)({
  display: 'flex',
  justifyContent: 'center',
  gap: '1rem',
  marginBottom: '2rem',
  paddingTop: '1rem',
});

const LanguageButton = styled(Button)<{ active?: boolean }>(({ active }) => ({
  padding: '0.5rem 1.5rem',
  borderRadius: '20px',
  fontWeight: 600,
  fontSize: '0.9rem',
  textTransform: 'none',
  backgroundColor: active ? novikTheme.colors.primary : 'transparent',
  color: active ? '#ffffff' : novikTheme.colors.text,
  border: `2px solid ${novikTheme.colors.primary}`,
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: active ? novikTheme.colors.primaryDark : 'rgba(136, 169, 78, 0.1)',
  },
}));

const EffectiveDate = styled(Typography)({
  fontSize: '0.9rem',
  color: novikTheme.colors.textMuted,
  fontStyle: 'italic',
  marginBottom: '1.5rem',
  textAlign: 'center',
  fontFamily: novikTheme.typography.fontFamily,
});

const SubTitle = styled(Typography)({
  fontSize: '1.2rem',
  fontWeight: 600,
  marginTop: '2rem',
  marginBottom: '1rem',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
});

const InfoBlock = styled(Box)({
  backgroundColor: novikTheme.colors.section,
  padding: '1.5rem',
  borderRadius: '8px',
  marginBottom: '1.5rem',
  marginTop: '1rem',
  border: `1px solid ${novikTheme.colors.border}`,
});

const ContactInfo = styled(Box)({
  '& p': {
    margin: '0.5rem 0',
    color: novikTheme.colors.textMuted,
    fontFamily: novikTheme.typography.fontFamily,
  },
  '& strong': {
    color: novikTheme.colors.text,
    fontWeight: 600,
  },
});

const LegalNotice = () => {
  const [language, setLanguage] = useState<'en' | 'es'>('en');

  useEffect(() => {
    window.scrollTo({ top: 0, behavior: 'instant' });
  }, []);

  const content = {
    en: {
      title: 'Legal Notice',
      documentType: 'Legal Documentation',
      effectiveDate: '17 August 2025',
      sections: [
        {
          title: 'Legal Notice',
          content:
            'In accordance with the duty of information set forth in Spanish Law 34/2002 on Information Society Services and Electronic Commerce (LSSI-CE), the following details are provided:',
        },
        {
          title: 'Company Information',
          type: 'info',
          details: [
            { label: 'Website Owner', value: 'BlueBioPlan S.L.' },
            { label: 'VAT / Tax ID', value: 'B-72913379' },
            {
              label: 'Registered Address',
              value: 'C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, Spain',
            },
            { label: 'Contact Email', value: 'info@novik.ai' },
          ],
        },
        {
          title: 'Purpose of the Website',
          content:
            'This website provides access to Novik, an AI-powered clinical decision support tool for dental professionals.',
        },
        {
          title: 'Intellectual and Industrial Property',
          content:
            'All content, trademarks, logos, designs, text, and software on this site are the property of BlueBioPlan S.L. or third parties licensed for use. Any reproduction, distribution, or modification without express consent is prohibited.',
        },
        {
          title: 'Disclaimer of Liability',
          content:
            'The information and recommendations provided by Novik are intended as clinical support and do not replace the professional judgment of dentists. BlueBioPlan S.L. is not liable for damages arising from misuse of the platform.',
        },
        {
          title: 'Jurisdiction and Applicable Law',
          content:
            'Any disputes or claims related to this website shall be governed by Spanish law and submitted to the courts of Murcia, Spain.',
        },
      ],
    },
    es: {
      title: 'Aviso Legal',
      documentType: 'Documentación Legal',
      effectiveDate: '17 de agosto de 2025',
      sections: [
        {
          title: 'Aviso Legal',
          content:
            'En cumplimiento con el deber de información establecido en la Ley 34/2002 de Servicios de la Sociedad de la Información y del Comercio Electrónico (LSSI-CE), se facilitan los siguientes datos identificativos:',
        },
        {
          title: 'Información de la Empresa',
          type: 'info',
          details: [
            { label: 'Titular del sitio web', value: 'BlueBioPlan S.L.' },
            { label: 'NIF/CIF', value: 'B-72913379' },
            {
              label: 'Domicilio social',
              value: 'C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, España',
            },
            { label: 'Correo electrónico de contacto', value: 'info@novik.ai' },
          ],
        },
        {
          title: 'Objeto del sitio web',
          content:
            'El presente sitio web tiene por finalidad dar a conocer la plataforma Novik, un asistente clínico con inteligencia artificial para profesionales de la odontología.',
        },
        {
          title: 'Propiedad intelectual e industrial',
          content:
            'Todos los contenidos, marcas, logotipos, diseños, textos y software presentes en este sitio web son propiedad de BlueBioPlan S.L. o de terceros con autorización. Queda prohibida su reproducción, distribución o modificación sin consentimiento expreso.',
        },
        {
          title: 'Exclusión de responsabilidad',
          content:
            'La información y recomendaciones ofrecidas por Novik tienen carácter de apoyo clínico y no sustituyen el juicio profesional del odontólogo. BlueBioPlan S.L. no se hace responsable de los daños derivados del uso indebido de la plataforma.',
        },
        {
          title: 'Jurisdicción y legislación aplicable',
          content:
            'Para la resolución de cualquier conflicto o cuestión relacionada con este sitio web, será de aplicación la legislación española y se someterán a los juzgados y tribunales de Murcia (España).',
        },
      ],
    },
  };

  const currentContent = content[language];

  return (
    <PageContainer>
      <LanguageToggle>
        <LanguageButton active={language === 'en'} onClick={() => setLanguage('en')}>
          English
        </LanguageButton>
        <LanguageButton active={language === 'es'} onClick={() => setLanguage('es')}>
          Español
        </LanguageButton>
      </LanguageToggle>

      <Typography variant="h1" component="h1" align="center">
        {currentContent.title}
      </Typography>

      <EffectiveDate>
        {currentContent.documentType}
        <br />
        {currentContent.effectiveDate}
      </EffectiveDate>

      <Divider sx={{ my: 3 }} />

      {currentContent.sections.map((section, index) => {
        if (section.type === 'info' && section.details) {
          return (
            <Box key={index} sx={{ mb: 3 }}>
              <InfoBlock>
                <ContactInfo>
                  {section.details.map((detail, dIndex) => (
                    <Typography key={dIndex} paragraph sx={{ mb: 1 }}>
                      <strong>{detail.label}:</strong> {detail.value}
                    </Typography>
                  ))}
                </ContactInfo>
              </InfoBlock>
            </Box>
          );
        }

        return (
          <Box key={index} sx={{ mb: 3 }}>
            <SubTitle>{section.title}</SubTitle>
            <Typography paragraph>{section.content}</Typography>
          </Box>
        );
      })}

      <Box sx={{ mt: 5, pt: 3, borderTop: `1px solid ${novikTheme.colors.border}` }}>
        <Typography
          sx={{
            textAlign: 'center',
            fontSize: '0.9rem',
            fontStyle: 'italic',
            color: 'text.secondary',
          }}
        >
          BlueBioPlan S.L. – Novik Legal Documentation
        </Typography>
      </Box>
    </PageContainer>
  );
};

export default LegalNotice;
