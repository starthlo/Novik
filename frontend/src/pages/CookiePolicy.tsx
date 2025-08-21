import { useState } from 'react';
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
  '& ul, & ol': {
    marginBottom: '1rem',
    paddingLeft: '2rem',
    '& li': {
      fontSize: '1rem',
      lineHeight: 1.7,
      marginBottom: '0.5rem',
      color: novikTheme.colors.textMuted,
      fontFamily: novikTheme.typography.fontFamily,
    },
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
  fontFamily: novikTheme.typography.fontFamily,
});

const SectionBox = styled(Box)({
  marginBottom: '2rem',
});

const CookiePolicy = () => {
  const [language, setLanguage] = useState<'en' | 'es'>('en');

  const content = {
    en: {
      title: 'Novik Cookie Policy',
      effectiveDate: 'Effective Date: 17 August 2025',
      intro:
        "This Cookie Policy explains how Novik (operated by BlueBioPlan S.L.) uses cookies and similar technologies on our website and application. We comply with global data-protection laws, including the EU ePrivacy Directive and GDPR, the UK GDPR, the California Consumer Privacy Act/CPRA, Brazil's LGPD, Canada's PIPEDA and Australia's Privacy Act 1988. Because Novik is used globally by healthcare professionals, our cookie practices are designed to meet the most stringent requirements and give you control over your data.",
      sections: [
        {
          title: '1. What Are Cookies?',
          content:
            'Cookies are small text files placed on your computer or device when you visit a website. They allow the site to recognize your device and remember information such as your login status, language preference or pages visited. We also use similar technologies like local storage, SDKs and web beacons; for simplicity, we refer to all of them as "cookies."',
        },
        {
          title: '2. Why We Use Cookies',
          content: `Novik uses cookies for several purposes:

• **Strictly necessary functions** – These cookies are essential to operate the platform and enable core features such as authenticating users, keeping you logged in, enabling secure payment transactions or remembering your cookie preferences. Without these cookies, certain services cannot be provided.

• **Functionality cookies** – With your consent, we use cookies to remember choices you make (e.g. language selection or interface customizations) and provide enhanced functionality.

• **Analytics cookies** – With your consent, we use analytics cookies to understand how you use our site (pages visited, time spent, errors encountered). This helps us improve the performance and content of Novik. For example, we use Google Analytics cookies such as _ga to collect anonymized statistics. All analytics data is aggregated and not used to identify you.

• **Advertising/targeting cookies** – We do not use advertising or targeting cookies. Novik does not serve third-party advertisements or track your online activity across other websites. Should we introduce advertising in the future, we will update this policy and obtain your explicit consent before setting any advertising cookies.

Some cookies are set by us (first-party cookies), while others are set by third parties we use (such as Stripe for payment processing or Google Analytics). We only place non-essential cookies with your consent.`,
        },
        {
          title: '3. Cookie Consent and Control',
          content: `Under the GDPR and the ePrivacy Directive, organisations must obtain users' consent before setting any cookies other than strictly necessary cookies. They must also provide clear information about what data each cookie collects, allow users to refuse non-essential cookies and make withdrawing consent as easy as giving it. We honour these requirements by using a cookie consent banner powered by a consent management platform (such as Cookiebot):

• **Accept all** – When you click "Accept all," you consent to the use of all cookie categories.
• **Reject all** – Rejects all non-essential cookies; only strictly necessary cookies will be set.
• **Customize** – Allows you to selectively enable functionality or analytics cookies while rejecting others.

You can change your preferences at any time by reopening the cookie settings (look for a cookie icon or link on the site) or by deleting cookies in your browser. We document and store your consent as required by law.`,
        },
        {
          title: '4. Legal Basis for Cookies',
          content: `• **Strictly necessary cookies** are used because they are required to deliver the service you request. We rely on our legitimate interests and contractual necessity to place them. In some jurisdictions, they do not require consent.

• **Functionality and analytics cookies** are set only with your consent (Art. 6(1)(a) GDPR). You may withdraw consent at any time, and doing so will not affect the lawfulness of prior processing. In jurisdictions such as California, rejecting these cookies is equivalent to opting out of the "sale" or "sharing" of personal information.

• We do not use advertising cookies or sell or share personal information.`,
        },
        {
          title: '5. Your Rights',
          content:
            'Because cookies may involve the processing of personal data, you have rights that vary depending on your jurisdiction. These include the right to be informed, the right of access, rectification or deletion, the right to object or withdraw consent, the right to portability and (where applicable) the right to opt out of the sale or sharing of personal information. Please refer to our Privacy Policy for a detailed description of your rights and how to exercise them. You can also exercise cookie-specific choices via our consent banner.',
        },
        {
          title: '6. Retention of Cookies',
          content:
            'Cookies remain on your device for varying periods. Session cookies expire when you close your browser, while persistent cookies have an expiration date set by the cookie provider. We endeavour to keep cookie lifespans no longer than necessary for their purpose. You can manually delete cookies at any time through your browser settings.',
        },
        {
          title: '7. Updates to This Policy',
          content:
            'We may update this Cookie Policy to reflect changes in our practices or legal obligations. We will post the updated policy with a new effective date. Significant changes will be notified via our site or, if appropriate, by email. Your continued use of Novik after the update constitutes acceptance of the new policy.',
        },
        {
          title: '8. Contact Us',
          content:
            'If you have questions or comments about this Cookie Policy, please email us at info@novik.ai or write to BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, Spain.',
        },
      ],
    },
    es: {
      title: 'Política de Cookies de Novik',
      effectiveDate: 'Fecha de entrada en vigor: 17 de agosto de 2025',
      intro:
        'Esta Política de Cookies explica cómo Novik (operado por BlueBioPlan S.L.) utiliza cookies y tecnologías similares en nuestro sitio web y aplicación. Cumplimos con las leyes globales de protección de datos, incluidas la Directiva de privacidad electrónica y el RGPD de la UE, el RGPD del Reino Unido, la Ley de Privacidad del Consumidor de California/CPRA, la LGPD de Brasil, la PIPEDA de Canadá y la Ley de Privacidad de 1988 de Australia. Debido a que Novik es utilizado globalmente por profesionales de la salud, nuestras prácticas de cookies están diseñadas para cumplir con los requisitos más estrictos y darle control sobre sus datos.',
      sections: [
        {
          title: '1. ¿Qué son las cookies?',
          content:
            'Las cookies son pequeños archivos de texto que los sitios web colocan en su ordenador o dispositivo cuando los visita. Permiten que el sitio reconozca su dispositivo y recuerde información como su estado de inicio de sesión, preferencia de idioma o páginas visitadas. También utilizamos tecnologías similares como el almacenamiento local, SDKs y balizas web; por simplicidad, nos referimos a todas ellas como «cookies».',
        },
        {
          title: '2. Por qué utilizamos cookies',
          content: `Novik utiliza cookies para varios fines:

• **Funciones estrictamente necesarias** – Estas cookies son esenciales para hacer funcionar la plataforma y permiten funciones básicas como autenticar a los usuarios, mantenerle conectado, habilitar transacciones de pago seguras o recordar sus preferencias de cookies. Sin estas cookies no se pueden prestar determinados servicios.

• **Cookies de funcionalidad** – Con su consentimiento, utilizamos cookies para recordar las opciones que realiza (p. ej., selección de idioma o personalizaciones de interfaz) y proporcionar funciones mejoradas.

• **Cookies de analítica** – Con su consentimiento, utilizamos cookies de analítica para comprender cómo utiliza nuestro sitio (páginas visitadas, tiempo de navegación, errores detectados). Esto nos ayuda a mejorar el rendimiento y el contenido de Novik. Por ejemplo, utilizamos cookies de Google Analytics como _ga para recopilar estadísticas anonimizadas. Todos los datos de analítica se agregan y no se utilizan para identificarle.

• **Cookies de publicidad/segmentación** – No utilizamos cookies de publicidad ni de segmentación. Novik no muestra anuncios de terceros ni rastrea su actividad en otros sitios web. Si en el futuro introducimos publicidad, actualizaremos esta política y obtendremos su consentimiento explícito antes de instalar cookies publicitarias.

Algunas cookies las establecemos nosotros (cookies de primera parte), mientras que otras las establecen terceros que utilizamos (como Stripe para pagos o Google Analytics). Solo colocamos cookies no esenciales con su consentimiento.`,
        },
        {
          title: '3. Consentimiento y control de cookies',
          content: `Según el GDPR y la Directiva de privacidad electrónica, las organizaciones deben obtener el consentimiento de los usuarios antes de colocar cualquier cookie que no sea estrictamente necesaria. También deben proporcionar información clara sobre los datos que recopila cada cookie, permitir a los usuarios rechazar las cookies no esenciales y facilitar la retirada del consentimiento. Cumplimos estos requisitos mediante un banner de consentimiento de cookies proporcionado por una plataforma de gestión de consentimiento (como Cookiebot):

• **Aceptar todo** – Al hacer clic en «Aceptar todo», consiente el uso de todas las categorías de cookies.
• **Rechazar todo** – Rechaza todas las cookies no esenciales; solo se instalarán las cookies estrictamente necesarias.
• **Personalizar** – Le permite habilitar selectivamente cookies de funcionalidad o analítica mientras rechaza otras.

Puede cambiar sus preferencias en cualquier momento volviendo a abrir la configuración de cookies (busque un icono o enlace de cookies en el sitio) o eliminando las cookies en su navegador. Documentamos y almacenamos su consentimiento según lo exige la ley.`,
        },
        {
          title: '4. Base legal de las cookies',
          content: `• **Las cookies estrictamente necesarias** se utilizan porque son necesarias para prestar el servicio que solicita. Nos basamos en nuestros intereses legítimos y en la necesidad contractual para colocarlas. En algunas jurisdicciones no se requiere consentimiento para estas cookies.

• **Las cookies de funcionalidad y analítica** se instalan solo con su consentimiento (art. 6(1)(a) GDPR). Puede retirar el consentimiento en cualquier momento; retirarlo no afecta a la legalidad del tratamiento anterior. En jurisdicciones como California, rechazar estas cookies equivale a optar por no participar en la «venta» o «compartición» de información personal.

• No utilizamos cookies publicitarias ni vendemos o compartimos información personal.`,
        },
        {
          title: '5. Sus derechos',
          content:
            'Como las cookies pueden implicar el tratamiento de datos personales, usted tiene derechos que varían según su jurisdicción. Estos incluyen el derecho a ser informado, el derecho de acceso, rectificación o eliminación, el derecho a oponerse o retirar el consentimiento, el derecho a la portabilidad y (cuando corresponda) el derecho a optar por no vender o compartir información personal. Consulte nuestra Política de Privacidad para obtener una descripción detallada de sus derechos y cómo ejercerlos. También puede ejercer opciones específicas de cookies mediante nuestro banner de consentimiento.',
        },
        {
          title: '6. Conservación de las cookies',
          content:
            'Las cookies permanecen en su dispositivo durante distintos períodos. Las cookies de sesión expiran cuando cierra el navegador, mientras que las cookies persistentes tienen una fecha de expiración establecida por el proveedor. Procuramos que la duración de las cookies no sea superior a la necesaria para su propósito. Puede eliminar las cookies manualmente en cualquier momento a través de la configuración de su navegador.',
        },
        {
          title: '7. Actualizaciones de esta política',
          content:
            'Podemos actualizar esta Política de Cookies para reflejar cambios en nuestras prácticas o obligaciones legales. Publicaremos la política actualizada con una nueva fecha de entrada en vigor. Los cambios significativos se notificarán a través de nuestro sitio o, en su caso, por correo electrónico. Su uso continuado de Novik después de la actualización constituye aceptación de la nueva política.',
        },
        {
          title: '8. Contacto',
          content:
            'Si tiene preguntas o comentarios sobre esta Política de Cookies, envíenos un correo electrónico a info@novik.ai o escríbanos a BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, España.',
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

      <EffectiveDate align="center">{currentContent.effectiveDate}</EffectiveDate>

      <SectionBox>
        <Typography paragraph>{currentContent.intro}</Typography>
      </SectionBox>

      <Divider sx={{ my: 3 }} />

      {currentContent.sections.map((section, index) => (
        <SectionBox key={index}>
          <Typography variant="h2" component="h2">
            {section.title}
          </Typography>
          {section.content.split('\n\n').map((paragraph, pIndex) => {
            // Check if paragraph contains bullet points
            if (paragraph.includes('•')) {
              const items = paragraph.split('\n').filter(item => item.trim().startsWith('•'));
              const introText = paragraph.split('\n')[0];

              return (
                <Box key={pIndex}>
                  {introText && !introText.startsWith('•') && (
                    <Typography paragraph>{introText}</Typography>
                  )}
                  <ul>
                    {items.map((item, iIndex) => {
                      const text = item.replace('•', '').trim();
                      // Parse bold text
                      const parts = text.split(/(\*\*.*?\*\*)/);
                      return (
                        <li key={iIndex}>
                          {parts.map((part, partIndex) => {
                            if (part.startsWith('**') && part.endsWith('**')) {
                              return <strong key={partIndex}>{part.slice(2, -2)}</strong>;
                            }
                            return part;
                          })}
                        </li>
                      );
                    })}
                  </ul>
                </Box>
              );
            } else {
              // Parse bold text in regular paragraphs
              const parts = paragraph.split(/(\*\*.*?\*\*)/);
              return (
                <Typography paragraph key={pIndex}>
                  {parts.map((part, partIndex) => {
                    if (part.startsWith('**') && part.endsWith('**')) {
                      return <strong key={partIndex}>{part.slice(2, -2)}</strong>;
                    }
                    return part;
                  })}
                </Typography>
              );
            }
          })}
        </SectionBox>
      ))}
    </PageContainer>
  );
};

export default CookiePolicy;
