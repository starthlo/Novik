import { useEffect, useState } from 'react';
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
  textAlign: 'center',
  fontFamily: novikTheme.typography.fontFamily,
});

const SectionBox = styled(Box)({
  marginBottom: '2rem',
});

const TermsOfUse = () => {
  const [language, setLanguage] = useState<'en' | 'es'>('en');

  useEffect(() => {
    window.scrollTo({ top: 0, behavior: 'instant' });
  }, []);

  const content = {
    en: {
      title: 'Novik Terms of Use',
      effectiveDate: 'Last Updated: 17 August 2025',
      intro:
        'These Terms of Use (the "Terms") govern access to and use of Novik, an AI-powered clinical decision support platform operated by BlueBioPlan S.L. ("BlueBioPlan", "we", "our" or "us"). By accessing or using Novik you agree to be bound by these Terms. If you do not agree, do not use Novik.',
      sections: [
        {
          title: '1. Introduction and Acceptance',
          content:
            'Novik provides AI-generated suggestions to assist licensed dental and medical professionals. It is intended solely for use by adults who are legally qualified healthcare professionals. Novik is not designed for patients, minors or the general public. By using Novik you represent that you are a licensed professional (or supervised trainee) of at least 18 years of age and that you have authority to enter into these Terms on your own behalf or on behalf of an organisation you represent. You also agree to our Privacy Policy and Cookie Policy, which form part of these Terms.',
        },
        {
          title: '2. Description of the Service',
          content:
            'Novik analyses anonymized clinical inputs and external sources to provide informational suggestions. The platform does not provide medical or dental advice, diagnoses or treatment. All outputs are for informational and educational purposes only. You remain solely responsible for verifying any information from Novik and for all professional decisions. Novik should never be used in emergencies or as the sole basis for diagnosing or treating any individual. Novik does not create a provider–patient relationship between BlueBioPlan and any person.',
        },
        {
          title: '3. Eligibility, Accounts and Security',
          content:
            'To use Novik you must be a licensed dental or medical professional (or a supervised student/trainee) and of legal adult age in your jurisdiction. You are responsible for maintaining the confidentiality of your login credentials and for all activity under your account. You must not allow patients, minors or other unauthorised individuals to access Novik. You may not transfer your account without our prior written consent. If you believe your account has been compromised, notify us immediately.',
        },
        {
          title: '4. Data Input and De-identification',
          content:
            'Novik is designed to process de-identified case data only. Health information that has been de-identified is not protected health information (PHI) under HIPAA. Do not input any patient identifiers or personal information. You are solely responsible for removing names, contact details, medical record numbers, facial images and other identifiers from the data you submit. If identifiable data is inadvertently provided, we will treat it as accidental and delete it. We do not collect or maintain PHI and are not a covered entity or business associate under HIPAA.',
        },
        {
          title: '5. Acceptable Use and Prohibited Activities',
          content: `You agree to use Novik for lawful purposes and in accordance with these Terms. In particular, you must not:

• **Use Novik for diagnosis or treatment.** Do not use Novik to diagnose or treat any specific patient or to provide direct medical advice.

• **Rely without verification.** Do not rely solely on Novik's outputs; always apply your professional judgment and consult authoritative sources.

• **Submit personal or sensitive data.** Do not enter any personal, identifying or sensitive information about patients or other individuals.

• **Reverse engineer or copy.** Do not attempt to access or discover Novik's source code, algorithms or underlying components, and do not use Novik outputs to train or develop competing products.

• **Resell or misappropriate outputs.** Do not resell, sublicense, publish or otherwise exploit Novik or its outputs outside your own internal professional use without our written consent.

• **Abuse the service.** Do not overload Novik with excessive requests, use bots or scripts, or attempt to disrupt or bypass security controls.

• **Submit illegal or unethical content.** Do not use Novik to generate or transmit content that is illegal, harmful, defamatory, discriminatory, obscene or otherwise unethical.

• **Violate laws.** Comply with all applicable laws, including data-protection, privacy, intellectual property and export control laws.

Violation of this section may result in suspension or termination of your access to Novik and may expose you to legal consequences.`,
        },
        {
          title: '6. User Content and Intellectual Property',
          content: `You retain ownership of the original content or data you provide ("User Content"). You grant BlueBioPlan a worldwide, royalty-free, transferable and sublicensable licence to use, reproduce, modify and create derivative works from your User Content as needed to operate and improve Novik. You represent that you have the rights and consents necessary to submit the User Content.

BlueBioPlan and its licensors own all rights in the Novik software, models, documentation and outputs. We grant you a limited, non-exclusive, non-transferable licence to use Novik outputs solely for your own internal clinical decision support. You may not reproduce, publish, distribute or commercialise Novik outputs without our consent. All Novik and BlueBioPlan trademarks and logos belong to us; you may not use them without permission.

If you provide feedback or suggestions about Novik, we may use them without obligation to you.`,
        },
        {
          title: '7. Suspension and Termination',
          content:
            'We may suspend or terminate your access to Novik if you violate these Terms, if required by law, or if we believe your use poses a risk to the platform, other users or patients. Upon termination, your right to use Novik ceases immediately, but the provisions relating to intellectual property, limitation of liability, indemnification and dispute resolution survive.',
        },
        {
          title: '8. Disclaimers',
          content: `Novik and all associated information are provided "as is" and "as available". To the fullest extent permitted by law, BlueBioPlan disclaims all warranties, express or implied, including warranties of merchantability, fitness for a particular purpose, accuracy and non-infringement. We do not warrant that Novik will meet your requirements, be uninterrupted or error-free, or provide accurate or up-to-date outputs. Your use of Novik is at your own risk. BlueBioPlan makes no guarantee regarding clinical outcomes and assumes no duty of care toward your patients.

Novik may link to third-party services or content. BlueBioPlan does not endorse or control those third parties and is not responsible for their content or practices.`,
        },
        {
          title: '9. Limitation of Liability',
          content:
            'To the extent permitted by law, BlueBioPlan and its affiliates, officers, directors, employees and agents will not be liable for any indirect, incidental, special, consequential or punitive damages (including loss of profits, data or goodwill) arising out of or relating to your use of Novik or these Terms, even if advised of the possibility of such damages. We are not liable for any injuries or damages to patients or third parties resulting from your reliance on Novik. Our total liability under these Terms will not exceed the amount you paid (if any) for access to Novik in the twelve months before the event giving rise to liability; if you use Novik for free, our liability will be zero.',
        },
        {
          title: '10. Indemnification',
          content:
            "You agree to indemnify, defend and hold harmless BlueBioPlan, its affiliates and their respective officers, directors, employees and agents from and against any claims, liabilities, damages, losses and expenses (including reasonable attorneys' fees) arising from: (a) your use or misuse of Novik; (b) your violation of these Terms; or (c) your violation of any law or the rights of any third party. BlueBioPlan may assume the exclusive defence and control of any matter subject to indemnification, in which case you agree to cooperate and not settle any matter without our prior written consent.",
        },
        {
          title: '11. Governing Law and Dispute Resolution',
          content: `These Terms are governed by the laws of Spain, excluding its conflict of law rules. The courts of Spain, particularly those located in Murcia, will have exclusive jurisdiction over any disputes arising from these Terms or your use of Novik. You agree to submit to the personal jurisdiction of those courts.

Before commencing formal legal proceedings, you agree to attempt to resolve any dispute with us informally by contacting info@novik.ai. If a dispute cannot be resolved informally, either party may pursue remedies in the courts described above. Notwithstanding the foregoing, BlueBioPlan may seek injunctive or other equitable relief in any jurisdiction to protect its intellectual property or prevent misuse of Novik.`,
        },
        {
          title: '12. Changes to These Terms',
          content:
            'We may modify these Terms at any time. We will post the updated Terms on Novik and update the "Last Updated" date. For material changes, we will provide additional notice (e.g. via email or in-app notification). Unless required otherwise by law, changes will take effect 15 days after notice. Your continued use of Novik after the effective date constitutes acceptance of the updated Terms. If you do not agree to the changes, you must stop using Novik.',
        },
        {
          title: '13. Miscellaneous',
          content: `These Terms constitute the entire agreement between you and BlueBioPlan regarding Novik and supersede all prior agreements. If any provision is held invalid or unenforceable, the remaining provisions remain in full force. Our failure to enforce any right does not waive it. We may assign these Terms without restriction; you may not assign your rights without our consent. These Terms do not create any partnership, joint venture or agency relationship between us. Headings are for convenience only.

If you have questions about these Terms, please contact info@novik.ai or write to BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, Spain.`,
        },
      ],
    },
    es: {
      title: 'Términos de Uso de Novik',
      effectiveDate: 'Fecha de última actualización: 17 de agosto de 2025',
      intro:
        'Estos Términos de Uso (los "Términos") rigen el acceso y uso de Novik, una plataforma de apoyo a la decisión clínica impulsada por IA operada por BlueBioPlan S.L. ("BlueBioPlan", "nosotros", "nuestro" o "nos"). Al acceder o utilizar Novik, usted acepta estar sujeto a estos Términos. Si no está de acuerdo, no utilice Novik.',
      sections: [
        {
          title: '1. Introducción y aceptación',
          content:
            'Novik ofrece sugerencias generadas por IA para ayudar a profesionales médicos y odontológicos colegiados. Está destinado exclusivamente a adultos que sean profesionales sanitarios legalmente cualificados. Novik no está diseñado para pacientes, menores ni público general. Al utilizar Novik, usted declara que es un profesional colegiado (o un estudiante/tesista supervisado) mayor de 18 años y que tiene autoridad para aceptar estos Términos en su propio nombre o en nombre de la organización que representa. También acepta nuestra Política de Privacidad y Política de Cookies, que forman parte de estos Términos.',
        },
        {
          title: '2. Descripción del servicio',
          content:
            'Novik analiza entradas clínicas anonimizadas y fuentes externas para proporcionar sugerencias informativas. La plataforma no proporciona asesoramiento médico ni odontológico, diagnósticos ni tratamientos. Todos los resultados son únicamente informativos y educativos. Usted sigue siendo el único responsable de verificar cualquier información de Novik y de todas las decisiones profesionales. Novik nunca debe utilizarse en emergencias ni como única base para diagnosticar o tratar a ninguna persona. Novik no crea una relación médico-paciente entre BlueBioPlan y ninguna persona.',
        },
        {
          title: '3. Elegibilidad, cuentas y seguridad',
          content:
            'Para utilizar Novik debe ser un profesional médico u odontológico colegiado (o un estudiante/tesista supervisado) y tener la mayoría de edad legal en su jurisdicción. Usted es responsable de mantener la confidencialidad de sus credenciales de inicio de sesión y de todas las actividades realizadas con su cuenta. No debe permitir que pacientes, menores u otras personas no autorizadas accedan a Novik. No puede transferir su cuenta sin nuestro consentimiento previo por escrito. Si cree que su cuenta se ha visto comprometida, notifíquenoslo de inmediato.',
        },
        {
          title: '4. Introducción de datos y desidentificación',
          content:
            'Novik está diseñado para procesar únicamente datos desidentificados. La información sanitaria que ha sido desidentificada no es información sanitaria protegida (PHI) según la normativa HIPAA. No introduzca identificadores de pacientes ni información personal. Usted es el único responsable de eliminar nombres, datos de contacto, números de historia clínica, imágenes de rostros y otros identificadores de los datos que envíe. Si se proporciona accidentalmente información identificable, la trataremos como accidental y la eliminaremos. No recopilamos ni mantenemos PHI y no somos una entidad cubierta ni un asociado comercial según HIPAA.',
        },
        {
          title: '5. Uso aceptable y actividades prohibidas',
          content: `Usted se compromete a utilizar Novik para fines lícitos y de acuerdo con estos Términos. En particular, no debe:

• **Usar Novik para diagnóstico o tratamiento.** No utilice Novik para diagnosticar o tratar a ningún paciente específico ni para proporcionar asesoramiento médico directo.

• **Confiar sin verificación.** No confíe únicamente en los resultados de Novik; aplique siempre su criterio profesional y consulte fuentes autorizadas.

• **Enviar datos personales o sensibles.** No introduzca información personal, identificativa o sensible sobre pacientes u otras personas.

• **Ingeniería inversa o copia.** No intente acceder o descubrir el código fuente, algoritmos o componentes subyacentes de Novik, ni utilice los resultados de Novik para entrenar o desarrollar productos competidores.

• **Revender o apropiarse indebidamente de los resultados.** No revenda, sublicencie, publique ni explote de otro modo Novik o sus resultados fuera de su propio uso profesional interno sin nuestro consentimiento escrito.

• **Abusar del servicio.** No sobrecargue Novik con solicitudes excesivas, utilice bots o scripts, ni intente interrumpir o sortear los controles de seguridad.

• **Enviar contenido ilegal o poco ético.** No utilice Novik para generar o transmitir contenido ilegal, dañino, difamatorio, discriminatorio, obsceno o poco ético.

• **Violar leyes.** Cumpla todas las leyes aplicables, incluidas las leyes de protección de datos, privacidad, propiedad intelectual y control de exportaciones.

La violación de esta sección puede dar lugar a la suspensión o terminación de su acceso a Novik y puede exponerle a consecuencias legales.`,
        },
        {
          title: '6. Contenido del usuario y propiedad intelectual',
          content: `Usted conserva la propiedad del contenido o datos originales que proporcione («Contenido del Usuario»). Otorga a BlueBioPlan una licencia mundial, gratuita, transferible y sublicenciable para utilizar, reproducir, modificar y crear obras derivadas de su Contenido del Usuario según sea necesario para operar y mejorar Novik. Usted declara que dispone de los derechos y consentimientos necesarios para enviar el Contenido del Usuario.

BlueBioPlan y sus licenciantes poseen todos los derechos sobre el software, los modelos, la documentación y los resultados de Novik. Le concedemos una licencia limitada, no exclusiva e intransferible para utilizar los resultados de Novik únicamente para su propio apoyo interno a la decisión clínica. No puede reproducir, publicar, distribuir ni comercializar los resultados de Novik sin nuestro consentimiento. Todas las marcas y logotipos de Novik y BlueBioPlan nos pertenecen; no puede utilizarlos sin permiso.

Si proporciona comentarios o sugerencias sobre Novik, podemos utilizarlos sin obligación hacia usted.`,
        },
        {
          title: '7. Suspensión y terminación',
          content:
            'Podemos suspender o cancelar su acceso a Novik si infringe estos Términos, si la ley lo exige o si creemos que su uso supone un riesgo para la plataforma, otros usuarios o pacientes. Tras la terminación, su derecho a utilizar Novik cesará inmediatamente, pero las disposiciones relativas a propiedad intelectual, limitación de responsabilidad, indemnización y resolución de disputas seguirán vigentes.',
        },
        {
          title: '8. Descargos de responsabilidad',
          content: `Novik y toda la información asociada se proporcionan «tal cual» y «según disponibilidad». En la máxima medida permitida por la ley, BlueBioPlan rechaza todas las garantías, expresas o implícitas, incluidas las garantías de comerciabilidad, idoneidad para un propósito particular, exactitud y no infracción. No garantizamos que Novik satisfaga sus requisitos, sea ininterrumpido o libre de errores, o proporcione resultados precisos o actualizados. Usted utiliza Novik bajo su propio riesgo. BlueBioPlan no garantiza resultados clínicos y no asume ningún deber de cuidado hacia sus pacientes.

Novik puede enlazar con servicios o contenido de terceros. BlueBioPlan no respalda ni controla a esos terceros y no es responsable de su contenido o prácticas.`,
        },
        {
          title: '9. Limitación de responsabilidad',
          content:
            'En la medida permitida por la ley, BlueBioPlan y sus afiliados, directivos, empleados y agentes no serán responsables de ningún daño indirecto, incidental, especial, consecuente o punitivo (incluidos pérdida de beneficios, datos o reputación) que surja de su uso de Novik o de estos Términos, incluso si se nos ha advertido de la posibilidad de tales daños. No somos responsables de ninguna lesión o daño a pacientes o terceros resultante de su confianza en Novik. Nuestra responsabilidad total en virtud de estos Términos no excederá la cantidad que usted haya pagado (si corresponde) por acceder a Novik en los doce meses anteriores al hecho que originó la responsabilidad; si utiliza Novik de forma gratuita, nuestra responsabilidad será cero.',
        },
        {
          title: '10. Indemnización',
          content:
            'Usted acepta indemnizar, defender y eximir de responsabilidad a BlueBioPlan, sus afiliados y sus respectivos directivos, empleados y agentes por cualquier reclamación, responsabilidad, daño, pérdida y gasto (incluidos honorarios razonables de abogados) que surjan de: (a) su uso o uso indebido de Novik; (b) su violación de estos Términos; o (c) su violación de cualquier ley o de los derechos de terceros. BlueBioPlan puede asumir la defensa exclusiva de cualquier asunto sujeto a indemnización, en cuyo caso usted acepta cooperar y no llegar a acuerdos sin nuestro consentimiento previo por escrito.',
        },
        {
          title: '11. Ley aplicable y resolución de disputas',
          content: `Estos Términos se rigen por las leyes de España, excluyendo sus normas sobre conflictos de leyes. Los tribunales de España, concretamente los situados en Murcia, tendrán jurisdicción exclusiva sobre cualquier disputa derivada de estos Términos o de su uso de Novik. Usted acepta someterse a la jurisdicción personal de dichos tribunales.

Antes de iniciar procedimientos legales formales, acepta intentar resolver cualquier disputa con nosotros de manera informal contactando con info@novik.ai. Si no se puede resolver una disputa de forma amistosa, cualquiera de las partes podrá presentar acciones ante los tribunales indicados. Sin perjuicio de lo anterior, BlueBioPlan podrá solicitar medidas cautelares u otra compensación equitativa en cualquier jurisdicción para proteger su propiedad intelectual o prevenir el uso indebido de Novik.`,
        },
        {
          title: '12. Cambios en estos Términos',
          content:
            'Podemos modificar estos Términos en cualquier momento. Publicaremos los Términos actualizados en Novik y actualizaremos la fecha de «Última actualización». Para cambios significativos, proporcionaremos un aviso adicional (por ejemplo, mediante correo electrónico o notificación en la aplicación). A menos que la ley requiera lo contrario, los cambios entrarán en vigor 15 días después del aviso. Su uso continuado de Novik después de la fecha de entrada en vigor constituye la aceptación de los Términos actualizados. Si no está de acuerdo con los cambios, debe dejar de utilizar Novik.',
        },
        {
          title: '13. Miscelánea',
          content: `Estos Términos constituyen el acuerdo completo entre usted y BlueBioPlan con respecto a Novik y sustituyen todos los acuerdos previos. Si alguna disposición se considera inválida o inaplicable, las disposiciones restantes permanecerán en pleno vigor. Nuestro incumplimiento a la hora de hacer valer cualquier derecho no implica renuncia a ese derecho. Podemos ceder estos Términos sin restricciones; usted no puede ceder sus derechos sin nuestro consentimiento. Estos Términos no crean ninguna relación de sociedad, empresa conjunta o agencia entre nosotros. Los encabezados se incluyen únicamente por comodidad.

Si tiene preguntas sobre estos Términos, póngase en contacto con info@novik.ai o escriba a BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, España.`,
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

      <EffectiveDate>{currentContent.effectiveDate}</EffectiveDate>

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

      <Box sx={{ mt: 5, pt: 3, borderTop: `1px solid ${novikTheme.colors.border}` }}>
        <Typography
          sx={{
            textAlign: 'center',
            fontSize: '0.9rem',
            fontStyle: 'italic',
            color: 'text.secondary',
          }}
        >
          {language === 'en'
            ? 'By using Novik, you acknowledge that you have read, understood and agree to be bound by these Terms of Use.'
            : 'Al utilizar Novik, usted reconoce que ha leído, comprendido y acepta estar sujeto a estos Términos de Uso.'}
        </Typography>
      </Box>
    </PageContainer>
  );
};

export default TermsOfUse;
