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

const PrivacyPolicy = () => {
  const [language, setLanguage] = useState<'en' | 'es'>('en');

  const content = {
    en: {
      title: 'Novik Privacy Policy',
      effectiveDate: 'Effective Date: 17 August 2025',
      intro:
        'Novik is an AI-driven clinical decision support tool operated by BlueBioPlan S.L. ("BlueBioPlan", "we" or "us"). We are committed to protecting your privacy. This policy describes how we collect, use, share and protect personal information when you use the Novik platform. We comply with the European Union\'s General Data Protection Regulation (GDPR), Spain\'s LOPD-GDD, the United States\' Health Insurance Portability and Accountability Act (HIPAA), the California Consumer Privacy Act/California Privacy Rights Act (CCPA/CPRA), Brazil\'s Lei Geral de Proteção de Dados (LGPD), Canada\'s Personal Information Protection and Electronic Documents Act (PIPEDA), Australia\'s Privacy Act 1988, and other applicable laws in the jurisdictions where Novik is used. BlueBioPlan is established in Spain and acts as the data controller for personal data collected through Novik.',
      sections: [
        {
          title: '1. Who We Are and Scope',
          content:
            'This Privacy Policy applies to all users of the Novik web and mobile platform worldwide. It covers personal information collected through the Novik application, our websites and related services. It does not apply to third-party websites or services that may be linked within Novik. If you do not agree with our practices, please do not use Novik.',
        },
        {
          title: '2. Data Controller and Contact Details',
          content:
            'The data controller is BlueBioPlan S.L., a company registered in Spain (CIF B-72913379) with its business address at C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, Spain. You can contact us by email at info@novik.ai or by telephone at +34 690 957 910. If required by law, we will designate a Data Protection Officer (DPO) and publish their contact details here.',
        },
        {
          title: '3. Information We Collect',
          content: `We strive to minimize the personal data we collect. Novik is designed to process only anonymized case data; we do not collect protected health information (PHI). The categories of data we collect include:

**Account Information** – When you register or sign in, we collect your email address and a hashed (encrypted) password. We may record your role (professional or student), your country or region, and your account creation date. We do not require your full name or physical address to use Novik.

**Technical and Usage Data** – We automatically collect technical data such as IP address, device identifiers, browser type, app version, operating system, user-agent string and language. We log usage information such as log-in timestamps, queries entered, features used, session duration and error reports. This helps us secure the platform and improve performance.

**User Inputs and Uploaded Files** – If you input text or upload files (e.g. case descriptions, PDFs, audio recordings or X-ray images) to Novik, our systems process that content to generate recommendations. Users must ensure that no personal identifiers are included. Under the U.S. HIPAA Privacy Rule, health information that has been de-identified is not PHI and thus is not subject to HIPAA. You are responsible for removing identifiers such as names, contact details, medical record numbers, facial images and any other elements that could identify an individual. We treat the content you submit as anonymous and do not link it to your identity.

**Cookies and Similar Technologies** – We use a minimal number of cookies to operate our website and app. Strictly necessary cookies enable core features (e.g. login, payment). With your consent, we may also use functionality or analytics cookies to remember preferences and analyze how users interact with Novik. We do not use advertising or tracking cookies. For details and to manage your choices, please see our separate Cookie Policy.

**Communications** – When you contact us (e.g. by email or through a support form), we collect the information you provide (such as your name, email address and message content). We use this to respond to you and maintain records of the communication.

We do not intentionally collect special categories of personal data (e.g. race, religion, genetic data) or payment/financial information through Novik. If payments are required in future, these will be handled by a secure third-party provider under separate terms.`,
        },
        {
          title: '4. How We Use Your Information',
          content: `We use the collected data for the following purposes, consistent with applicable law:

• **Providing and improving the service.** We process your account details to authenticate you and provide access to Novik. We analyze your inputs (in anonymized form) to generate clinical recommendations and to monitor and improve the accuracy, security and performance of Novik. We maintain logs to debug issues and prevent misuse.

• **Personalized recommendations.** Novik tailors its responses based on the anonymized case information you provide. For example, if you indicate that a patient has an allergy (without identifying the patient), the system will adjust its suggestions accordingly. These recommendations are for your professional reference only and are not connected to your identity or the patient's identity.

• **Analytics and development.** We analyze aggregated usage patterns to understand how clinicians interact with Novik. This helps us improve features and train our algorithms. We may create de-identified statistics or reports from aggregated data.

• **Communications.** We use your contact information to respond to your inquiries, notify you about changes to our policies or terms, send technical notices and, if applicable, send marketing or research communications (you may opt out of marketing at any time).

• **Compliance, safety and legal.** We may use information to comply with legal obligations, enforce our Terms of Use, protect the rights and safety of BlueBioPlan and our users, and detect and prevent fraud, security breaches or other harmful activity.`,
        },
        {
          title: '5. Legal Basis for Processing',
          content: `Under the GDPR and equivalent laws, we rely on the following legal bases:

• **Performance of a contract.** We process account information and user inputs to provide the service you request. Without this data, we cannot deliver Novik.

• **Legitimate interests.** We process technical data, usage logs and minimal contact details to secure and improve Novik. We have a legitimate interest in understanding how our service is used and ensuring its reliability and safety. These interests are balanced against your rights and freedoms.

• **Consent.** We rely on your consent to place non-essential cookies and to send you optional communications. You can withdraw consent at any time.

• **Legal compliance.** We may process data as necessary to comply with applicable laws and regulations (for example, to respond to lawful requests from authorities).

BlueBioPlan does not process identifiable health data and is not a "covered entity" or "business associate" under HIPAA because de-identified health information is not protected health information (PHI) under the HIPAA Privacy Rule. If we ever process PHI in collaboration with U.S. healthcare providers, we will do so in compliance with HIPAA (including signing Business Associate Agreements and implementing required safeguards).`,
        },
        {
          title: '6. Data Sharing and Disclosure',
          content: `We do not sell personal information and do not share it for advertising or cross-context behavioural marketing. We only share personal data as follows:

• **Service providers.** We engage trusted third parties to provide hosting, analytics, payment processing or customer support. These providers may access personal information only to perform services on our behalf and are contractually bound to protect it.

• **Legal and regulatory purposes.** We may disclose information when required to comply with a legal obligation, court order or governmental request, to enforce our Terms of Use, to protect our rights or to prevent harm to users or the public.

• **Business transfers.** In the event of a merger, acquisition or sale of assets, personal data may be transferred as part of that transaction. We will notify you of such a change and any choices you may have.

• **Aggregated or de-identified data.** We may share aggregated statistics or anonymized reports that do not identify any individual, for research or statistical purposes.`,
        },
        {
          title: '7. International Data Transfers',
          content:
            "We store personal data on servers located in the European Union. However, our service providers may operate in other countries. When we transfer data outside the European Economic Area (EEA), we use legal mechanisms such as the European Commission's Standard Contractual Clauses or rely on adequacy decisions to ensure that your data receives an adequate level of protection. By using Novik, you consent to the transfer of your information to countries outside your country of residence, which may have different data protection laws than your own.",
        },
        {
          title: '8. Data Retention',
          content:
            'We retain personal data only for as long as necessary for the purposes described in this policy, including to provide the service, comply with legal obligations, resolve disputes and enforce our agreements. When data is no longer required, we will either anonymize or securely delete it. Aggregated usage data may be retained for longer for analytics and research purposes, but it will not identify any individual.',
        },
        {
          title: '9. Security Measures',
          content:
            'We implement appropriate technical and organisational measures to protect personal information from unauthorized access, disclosure, alteration or destruction. These measures include encryption, access controls, secure development practices and regular monitoring. However, no method of transmission or storage is completely secure, and we cannot guarantee absolute security. You are responsible for safeguarding your account credentials and for promptly notifying us of any suspected security incident.',
        },
        {
          title: '10. Your Rights',
          content: `Depending on your location, you have specific rights regarding your personal data.

**Users in the European Union and United Kingdom**
Under the GDPR/UK GDPR you have the right to: (i) access the personal data we hold about you; (ii) rectify inaccurate or incomplete data; (iii) erase your data in certain circumstances ("right to be forgotten"); (iv) restrict or object to processing; (v) data portability, i.e. obtain your data in a structured, commonly used format; and (vi) withdraw consent at any time for processing based on consent. You also have the right to lodge a complaint with your local data protection authority.

**Users in California (CCPA/CPRA)**
California residents have the right to: (i) request disclosure of the categories and specific pieces of personal information we collect, use and disclose; (ii) request deletion of personal information (subject to exceptions); (iii) request correction of inaccurate personal information; (iv) opt out of the sale or sharing of personal information; and (v) not be discriminated against for exercising privacy rights. We do not sell personal information and do not share it for cross-context behavioural advertising. To exercise your California rights, contact us using the details below.

**Users in Brazil (LGPD)**
Under the LGPD you may: (i) confirm whether we process your personal data; (ii) request access to your data; (iii) correct incomplete, inaccurate or outdated data; (iv) anonymize, block or delete unnecessary or excessive data; (v) request portability to another provider; (vi) request information about public and private entities with which we share data; (vii) revoke consent; and (viii) lodge a complaint with the Brazilian data protection authority (ANPD).

**Users in Canada and Australia**
Canadian users have the right to access the personal information we hold and to challenge its accuracy under PIPEDA. Australian users may request access to and correction of their personal information and may complain to the Office of the Australian Information Commissioner (OAIC) if unsatisfied.

**Exercising Your Rights**
You may exercise your rights by contacting us at info@novik.ai. We may need to verify your identity before processing your request. We will respond within the timeframe required by applicable law (e.g. 30 days under GDPR, 45 days under CCPA/CPRA). We may refuse requests where permitted by law (for example, if fulfilling the request would infringe the rights of others or reveal confidential information). There is no charge for making a request unless it is unfounded or excessive, in which case we may charge a reasonable fee or refuse to comply.`,
        },
        {
          title: "11. Children's Privacy",
          content:
            'Novik is intended for use only by adults. We do not knowingly collect personal information from children under 18 years of age. If we become aware that a minor has provided us with personal data, we will delete it. Parents or guardians who believe that their child has provided personal data to us may contact us to request removal.',
        },
        {
          title: '12. Changes to This Policy',
          content:
            'We may update this Privacy Policy from time to time to reflect changes in our practices, technologies or legal requirements. When we make changes, we will revise the Effective Date at the top of the policy. For material changes, we may provide a more prominent notice or seek your consent where required. We encourage you to review this policy periodically. Your continued use of Novik after an update signifies acceptance of the revised policy.',
        },
        {
          title: '13. Contact Us',
          content:
            'If you have any questions, concerns or requests regarding this Privacy Policy or your personal data, please contact us at info@novik.ai or by postal mail at BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, Spain. If you are not satisfied with our response, you may lodge a complaint with your local data protection authority.',
        },
      ],
    },
    es: {
      title: 'Política de Privacidad de Novik',
      effectiveDate: 'Fecha de entrada en vigor: 17 de agosto de 2025',
      intro:
        'Novik es una herramienta de apoyo a la decisión clínica impulsada por IA operada por BlueBioPlan S.L. ("BlueBioPlan", "nosotros" o "nos"). Estamos comprometidos a proteger su privacidad. Esta política describe cómo recopilamos, usamos, compartimos y protegemos la información personal cuando utiliza la plataforma Novik. Cumplimos con el Reglamento General de Protección de Datos de la Unión Europea (RGPD), la LOPD-GDD de España, la Ley de Portabilidad y Responsabilidad del Seguro Médico de Estados Unidos (HIPAA), la Ley de Privacidad del Consumidor de California/Ley de Derechos de Privacidad de California (CCPA/CPRA), la Lei Geral de Proteção de Dados de Brasil (LGPD), la Ley de Protección de Información Personal y Documentos Electrónicos de Canadá (PIPEDA), la Ley de Privacidad de 1988 de Australia, y otras leyes aplicables en las jurisdicciones donde se utiliza Novik. BlueBioPlan está establecida en España y actúa como responsable del tratamiento de los datos personales recopilados a través de Novik.',
      sections: [
        {
          title: '1. Quiénes somos y alcance',
          content:
            'Esta Política de Privacidad se aplica a todos los usuarios de la plataforma web y móvil de Novik en todo el mundo. Cubre la información personal recopilada a través de la aplicación Novik, nuestros sitios web y los servicios relacionados. No se aplica a sitios o servicios de terceros que puedan estar vinculados dentro de Novik. Si no está de acuerdo con nuestras prácticas, no utilice Novik.',
        },
        {
          title: '2. Responsable del tratamiento y datos de contacto',
          content:
            'El responsable del tratamiento es BlueBioPlan S.L., una empresa registrada en España (CIF B-72913379) con domicilio social en C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, España. Puede contactarnos por correo electrónico en info@novik.ai o por teléfono en +34 690 957 910. Si la ley lo exige, designaremos un Delegado de Protección de Datos (DPD) y publicaremos aquí sus datos de contacto.',
        },
        {
          title: '3. Información que recopilamos',
          content: `Nos esforzamos por minimizar los datos personales que recopilamos. Novik está diseñado para procesar únicamente datos de casos anonimizados; no recopilamos información sanitaria protegida (PHI). Las categorías de datos que recopilamos incluyen:

**Información de la cuenta** – Cuando se registra o inicia sesión, recopilamos su dirección de correo electrónico y una contraseña cifrada. Podemos registrar su rol (profesional o estudiante), su país o región y la fecha de creación de su cuenta. No requerimos su nombre completo ni su dirección física para utilizar Novik.

**Datos técnicos y de uso** – Recopilamos automáticamente datos técnicos como la dirección IP, identificadores de dispositivo, tipo de navegador, versión de la aplicación, sistema operativo, cadena de agente de usuario e idioma. Registramos información de uso como marcas de tiempo de inicio de sesión, consultas introducidas, funciones utilizadas, duración de la sesión e informes de errores. Esto nos ayuda a proteger la plataforma y mejorar su rendimiento.

**Entradas del usuario y archivos cargados** – Si introduce texto o carga archivos (por ejemplo, descripciones de casos, PDFs, grabaciones de audio o imágenes de radiografías) en Novik, nuestros sistemas procesan ese contenido para generar recomendaciones. Los usuarios deben asegurarse de que no se incluyan identificadores personales. Según la normativa de EE. UU., la información sanitaria que se ha desidentificado no se considera PHI y, por tanto, no está sujeta a la HIPAA. Usted es responsable de eliminar identificadores como nombres, datos de contacto, números de historia clínica, imágenes de rostros y cualquier otro elemento que pueda identificar a una persona. Tratamos el contenido que envía como anónimo y no lo vinculamos a su identidad.

**Cookies y tecnologías similares** – Utilizamos un número mínimo de cookies para hacer funcionar nuestro sitio web y la aplicación. Las cookies estrictamente necesarias habilitan funciones básicas (por ejemplo, inicio de sesión, pago). Con su consentimiento, también podemos utilizar cookies de funcionalidad o analítica para recordar preferencias y analizar cómo interactúan los usuarios con Novik. No utilizamos cookies publicitarias ni de seguimiento. Para obtener más detalles y administrar sus opciones, consulte nuestra Política de Cookies.

**Comunicaciones** – Cuando se comunica con nosotros (por ejemplo, por correo electrónico o mediante un formulario de soporte), recopilamos la información que proporciona (como su nombre, dirección de correo electrónico y el contenido del mensaje). Utilizamos esto para responderle y mantener un registro de la comunicación.

No recopilamos intencionadamente categorías especiales de datos personales (p. ej., raza, religión, datos genéticos) ni información financiera/pagos a través de Novik. Si en el futuro se requieren pagos, estos serán gestionados por un proveedor tercero seguro en virtud de términos separados.`,
        },
        {
          title: '4. Cómo utilizamos su información',
          content: `Utilizamos los datos recopilados para los siguientes fines, de conformidad con la ley aplicable:

• **Proporcionar y mejorar el servicio.** Tratamos sus datos de cuenta para autenticarle y darle acceso a Novik. Analizamos sus entradas (de forma anonimizada) para generar recomendaciones clínicas y para supervisar y mejorar la precisión, seguridad y rendimiento de Novik. Mantenemos registros para depurar problemas y prevenir usos indebidos.

• **Recomendaciones personalizadas.** Novik adapta sus respuestas en función de la información clínica anonimizada que usted proporciona. Por ejemplo, si indica que un paciente tiene una alergia (sin identificar al paciente), el sistema ajustará sus sugerencias en consecuencia. Estas recomendaciones son solo para su referencia profesional y no están conectadas a su identidad ni a la del paciente.

• **Analítica y desarrollo.** Analizamos patrones de uso agregados para comprender cómo interactúan los clínicos con Novik. Esto nos ayuda a mejorar las funciones y entrenar nuestros algoritmos. Podemos crear estadísticas o informes desidentificados a partir de datos agregados.

• **Comunicaciones.** Utilizamos su información de contacto para responder a sus consultas, notificarle sobre cambios en nuestras políticas o términos, enviar avisos técnicos y, si corresponde, enviar comunicaciones de marketing o investigación (puede optar por no recibir marketing en cualquier momento).

• **Cumplimiento, seguridad y legalidad.** Podemos utilizar la información para cumplir con obligaciones legales, hacer cumplir nuestros Términos de Uso, proteger los derechos y la seguridad de BlueBioPlan y nuestros usuarios, y detectar y prevenir fraude, fallos de seguridad u otras actividades dañinas.`,
        },
        {
          title: '5. Base legal del tratamiento',
          content: `De acuerdo con el GDPR y leyes equivalentes, nos basamos en las siguientes bases legales:

• **Ejecución de un contrato.** Tratamos información de la cuenta y las entradas del usuario para proporcionar el servicio que solicita. Sin estos datos, no podemos prestar Novik.

• **Intereses legítimos.** Tratamos datos técnicos, registros de uso y datos de contacto mínimos para proteger y mejorar Novik. Tenemos un interés legítimo en comprender cómo se utiliza nuestro servicio y garantizar su fiabilidad y seguridad. Estos intereses se equilibran con sus derechos y libertades.

• **Consentimiento.** Dependemos de su consentimiento para establecer cookies no esenciales y para enviarle comunicaciones opcionales. Puede retirar su consentimiento en cualquier momento.

• **Cumplimiento legal.** Podemos tratar datos según sea necesario para cumplir con leyes y regulaciones aplicables (por ejemplo, para responder a solicitudes de autoridades).

BlueBioPlan no procesa datos de salud identificables y no es una «entidad cubierta» ni un «asociado comercial» según HIPAA, porque la información de salud desidentificada no es información sanitaria protegida (PHI) bajo la normativa HIPAA. Si en el futuro procesamos PHI en colaboración con proveedores de EE. UU., lo haremos en cumplimiento de HIPAA (incluida la firma de acuerdos de asociado comercial y la implementación de medidas de seguridad requeridas).`,
        },
        {
          title: '6. Compartición y divulgación de datos',
          content: `No vendemos información personal ni la compartimos con fines de publicidad o marketing conductual. Solo compartimos datos personales en los siguientes casos:

• **Proveedores de servicios.** Contratamos a terceros de confianza para proporcionar alojamiento, analítica, procesamiento de pagos o soporte al cliente. Estos proveedores pueden acceder a la información personal solo para prestar servicios en nuestro nombre y están contractualmente obligados a protegerla.

• **Fines legales y regulatorios.** Podemos divulgar información cuando sea necesario para cumplir con una obligación legal, una orden judicial o una solicitud gubernamental, para hacer cumplir nuestros Términos de Uso, para proteger nuestros derechos o para prevenir daños a usuarios o al público.

• **Transferencias empresariales.** En caso de fusión, adquisición o venta de activos, los datos personales pueden transferirse como parte de la transacción. Le notificaremos dicho cambio y cualquier opción que pueda tener.

• **Datos agregados o desidentificados.** Podemos compartir estadísticas agregadas o informes anonimizados que no identifiquen a ninguna persona, para fines de investigación o estadísticos.`,
        },
        {
          title: '7. Transferencias internacionales de datos',
          content:
            'Almacenamos los datos personales en servidores ubicados en la Unión Europea. Sin embargo, nuestros proveedores pueden operar en otros países. Cuando transferimos datos fuera del Espacio Económico Europeo (EEE), utilizamos mecanismos legales como las Cláusulas Contractuales Tipo de la Comisión Europea o dependemos de decisiones de adecuación para garantizar que sus datos reciben un nivel de protección adecuado. Al utilizar Novik, usted acepta la transferencia de su información a países fuera de su país de residencia, que pueden tener leyes de protección de datos diferentes a las suyas.',
        },
        {
          title: '8. Conservación de datos',
          content:
            'Conservamos los datos personales solo durante el tiempo necesario para los fines descritos en esta política, incluidos para proporcionar el servicio, cumplir con obligaciones legales, resolver disputas y hacer cumplir nuestros acuerdos. Cuando los datos ya no sean necesarios, los anonimizaremos o eliminaremos de forma segura. Los datos de uso agregados pueden conservarse por más tiempo para análisis e investigación, pero no identificarán a ninguna persona.',
        },
        {
          title: '9. Medidas de seguridad',
          content:
            'Implementamos medidas técnicas y organizativas apropiadas para proteger la información personal contra el acceso no autorizado, divulgación, alteración o destrucción. Estas medidas incluyen cifrado, controles de acceso, prácticas de desarrollo seguro y supervisión regular. Sin embargo, ningún método de transmisión o almacenamiento es completamente seguro y no podemos garantizar una seguridad absoluta. Usted es responsable de proteger las credenciales de su cuenta y de notificarnos de inmediato cualquier incidente de seguridad sospechoso.',
        },
        {
          title: '10. Sus derechos',
          content: `Dependiendo de su ubicación, tiene derechos específicos respecto a sus datos personales.

**Usuarios en la Unión Europea y el Reino Unido**
Bajo el GDPR/UK GDPR tiene derecho a: (i) acceder a los datos personales que tenemos sobre usted; (ii) rectificar datos inexactos o incompletos; (iii) eliminar sus datos en determinadas circunstancias («derecho al olvido»); (iv) restringir u oponerse al tratamiento; (v) portabilidad de datos, es decir, obtener sus datos en un formato estructurado de uso común; y (vi) retirar el consentimiento en cualquier momento para el tratamiento basado en el consentimiento. También tiene derecho a presentar una reclamación ante su autoridad de protección de datos local.

**Usuarios en California (CCPA/CPRA)**
Los residentes de California tienen derecho a: (i) solicitar la divulgación de las categorías y piezas específicas de información personal que recopilamos, usamos y divulgamos; (ii) solicitar la eliminación de la información personal (sujeto a excepciones); (iii) solicitar la corrección de información personal inexacta; (iv) optar por no vender ni compartir información personal; y (v) no ser discriminados por ejercer sus derechos de privacidad. No vendemos información personal ni la compartimos para publicidad conductual. Para ejercer sus derechos en California, póngase en contacto con nosotros.

**Usuarios en Brasil (LGPD)**
En virtud de la LGPD, puede: (i) confirmar si tratamos sus datos personales; (ii) solicitar acceso a sus datos; (iii) corregir datos incompletos, inexactos o desactualizados; (iv) anonimizar, bloquear o eliminar datos innecesarios o excesivos; (v) solicitar la portabilidad a otro proveedor; (vi) solicitar información sobre las entidades públicas y privadas con las que compartimos datos; (vii) revocar el consentimiento; y (viii) presentar una queja ante la autoridad brasileña de protección de datos (ANPD).

**Usuarios en Canadá y Australia**
Los usuarios canadienses tienen derecho a acceder a la información personal que tenemos y a impugnar su exactitud en virtud de PIPEDA. Los usuarios australianos pueden solicitar el acceso y la corrección de su información personal y pueden presentar una queja ante la Oficina del Comisionado de Información de Australia (OAIC) si no están satisfechos.

**Ejercicio de sus derechos**
Puede ejercer sus derechos poniéndose en contacto con nosotros en info@novik.ai. Es posible que necesitemos verificar su identidad antes de procesar su solicitud. Responderemos dentro del plazo requerido por la ley aplicable (por ejemplo, 30 días según el GDPR, 45 días según CCPA/CPRA). Podemos rechazar solicitudes en los casos permitidos por la ley (por ejemplo, si cumplir la solicitud infringiría los derechos de otras personas o revelaría información confidencial). No se cobrará ninguna tarifa por realizar una solicitud a menos que sea infundada o excesiva; en ese caso, podemos cobrar una tarifa razonable o negarnos a cumplir.`,
        },
        {
          title: '11. Privacidad de los menores',
          content:
            'Novik está destinado solo a adultos. No recopilamos intencionadamente información personal de menores de 18 años. Si nos damos cuenta de que un menor nos ha proporcionado datos personales, los eliminaremos. Los padres o tutores que crean que su hijo nos ha proporcionado datos personales pueden ponerse en contacto con nosotros para solicitar su eliminación.',
        },
        {
          title: '12. Cambios en esta política',
          content:
            'Podemos actualizar esta Política de Privacidad periódicamente para reflejar cambios en nuestras prácticas, tecnologías o requisitos legales. Cuando realicemos cambios, revisaremos la Fecha de entrada en vigor al principio de la política. Para cambios significativos, podremos proporcionar un aviso más destacado o solicitar su consentimiento si así se requiere. Le animamos a revisar esta política periódicamente. Su uso continuado de Novik tras una actualización implica la aceptación de la política revisada.',
        },
        {
          title: '13. Contacto',
          content:
            'Si tiene alguna pregunta, inquietud o solicitud respecto a esta Política de Privacidad o sus datos personales, póngase en contacto con nosotros en info@novik.ai o por correo postal en BlueBioPlan S.L., C/ La Cruz 2, Entresuelo, 30820 Alcantarilla, Murcia, España. Si no está satisfecho con nuestra respuesta, puede presentar una queja ante su autoridad de protección de datos local.',
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
            if (paragraph.startsWith('•')) {
              return (
                <ul key={pIndex}>
                  {paragraph.split('\n').map((item, iIndex) => {
                    if (item.trim().startsWith('•')) {
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
                    }
                    return null;
                  })}
                </ul>
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

export default PrivacyPolicy;
