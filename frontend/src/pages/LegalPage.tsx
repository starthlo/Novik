import { useEffect, useState } from 'react';
import { useLocation } from 'react-router-dom';
import Header from '../components/Common/Header';

const sections = [
  { id: 'terms', label: 'Terms of Service' },
  { id: 'privacy', label: 'Privacy Policy' },
  { id: 'cookies', label: 'Cookie Policy' },
];

const LegalPage = () => {
  const [active, setActive] = useState('terms');
  const location = useLocation();

  useEffect(() => {
    const hash = location.hash.replace('#', '');
    if (sections.some(s => s.id === hash)) {
      setActive(hash);
    }
  }, [location]);

  const content = {
    terms: (
      <div className="space-y-4 text-gray-700 text-md leading-relaxed">
        <h2 className="text-xl font-semibold">I. GENERAL INFORMATION</h2>
        <p>
          In compliance with the duty of information set forth by Spanish Law 34/2002 of July 11 on
          Services of the Information Society and Electronic Commerce (LSSI-CE), the following
          general information about this website is provided:
        </p>
        <p>
          The ownership of the website www.novik.ai is held by BlueBioPlan S.L.U., with Tax ID (NIF)
          B-72913379, registered in the Commercial Registry of Murcia (Spain) in Volume 3615, Book
          0, Page 218, Sheet MU-109359, Entry 1. The company is represented by Mr. Vicente Ferrer
          Pérez. The contact details of the owner are:
        </p>
        <ul>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Registered Address:</strong> C/ La Cruz 2, Entresuelo,
            30820 Alcantarilla, Murcia, Spain.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Contact Phone:</strong> +34 690 957 910.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Contact Email:</strong>{' '}
            <a href="mailto:info@novik.ai">info@novik.ai</a>
          </li>
        </ul>
        <p>
          <strong>Privacy Policy</strong>Details regarding the processing of personal data of users
          are governed by our Website’s Privacy Policy, which adheres to the EU General Data
          Protection Regulation (GDPR) and other applicable data protection laws.
        </p>
        <h2 className="text-xl font-semibold mt-6">II. GENERAL TERMS AND CONDITIONS OF USE</h2>
        <p>
          Purpose of the Conditions – The Website: These General Terms and Conditions of Use
          (hereinafter, Terms) regulate the access to and use of the Website. For the purpose of
          these Terms, Website is understood as the external display of the content interface (both
          static and dynamic, i.e., the navigation tree), and all content and elements integrated
          therein (hereinafter, Content), as well as all online services or resources that might be
          offered to users (hereinafter, Services).
        </p>
        <p>
          Description of Service: Novik is a global AI-powered tool intended for medical and dental
          professionals. It provides informational responses and recommendations that do not replace
          clinical judgment or a personalized professional evaluation. The use of Novik is meant to
          support, not substitute, the decision-making of qualified healthcare providers, and under
          no circumstances should it be seen as a final diagnosis or treatment plan.
        </p>
        <p>
          Novik reserves the right to modify, at any time and without prior notice, the presentation
          and configuration of the Website, as well as the Content and Services included on it. The
          User acknowledges and accepts that Novik may, at any time, interrupt, deactivate, or
          cancel any of these elements or the access to them.
        </p>
        <p>
          Access to the Website by the User is free and generally does not require prior
          registration. Simply browsing the public sections of the site does not require payment or
          user registration (aside from the cost of the User’s Internet connection). However, the
          use of certain Content or Services (for example, full access to Novik’s AI consultation
          features or professional tools) may be subject to prior subscription, registration, or
          payment, which will be indicated in those specific cases.
        </p>
        <p>
          The User: Accessing, browsing, or using the Website confers the status of User. From the
          moment of accessing the Website, the User accepts all the Terms set forth herein
          (including subsequent updates), without prejudice to the application of any mandatory
          legal regulations that may apply. Given the importance of these Terms, the User is advised
          to read them carefully each time they visit the Website.
        </p>
        <p>
          Novik’s Website provides a wide range of information, services, and data. The User assumes
          responsibility for proper useof the Website. This responsibility extends to, but is not
          limited to:
        </p>
        <ul>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>RLawful and Appropriate Use:</strong> Using the
            information, Content, and/or Services offered by Novik in accordance with these Terms,
            the law, good morals, and public order. The User agrees not to use the Website, Content,
            or Services for purposes or effects that are illicit, harmful to the rights and
            interests of third parties, or that may damage, disable, overburden, or impair the
            normal use of the Services, servers, or any networks connected to the Services.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Accuracy of User Information:</strong> Ensuring that all
            information provided by the User (for instance, when filling out any registration or
            contact form) is truthful, accurate, and lawful. The User is responsible for keeping
            such information updated. If the User provides or disseminates information that is
            false, incomplete, or not updated, Novik assumes no responsibility for any consequences
            that may ensue. The User also undertakes to immediately notify Novik of any unauthorized
            use of their account or credentials, or any security breach related to the Website, of
            which they become aware, in order to enable Novik to take appropriate measures.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Age and Capacity:</strong>
            The User declares to be of legal age (18 or older) and to have the legal capacity to be
            bound by these Terms. This Website is not intended for minors. Novik disclaims any
            liability for use of the Website by persons who do not meet this age requirement.
            Moreover, Novik is intended for professional use by qualified healthcare practitioners;
            use of the platform by individuals who are not healthcare professionals is allowed only
            under their own responsibility and does not alter the disclaimers and limitations
            outlined in these Terms.
          </li>
        </ul>
        <p>
          Mere access to the Website does not imply the establishment of any commercial or
          professional relationship between Novik (BlueBioPlan S.L.U.) and the User.
        </p>
        <h2 className="text-xl font-semibold mt-6">
          III. ACCESS AND BROWSING ON THE WEBSITE: DISCLAIMER OF WARRANTIES AND LIABILITY
        </h2>
        <p>
          Novik does not warrant the continuity, availability, or usefulness of the Website, its
          Content, or its Services. Novik will make reasonable efforts to ensure that the Website is
          accessible and functions properly, but it does not guarantee that access will be
          uninterrupted or free from errors.
        </p>
        <p>
          Similarly, Novik does not guarantee that the content or software accessible through the
          Website is error-free or free of harmful components. The User understands that any content
          downloaded or otherwise obtained through the use of the Website is obtained at their own
          discretion and risk. Under no circumstances will Novik be responsible for any damage,
          loss, or harm of any kind arising from accessing, browsing, and using the Website,
          including but not limited to those caused to computer systems or those caused by viruses
          or other malicious code.
        </p>
        <p>
          Novik shall not be liable for any damages that may arise to Users due to improper use of
          the Website. In particular, it is not responsible for any damage or harm caused by
          telecommunications failures, interruptions, computer viruses, malware or system outages
          that are beyond Novik’s control.
        </p>
        <p>
          The Website and its Content are provided on an “as is” and “as available” basis, without
          warranties of any kind, either express or implied, including, but not limited to, implied
          warranties of merchantability, fitness for a particular purpose, or non-infringement. The
          User assumes all responsibility and risk for the use of the Website.
        </p>
        <h2 className="text-xl font-semibold mt-6">IV. LINKING POLICY</h2>
        <p>
          The Website may include links (such as hyperlinks, banners, buttons, directories, or
          search tools) to third-party websites or resources which are not operated by Novik. These
          links are provided solely as a convenience to the User, to facilitate access to
          information, content, and services that may be of interest available on the Internet. The
          inclusion of these external links does not imply any kind of association, merger,
          endorsement, or recommendation by Novik of the linked websites or their content.
        </p>
        <p>
          Novik has no control over and does not manage any of these external sites or resources.
          Therefore, Novik does not assume any responsibility for the availability, contents,
          advertising, products, services or any other material on or available from such external
          sites or resources. The User acknowledges and agrees that Novik shall not be responsible
          or liable, directly or indirectly, for any damage or loss caused by or in connection with
          the use of or reliance on any such content, goods, or services available on or through any
          such linked site or resource.
        </p>
        <p>
          Novik also does not guarantee the technical availability, accuracy, truthfulness, validity
          or legality of external sites accessible through links. Novik will not, and is under no
          obligation to, monitor the content of other websites.
        </p>
        <p>
          Outgoing links (from Novik to other sites): The Website may contain links to third-party
          websites. If the User follows a link to any of these websites, they should note that these
          websites have their own terms and privacy policies, and Novik does not accept any
          responsibility or liability for those policies or for any content on those websites.
        </p>
        <p>
          Incoming links (from other sites to Novik): Third parties who intend to include a
          hyperlink to www.novik.ai on their own websites or digital platforms must comply with the
          following conditions:
        </p>
        <ul>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;No prior authorization is required for linking to the homepage
            www.novik.ai. However, deep linking or framing (displaying Novik content within another
            site) is prohibited without express written consent from Novik.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;The link must be full and complete (directed to the full URL of
            Novik’s homepage). Under no circumstance may the linking website reproduce, in any form,
            the Content of the Website or parts thereof, nor may it establish a browser or border
            environment on the Content of the Website.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;False, inaccurate, or incorrect statements about Novik, the
            Website, or its services and/or content are not permitted. The linking site shall not
            state or imply that Novik has authorized the link or has any relationship or affiliation
            with that site without Novik’s express authorization.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;The linking website shall not contain any content that could be
            considered as distasteful, offensive or controversial, and should contain only content
            that is appropriate for all age groups. Specifically, it must not host content that is
            illicit, or contrary to morals, good customs, public order, or any third-party rights.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;The existence of a hyperlink on a third-party site does not
            imply in any case the existence of a relationship between Novik and the owner of that
            website, nor the acceptance or approval by Novik of its contents or services. Novik
            reserves the right to request the removal of links to its Website that it deems not to
            comply with the conditions set forth in this section.
          </li>
        </ul>
        <h2 className="text-xl font-semibold mt-6">V. INTELLECTUAL AND INDUSTRIAL PROPERTY</h2>
        <p>
          Novik, either as the owner or as an assignee, holds all intellectual and industrial
          property rights of the Website, as well as of all the Content and elements included
          therein (by way of example, but not limited to: images, audio, video, software, text,
          trademarks, logos, color combinations, structure and design, the selection of materials
          used, computer programs necessary for its operation, underlying source code, etc.). These
          works are protected by intellectual property laws of Spain, EU regulations, and
          international treaties. All rights are reserved.
        </p>
        <p>
          Pursuant to the provisions of Intellectual Property Law, the reproduction, distribution
          and public communication, including the mode of making available, of all or part of the
          Content of this Website, for commercial purposes, in any medium and by any technical
          means, without the authorization of Novik, are expressly prohibited. The User agrees to
          respect Novik’s intellectual and industrial property rights.
        </p>
        <p>
          The User is only permitted to view and obtain a temporary private copy of the Content for
          their exclusive personal and private use in their computer systems (software and
          hardware), and not subsequently communicated to third parties. Provided that it is not for
          commercial or professional purposes, and solely to allow the use of the Service, the User
          may print, download and store portions of the Content of the Website for personal use.
          Under no circumstances shall this mean an authorization or license to the User to exploit
          or use the Content for any public or commercial purpose, or any right of alteration,
          modification or decompilation of the Content beyond what is allowed by mandatory law.
        </p>
        <p>
          The User must refrain from deleting, altering, evading or manipulating any protection
          device or security system installed on the Website.
        </p>
        <p>
          In the event that any User or third party considers that any of the Content on the Website
          infringes their intellectual or industrial property rights, they should immediately notify
          Novik via the contact details provided in the General Information section of this Legal
          Notice and Terms of Use, providing sufficient detail to enable Novik to identify and, if
          applicable, verify the alleged infringement and take appropriate measures.
        </p>
        <h2 className="text-xl font-semibold mt-6">
          VI. LEGA L ACTIONS, APP LICAB LE LAW AND JURISDICTION
        </h2>
        <p>
          Novik reserves the right to bring any civil or criminal actions that it deems necessary
          for the improper use of the Website or for breach of these Terms.
        </p>
        <p>
          The relationship between the User and Novik shall be governed by the laws in force in
          Spain. Should any dispute arise in relation to the interpretation or application of these
          Terms, the parties expressly submit themselves (with waiver of any other forum that may
          correspond to them) to the jurisdiction of the Courts and Tribunals of the city of Murcia
          (Spain), except in those cases where the user is legally considered a consumer, in which
          case any dispute shall be resolved in the courts of the consumer’s place of residence or
          domicile.
        </p>
        <p>
          International Use: Novik’s platform is accessible globally. While these Terms choose
          Spanish law and jurisdiction, Novik endeavors to align with relevant legal requirements
          across different jurisdictions to the extent applicable (especially concerning user data
          protection and privacy rights). However, nothing in this clause shall be interpreted as an
          explicit or implicit submission of Novik to any jurisdiction other than the one stated
          above, except where mandatory applicable laws provide otherwise.
        </p>
        <h2 className="text-xl font-semibold mt-6">
          VI. LEGA L ACTIONS, APP LICAB LE LAW AND JURISDICTION
        </h2>
        <p>
          Novik is an AI-assisted tool designed to provide information and recommendations in the
          dental, stomatological, and medical fields. It is crucial to understand that Novik does
          not perform clinical diagnoses nor make final decisionsregarding medical or dental
          treatments. The information provided by Novik is for guidance purposes only and should
          never be used as a substitute for consultation, diagnosis, or treatment by a qualified
          healthcare professional.
        </p>
        <p>
          In view of the above, the following clarifications apply to the use of Novik as a clinical
          decision support tool:
        </p>
        <ul>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>User’s Responsibility:</strong> The User of Novik is
            solely responsible for interpreting the information provided by the platform and for
            making decisions based on their own professional judgment or after consulting with a
            qualified healthcare professional. Novik assumes no liability for any actions the User
            takes or fails to take based on the answers or information obtained through the
            platform. Any clinical or professional decisions made by the User are made at the User’s
            own risk and responsibility.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Professional Consultation Recommended:</strong>It is
            strongly recommended that any clinical action or decision on patient treatment be
            confirmed and validated by a licensed dentist, physician, or appropriate healthcare
            professional. Novik’s responses are intended to assist healthcare professionals but do
            not replace medical/dental judgment and should not be considered a definitive diagnosis
            or sole basis for a treatment plan. Users should use the outputs as supplementary
            information and always consider the individual patient’s situation.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>No Medical Advice or Guarantee:</strong> The information
            provided by Novik does not constitute medical or dental advice. Novik disclaims any
            liability for consequences arising from the application of recommendations or
            suggestions provided by the AI. The User acknowledges that any reliance on such
            information is at their own discretion. Any clinical outcomes (positive or negative)
            resulting from following advice obtained via Novik are beyond Novik’s control and remain
            the responsibility of the User.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Use Under Professional Supervision:</strong> Novik is
            intended to be used exclusively by healthcare professionals or under their direct
            supervision. Use of the platform by any other person, or in any context outside a
            professional healthcare environment, is undertaken at the User’s sole risk. Novik will
            not be liable for any harm or damage resulting from improper use of the platform,
            including use by individuals who lack the necessary medical or dental training or who
            apply the information without proper professional oversight.
          </li>
        </ul>
        <h2 className="text-xl font-semibold mt-6">
          VIII. ADDITIONAL DISCLAIMER OF WARRANTIES AND LIABILITY
        </h2>
        <p>
          While the Novik team strives to provide accurate and up-to-date information, the platform
          and all its Services are provided on an “as is” basis. Novik makes no express or implied
          warranties or representations regarding the information, content, and Services, including
          but not limited to warranties of accuracy, reliability, completeness, or timeliness of the
          information provided. Novik does not guarantee that the information will be appropriate or
          effective for every particular scenario or that it will be error-free.
        </p>
        <p>
          The User understands and agrees that, to the fullest extent permitted by applicable law,
          in no event shall Novik (BlueBioPlan S.L.U.) be liable for any direct, indirect,
          incidental, special, consequential, or exemplary damages – including but not limited to,
          damages for loss of profits, goodwill, data, or other intangible losses – resulting from
          the use of or inability to use the platform, even if Novik has been advised of the
          possibility of such damages.
        </p>
        <p>
          Specifically, Novik will not be liable for damages arising from: (i) any errors or
          omissions in the content provided by the AI, (ii) any decision made or action taken by the
          User in reliance upon any information or advice provided by Novik, (iii) downtime or
          unavailability of the Service, (iv) unauthorized access to or alteration of User
          transmissions or data, or (v) any other matter relating to use of the Website or Services.
        </p>
        <p>
          Nothing in these Terms shall limit or exclude any liability that cannot be limited or
          excluded by law. If the User is dissatisfied with the Service or the Terms, their sole and
          exclusive remedy is to discontinue using Novik.
        </p>
        <h2 className="text-xl font-semibold mt-6">IX. USER RESPONSIBILITY FOR UPLOADED CONTENT</h2>
        <p>
          If the User uploads or submits any content or data to the Novik platform (for example,
          patient case information, images such as X-rays or photographs, medical/dental history
          details, or other data for analysis), the User is solely responsible for ensuring that all
          such content is thoroughly anonymized or de-identified prior to upload. This obligation
          includes, but is not limited to, removing any personally identifiable information such as
          names, addresses, contact information, identification numbers, facial images, or any other
          data that could directly or indirectly identify an individual (whether the User or a third
          party, such as a patient).
        </p>
        <p>
          Novik does not proactively monitor content uploaded by Users and operates under the
          assumption that Users comply with applicable data protection laws and ethical guidelines
          by anonymizing any personal data. Any failure to anonymize data before uploading is the
          sole responsibility of the User. Novik shall not be liable for any consequences, including
          legal penalties or claims for privacy breaches, that may result from the User’s failure to
          properly anonymize or obtain necessary consents for the data they provide to the platform.
        </p>
        <p>
          The User agrees to indemnify and hold harmless Novik (BlueBioPlan S.L.U.) and its
          affiliates, officers, employees, and agents from and against any claim, demand,
          investigation, or legal action arising from personal data or sensitive information that
          the User uploads to the platform without appropriate anonymization or authorization. The
          User will bear any costs (including reasonable attorneys’ fees) and damages incurred by
          Novik as a result of such a breach of the User’s obligations in this section.
        </p>
        <p>
          In simpler terms, the User should only upload information that is stripped of personal
          identifiers. If any personal data is inadvertently included, the User must have obtained
          all necessary consents from the data subject, and even then, should exercise extreme
          caution considering the sensitive nature of health information. Novik’s liability in
          handling user-provided data is limited to acting as a data processor under the User’s
          direction, and the User, as the data provider, assumes the role of data controller for any
          personal data they input, with all associated legal responsibilities under applicable data
          protection laws (GDPR, HIPAA, etc., as relevant).
        </p>
        <h2 className="text-xl font-semibold mt-6">X. USE OF UPLOADED DATA</h2>
        <p>
          By uploading content to Novik, the User grants Novik (BlueBioPlan S.L.U.) a non-exclusive,
          irrevocable, worldwide, royalty-free, and sublicensable license to use that content
          strictly for purposes related to the provision, maintenance, and improvement of Novik’s
          services. Such purposes include, for example: performing aggregated statistical analysis,
          refining algorithms and AI models, improving answer accuracy, developing new features, and
          conducting internal research to enhance service quality.
        </p>
        <p>
          This license is granted for the maximum duration permitted by applicable law. Importantly,
          this license to use User-uploaded content is limited to internal and service- related
          objectives. Under no circumstances will Novik use or exploit anonymized user-uploaded data
          for commercial purposes (such as marketing or selling data to third parties) without
          obtaining the User’s explicit informed consent.
        </p>
        <p>
          Novik commits to strict compliance with data protection regulations when handling any data
          provided by Users. Specifically, Novik will not use any uploaded data that qualifies as
          personal data (i.e., information relating to an identified or identifiable natural person)
          except in accordance with the applicable data protection laws (such as the EU’s GDPR,
          Spain’s LOPDGDD, California’s CCPA, or any other relevant jurisdictional law). In
          practice, Novik’s platform is designed so that Users should only provide anonymized data,
          meaning that Novik generally processes data that is not identifiable. If, however, any
          personal data is processed, Novik will ensure all required legal bases for processing are
          present (e.g., User’s consent or other lawful grounds), and that appropriate security and
          confidentiality measures are in place as required by law. (For further information on how
          personal data is handled, refer to our Privacy Policy.)
        </p>
        <p>
          Multijurisdictional Compliance: Novik operates globally and strives to follow best
          practices in data privacy across different jurisdictions. For example, the table below
          outlines some key differences between the European Union’s General Data Protection
          Regulation (GDPR) and the California Consumer Privacy Act (CCPA), which are two prominent
          data protection laws that reflect the evolving standards in their respective regions. This
          comparison highlights how Novik aligns its data handling with these frameworks:
        </p>
        <p>
          In essence, Novik adheres to the principle of applying the most protective standards among
          the applicable laws to safeguard user data. When Users provide data (even if anonymized as
          expected), Novik’s internal policies ensure compliance with global data protection
          requirements such as GDPR and CCPA, among others. This means respecting rights of access
          and deletion, maintaining transparency about data use, and ensuring data security. For
          detailed information on how Novik handles personal information and to learn about user
          privacy rights, please see our Privacy Policy
        </p>
        <h2 className="text-xl font-semibold mt-6">XI. SPECIFIC PROHIBITIONS FOR USERS</h2>
        <p>
          When using Novik’s Services, the User agrees NOT to engage in any conduct that could
          damage, disable, overburden, or impair the platform or interfere with any other party’s
          use of the Services. In particular, and without limitation, the User shall not:
        </p>
        <ul>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Upload Unauthorized Content:</strong>Upload, transmit,
            or share via Novik any content for which the User does not have the necessary rights or
            permissions, or which infringes any intellectual property rights, copyrights,
            trademarks, trade secrets or privacy/publicity rights of any third party. This includes,
            for example, clinical images, documents, or other material that the User is not
            authorized to disclose or that violate confidentiality obligations. The User is
            responsible for ensuring that any case data or materials uploaded either belong to the
            User or the User has obtained all relevant consents and authorizations to use them on
            the platform.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>System Manipulation:</strong> Use automated scripts,
            bots, scrapers, or other technological tools to access Novik’s Services, collect data
            from the platform, or otherwise interact with or monitor the Services in a manner not
            authorized by Novik. Any attempt to interfere with the normal functioning of the
            platform, including by hacking, decrypting, injecting malicious code, or other illicit
            means, is strictly forbidden. The User shall not attempt to circumvent any content
            filtering techniques or security measures that Novik employs, nor attempt to access
            areas or features of the Service for which they have no authorization.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Distribution of Malicious Content:</strong> Upload,
            share, or otherwise introduce any viruses, worms, malware, spyware, Trojan horses, or
            any other malicious or harmful software or content that could disrupt the proper
            functioning of the platform or compromise its security or that of Users. The User must
            refrain from using the Service to transmit unsolicited mass communications (spam),
            phishing content, or any material that contains harmful or deleterious programs.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Unauthorized Commercial Use:</strong> Exploit the
            platform or any portion of its Content for any commercial purposenot expressly permitted
            by Novik. This includes, without limitation: selling, reselling or renting the Service
            or access to the Service; using the output of Novik (such as AI-generated reports or
            answers) as part of a for-profit service to third parties without Novik’s consent; or
            any form of unauthorized advertising, marketing, or promotional activity through the
            platform (e.g., posting advertisements in query inputs or outputs, or using Novik’s
            responses to generate content for commercial distribution without a license). The
            Service is provided for the User’s own professional use with their patients or practice;
            any broader commercial utilization requires a separate agreement.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Impersonation:</strong> Impersonate or misrepresent
            affiliation with any person or entity. The User shall not pretend to be another user, a
            Novik administrator, an employee of BlueBioPlan, a medical expert that they are not, or
            any other person. Similarly, the User must not falsely state or imply a relationship
            with Novik or BlueBioPlan (for instance, claiming to be an “official partner” or
            “certified” by Novik) without prior written authorization. This prohibition extends to
            the misuse of credentials — for example, sharing login details to let unauthorized
            persons access the platform under the User’s account is not allowed.
          </li>
          <li>
            &nbsp;&nbsp;&nbsp;&nbsp;<strong>Unauthorized Data Collection:</strong> Collect or
            harvest any personally identifiable information from other users of the Service or from
            the platform without proper consent. This includes attempting to access or search the
            Service by any means (automated or otherwise) other than through our currently
            available, published interfaces that are provided by Novik. Users must not probe, scan,
            or test the vulnerability of any Novik system or network, or breach any security or
            authentication measures. Gathering data such as other users’ content, profiles, or usage
            information through techniques like scraping or via any form of packet sniffer or
            network monitoring tool is strictly prohibited.
          </li>
        </ul>
        <p>
          Violation of the above prohibitions constitutes a material breach of these Terms and may
          result in immediate suspension or termination of the User’s account and access to the
          Services. Additionally, Novik (BlueBioPlan S.L.U.) reserves the right to take appropriate
          legal action against the User for any illegal or unauthorized use of the platform,
          including seeking injunctive relief, damages, and reporting to law enforcement authorities
          if applicable. The User understands that Novik may monitor compliance with these Terms
          (subject to Novik’s Privacy Policy and applicable laws) and that Novik has the right to
          investigate any suspected violation.
        </p>
        <h2 className="text-xl font-semibold mt-6">XII. LANGUAGE AND PREVAILING VERSION</h2>
        <p>
          This Legal Notice and Terms of Use is available in Spanish and may be provided in English
          (or other languages) for the convenience of international users. In the event of any
          discrepancy or conflict between the Spanish version and an English (translated) version,
          the Spanish version shall prevail and be considered the authoritative text for legal
          purposes. The English text is provided as a courtesy translation. By agreeing to these
          Terms, the User acknowledges that the original Terms were drafted in Spanish and that
          Spanish law is the governing law as specified above.
        </p>
      </div>
    ),
    privacy: (
      <div className="space-y-4 text-gray-700 text-md leading-relaxed">
        <h2>www.novik.ai</h2>

        <h2>I. Privacy Policy and Data Protection</h2>
        <p>
          In compliance with current legislation, <strong>Novik</strong> (hereinafter also referred
          to as the Website) commits to adopting the necessary technical and organizational measures
          according to the appropriate security level for the risk of the data collected.
        </p>

        <h2>Laws Incorporated in this Privacy Policy</h2>
        <p>
          This privacy policy is adapted to the current Spanish and European regulations on personal
          data protection on the internet. Specifically, it complies with the following regulations:
        </p>
        <ul>
          <li>
            <strong>Regulation (EU) 2016/679</strong> (GDPR)
          </li>
          <li>
            <strong>Organic Law 3/2018</strong> (LOPD-GDD)
          </li>
          <li>
            <strong>Royal Decree 1720/2007</strong> (RDLOPD)
          </li>
          <li>
            <strong>Law 34/2002</strong> (LSSI-CE)
          </li>
        </ul>

        <h2>Identity of the Data Controller</h2>
        <p>
          The data controller is <strong>BlueBioPlan SL</strong>, with NIF: B-72913379.
        </p>
        <p>
          <strong>Address:</strong> C/ La Cruz 2 Entlo. Alcantarilla 30820 Murcia (Spain)
        </p>
        <p>
          <strong>Contact phone:</strong> +34 690 957 910
        </p>
        <p>
          <strong>Contact email:</strong> info@novik.ai
        </p>

        <h2>Personal Data Record</h2>
        <p>
          In compliance with the GDPR and LOPD-GDD, we inform you that the personal data collected
          by Novik through the forms on its pages will be included and processed in our file for the
          purpose of facilitating, speeding up, and fulfilling the commitments established between
          Novik and the User or maintaining the relationship established in the forms that they
          complete, or to address a request or inquiry.
        </p>
        <p>
          Unless the exception provided in Article 30.5 of the GDPR applies, a record of processing
          activities is maintained that specifies, according to its purposes, the processing
          activities carried out and other circumstances established by the GDPR.
        </p>

        <h2>Principles Applicable to the Processing of Personal Data</h2>
        <p>
          The processing of User's personal data will be subject to the following principles as
          stated in Article 5 of the GDPR and Article 4 and following of Organic Law 3/2018:
        </p>
        <ul>
          <li>
            <strong>Lawfulness, fairness, and transparency:</strong> User consent will be required
            at all times with complete transparency regarding the purposes.
          </li>
          <li>
            <strong>Purpose limitation:</strong> Data collected will be for specific, explicit, and
            legitimate purposes.
          </li>
          <li>
            <strong>Data minimization:</strong> Only necessary data will be collected.
          </li>
          <li>
            <strong>Accuracy:</strong> Data must be accurate and kept up to date.
          </li>
          <li>
            <strong>Storage limitation:</strong> Data will only be kept as long as necessary for
            processing purposes.
          </li>
          <li>
            <strong>Integrity and confidentiality:</strong> Data will be processed securely and
            confidentially.
          </li>
          <li>
            <strong>Proactive accountability:</strong> The controller will ensure compliance with
            these principles.
          </li>
        </ul>

        <h2>Categories of Personal Data</h2>
        <p>
          The categories of data processed by Novik are solely identifying data. No special
          categories of data as defined in Article 9 of the GDPR are processed.
        </p>
        <h2>Legal Basis for Processing Personal Data</h2>
        <p>
          The legal basis for the processing of personal data is <strong>consent</strong>. Novik
          commits to obtaining the express and verifiable consent of the User for the processing of
          their personal data for one or more specific purposes.
        </p>
        <p>
          The User will have the right to <strong>withdraw their consent</strong> at any time.
          Withdrawal of consent will be as easy as giving it. As a general rule, the withdrawal of
          consent will not condition the use of the Website.
        </p>
        <p>
          Whenever the User must or can provide their data through forms to make inquiries, request
          information, or for reasons related to the content of the Website, they will be informed
          if the completion of any of them is mandatory due to their being essential for the correct
          execution of the operation performed.
        </p>

        <h2>Purposes of Processing Personal Data</h2>
        <p>
          Personal data is collected and managed by Novik to facilitate, expedite, and fulfill the
          commitments established between the Website and the User or to maintain the relationship
          established in the forms that they complete or to address a request or inquiry.
        </p>
        <p>
          Likewise, the data may be used for{' '}
          <strong>
            commercial purposes, personalization, operational and statistical purposes
          </strong>
          , and activities related to Novik's corporate purpose, as well as for the extraction,
          storage of data, and marketing studies to tailor the Content offered to the User, as well
          as to improve the quality, operation, and navigation of the Website.
        </p>
        <p>
          At the time personal data is collected, the User will be informed of the{' '}
          <strong>specific purpose(s)</strong> of the processing to which the personal data will be
          used.
        </p>

        <h2>Periods of Retention of Personal Data</h2>
        <p>
          Personal data will only be retained for the minimum time necessary for the purposes of its
          processing and, in any case, only for the following period: <strong>36 months</strong>, or
          until the User requests its deletion.
        </p>
        <p>
          At the time personal data is collected, the User will be informed of the period during
          which personal data will be retained or, when that is not possible, the criteria used to
          determine this period.
        </p>

        <h2>Recipients of Personal Data</h2>
        <p>
          The User's personal data will be shared with the following recipients or categories of
          recipients:
        </p>
        <ul>
          <li>Google Analytics</li>
          <li>Other partner companies</li>
        </ul>
        <p>
          If the Data Controller intends to transfer personal data to a third country or
          international organization, the User will be informed at the time the personal data is
          obtained about the third country or international organization to which the data is
          intended to be transferred, as well as the existence or absence of an adequacy decision by
          the Commission.
        </p>

        <h2>Personal Data of Minors</h2>
        <p>
          In accordance with Articles 8 of the GDPR and 7 of Organic Law 3/2018, only those over the
          age of <strong>14</strong> may lawfully consent to the processing of their personal data
          by Novik.
        </p>
        <p>
          In the case of minors under 14 years of age, the{' '}
          <strong>consent of parents or guardians</strong> will be required, and the processing will
          only be lawful to the extent that they have authorized it.
        </p>

        <h2>Confidentiality and Security of Personal Data</h2>
        <p>
          Novik commits to adopting the necessary technical and organizational measures, according
          to the appropriate security level for the risk of the data collected, to ensure the{' '}
          <strong>security of personal data</strong> and prevent accidental or unlawful destruction,
          loss, or alteration of personal data transmitted, stored, or otherwise processed, or
          unauthorized disclosure of or access to such data.
        </p>
        <p>
          The Website is equipped with an <strong>SSL (Secure Socket Layer) certificate</strong>,
          which ensures that personal data is transmitted securely and confidentially.
        </p>
        <p>
          However, since Novik cannot guarantee the invulnerability of the internet, the Data
          Controller commits to <strong>notifying the User without undue delay</strong> when a
          personal data breach occurs that is likely to result in a high risk to the rights and
          freedoms of natural persons.
        </p>
        <p>
          Personal data will be treated as <strong>confidential</strong> by the Data Controller, who
          commits to ensuring and guaranteeing by means of a legal or contractual obligation that
          such confidentiality is respected by its employees, partners, and any person to whom the
          information is made accessible.
        </p>
        <h2>Rights Derived from the Processing of Personal Data</h2>
        <p>
          The User has the following rights over Novik and may therefore exercise the following
          rights against the Data Controller as recognized by the GDPR and Organic Law 3/2018, of 5
          December, on Personal Data Protection and Guarantee of Digital Rights:
        </p>
        <ul>
          <li>
            <strong>Right of access:</strong> The right of the User to obtain confirmation of
            whether Novik is processing their personal data or not, and if so, to obtain information
            about their specific personal data and the processing that Novik has carried out or is
            carrying out.
          </li>
          <li>
            <strong>Right of rectification:</strong> The right of the User to have their personal
            data modified if it proves to be inaccurate or incomplete.
          </li>
          <li>
            <strong>Right of erasure ("the right to be forgotten"):</strong> The right to obtain the
            erasure of personal data under certain conditions, including when it's no longer needed,
            unlawfully processed, or consent is withdrawn.
          </li>
          <li>
            <strong>Right to restriction of processing:</strong> The right to limit processing under
            specific scenarios such as when accuracy is contested or processing is unlawful.
          </li>
          <li>
            <strong>Right to data portability:</strong> The right to receive personal data in a
            structured, machine-readable format and transmit it to another controller, when
            processing is automated.
          </li>
          <li>
            <strong>Right of opposition:</strong> The right to object to the processing of personal
            data or request cessation.
          </li>
          <li>
            <strong>Right not to be subject to automated decision-making:</strong> The right not to
            be subject to decisions based solely on automated processing, including profiling.
          </li>
        </ul>

        <p>
          The User may exercise these rights by written communication addressed to the Data
          Controller with the reference <strong>"GDPR-www.novik.ai"</strong>, specifying the
          following:
        </p>
        <ul>
          <li>
            Name, surname of the User, and a copy of the DNI or legal ID. If represented, the
            representative's ID and proof of representation are required.
          </li>
          <li>A specific request outlining the information or rights to be exercised.</li>
          <li>Address for notifications.</li>
          <li>Date and signature.</li>
          <li>Any supporting documents.</li>
        </ul>

        <p>This request and any other attached documents may be sent to:</p>
        <ul>
          <li>
            <strong>Postal address:</strong> C/ La Cruz 2 Entlo. Alcantarilla 30820 Murcia (Spain)
          </li>
          <li>
            <strong>Email:</strong> info@novik.ai
          </li>
        </ul>

        <h2>Links to Third-Party Websites</h2>
        <p>
          The Website may include hyperlinks to third-party websites that are not operated by Novik.
          The owners of those websites have their own data protection policies and are solely
          responsible for their content and practices.
        </p>
        <h2>Claims Before the Supervisory Authority</h2>
        <p>
          In the event that the User considers that there is a problem or violation of the current
          regulations in the way their personal data is being processed, they will have the right to
          effective judicial protection and to file a complaint with a supervisory authority,
          particularly in the State in which they have their habitual residence, place of work, or
          place of the alleged violation.
        </p>
        <p>
          In the case of Spain, the supervisory authority is the
          <a href="https://www.aepd.es/" target="_blank" rel="noopener noreferrer">
            Spanish Data Protection Agency
          </a>
          .
        </p>

        <h2>II. Acceptance and Changes to This Privacy Policy</h2>
        <p>
          The User must have read and be in agreement with the conditions on personal data
          protection contained in this Privacy Policy, as well as accept the processing of their
          personal data so that the Data Controller can proceed with it in the manner, during the
          periods, and for the purposes stated. The use of the Website will imply acceptance of its
          Privacy Policy.
        </p>
        <p>
          Novik reserves the right to modify its Privacy Policy in accordance with its own criteria,
          or due to a legislative, jurisprudential, or doctrinal change by the Spanish Data
          Protection Agency. Changes or updates to this Privacy Policy will not be explicitly
          notified to the User. It is recommended that the User consult this page periodically to
          stay informed of the latest changes or updates.
        </p>
        <p>
          This Privacy Policy was updated to comply with Regulation (EU) 2016/679 and Organic Law
          3/2018, on Personal Data Protection and Guarantee of Digital Rights.
        </p>

        <h2>Annex 1: Collection of Registration Data</h2>
        <p>
          <strong>Personal Data Collected:</strong>
        </p>
        <ul>
          <li>
            <strong>Email:</strong> Required for the creation and management of the user's account,
            and for communication.
          </li>
          <li>
            <strong>Country:</strong> Collected for statistical purposes and personalization based
            on local regulations.
          </li>
        </ul>

        <p>
          <strong>Purpose of Processing:</strong>
        </p>
        <ul>
          <li>Managing user accounts on the platform.</li>
          <li>Sending notifications about updates, events, or changes to the platform.</li>
          <li>Responding to inquiries and providing technical support.</li>
          <li>Personalizing user experience by region.</li>
          <li>Conducting usage statistics for service improvement.</li>
        </ul>

        <h2>Annex 2: Processing of Anonymized Patient Data</h2>
        <p>
          <strong>Data Anonymization:</strong>
        </p>
        <ul>
          <li>Novik only processes anonymized patient data.</li>
          <li>
            Users are responsible for ensuring data is anonymized in line with applicable
            regulations (e.g., GDPR).
          </li>
        </ul>

        <p>
          <strong>Purpose of Processing:</strong>
        </p>
        <ul>
          <li>Generate clinical recommendations based on the provided medical information.</li>
          <li>
            Data is not shared with third parties or used for commercial or advertising purposes
            outside Novik's scope.
          </li>
        </ul>

        <h2>Annex 3: Cookies and Tracking</h2>
        <p>
          <strong>Use of Cookies:</strong>
        </p>
        <ul>
          <li>Enhance the user experience.</li>
          <li>Store user preferences.</li>
          <li>Analyze website traffic.</li>
          <li>Remember the user's session.</li>
          <li>Provide personalized content.</li>
        </ul>
        <h2>Cookie Management</h2>
        <ul>
          <li>
            Users can configure their browser to reject cookies. However, certain platform
            functionalities may be limited if they are disabled.
          </li>
        </ul>

        <h2>Annex 4: User Rights</h2>
        <p>
          <strong>Access and Rectification:</strong>
        </p>
        <ul>
          <li>
            Users have the right to access and modify the information provided during registration
            (such as their email and country).
          </li>
        </ul>
        <p>
          <strong>Account Deletion:</strong>
        </p>
        <ul>
          <li>
            Users may request the deletion of their account and personal data stored on the platform
            at any time. After deletion, only anonymized data necessary for statistical purposes or
            to comply with legal obligations will be retained.
          </li>
        </ul>

        <h2>Annex 5: International Data Transfers</h2>
        <p>
          <strong>Data Transfers Outside the EU:</strong>
        </p>
        <ul>
          <li>
            If Novik's servers or service providers are located outside the European Union, it will
            be ensured that the transfer of personal data is carried out in compliance with
            appropriate protective measures, such as the GDPR's standard contractual clauses.
          </li>
        </ul>

        <h2>Annex 6: Data Sharing with Third Parties</h2>

        <p>
          <strong>1. Sharing Personal Data with Third Parties:</strong>
        </p>
        <ul>
          <li>
            Novik does not share, sell, or rent users' personal data (such as email or country) to
            third parties for commercial purposes without the user's explicit consent.
          </li>
          <li>
            Personal data will only be shared with third parties in the following cases:
            <ul>
              <li>
                To comply with a legal obligation, such as a court order or request from a competent
                authority.
              </li>
              <li>
                To provide the contracted service by the user, in which case the data will be shared
                with service providers that assist Novik in developing and managing the platform.
                These service providers will be subject to the same confidentiality and data
                protection obligations as Novik.
              </li>
              <li>
                If the user provides their explicit consent, in which case the specific purpose of
                sharing the data with third parties will be detailed.
              </li>
            </ul>
          </li>
        </ul>

        <p>
          <strong>2. Third-Party Processing and Contractual Obligations:</strong>
        </p>
        <ul>
          <li>
            If Novik contracts third parties for the processing of personal data, such as hosting
            providers or technical services, these third parties will be contractually obligated to:
            <ul>
              <li>Use the personal data only for the purposes established by Novik.</li>
              <li>
                Protect the confidentiality and security of personal data in accordance with the
                GDPR or any applicable regulation.
              </li>
              <li>
                Not use or share the data for other purposes without the user's explicit consent.
              </li>
            </ul>
          </li>
        </ul>

        <p>
          <strong>3. Sharing Anonymized Data:</strong>
        </p>
        <ul>
          <li>
            Novik may share anonymized data with third parties for statistical, research, or service
            improvement purposes, provided that:
            <ul>
              <li>This data cannot be used to identify any person.</li>
              <li>
                It is ensured that the data processing complies with the highest standards of
                security and anonymization.
              </li>
              <li>
                Third parties receiving this data may not, under any circumstances, attempt to
                re-identify individuals from the information received.
              </li>
            </ul>
          </li>
        </ul>

        <p>
          <strong>4. International Data Transfers:</strong>
        </p>
        <ul>
          <li>
            If users' personal data or anonymized data are transferred to third parties outside the
            EU, Novik will ensure that such transfers comply with applicable data protection
            regulations by using tools such as:
            <ul>
              <li>Standard contractual clauses approved by the European Commission.</li>
            </ul>
          </li>
        </ul>
        <ul>
          <li>
            <strong>Certifications under international frameworks:</strong> Novik ensures compliance
            with adequate data protection standards for international transfers by relying on
            certifications such as the Privacy Shield (where applicable).
          </li>
        </ul>

        <h2>5. Transparency and User Consent</h2>
        <ul>
          <li>
            <strong>Transparency:</strong> Novik will maintain complete transparency regarding any
            sharing of personal data with third parties. Users will be informed before any data
            sharing not covered by this privacy policy.
          </li>
          <li>
            <strong>User Consent:</strong> The User may revoke their consent for sharing their data
            with third parties at any time, except in cases where it is necessary to comply with
            legal or contractual obligations.
          </li>
        </ul>
      </div>
    ),
    cookies: <p>Coming soon: Cookie Policy</p>,
  };

  return (
    <div className="bg-dental w-screen h-screen flex flex-col">
      <Header />
      <div className="flex overflow-y-auto">
        <aside className="w-64 border-r p-4 space-y-2">
          {sections.map(({ id, label }) => (
            <button
              key={id}
              onClick={() => setActive(id)}
              className={`block w-full text-left p-2 rounded ${
                active === id ? 'bg-gray-100 font-semibold text-blue-600' : ''
              }`}
            >
              {label}
            </button>
          ))}
        </aside>

        <main className="flex-1 p-8 overflow-y-auto shadow-inner">
          <h1 className="text-2xl font-bold mb-4">{sections.find(s => s.id === active)?.label}</h1>
          <div className="text-gray-700 leading-relaxed">
            {content[active as keyof typeof content]}
          </div>
        </main>
      </div>
    </div>
  );
};

export default LegalPage;
