import {
  FaEnvelope,
  FaPhone,
  FaWhatsapp,
  FaFacebookF,
  FaInstagram,
  FaLinkedinIn,
} from "react-icons/fa";
import { Link } from "react-router-dom";

function Footer() {
  return (
    <footer className="bg-white text-gray-500 border-t pt-8 pb-4 text-sm mt-0">
      <div className="max-w-7xl mx-auto px-4 grid grid-cols-1 md:grid-cols-3 gap-8 items-start">
        <div>
          <h2 className="text-2xl font-semibold">
            Novik<span className="text-orange-600">.</span>
          </h2>
          <p className="mt-2 text-base">© 2024 Novik. All rights reserved.</p>
          <div className="flex flex-wrap gap-2 mt-2 text-base text-orange-500">
            <Link to="/legal#terms">Terms of Service</Link> |
            <Link to="/legal#privacy">Privacy Policy</Link> |
            <Link to="/legal#cookies">Cookie Policy</Link>
          </div>
        </div>

        <div className="flex flex-col gap-4 text-center md:text-left">
          <div>
            <h3 className="text-base font-semibold">Contact Us</h3>
            <div className="flex justify-center md:justify-start gap-4 text-orange-500 mt-2">
              <a href="mailto:info@novik.ai" title="Email">
                <FaEnvelope className="text-xl hover:scale-110 transition" />
              </a>
              <a href="tel:+34690957910" title="Call">
                <FaPhone className="text-xl hover:scale-110 transition" />
              </a>
              <a
                href="https://wa.me/34690957910"
                target="_blank"
                rel="noopener noreferrer"
                title="WhatsApp"
              >
                <FaWhatsapp className="text-xl hover:scale-110 transition" />
              </a>
            </div>
          </div>
          <div>
            <h3 className="text-base font-semibold">Social Links</h3>
            <div className="flex justify-center md:justify-start gap-4 text-orange-500 mt-2">
              <a
                href="https://www.facebook.com/profile.php?id=61567745501156"
                target="_blank"
                rel="noopener noreferrer"
                title="Facebook"
              >
                <FaFacebookF className="text-xl hover:scale-110 transition" />
              </a>
              <a
                href="https://www.instagram.com/dentalnovik/"
                target="_blank"
                rel="noopener noreferrer"
                title="Instagram"
              >
                <FaInstagram className="text-xl hover:scale-110 transition" />
              </a>
              <a
                href="https://www.linkedin.com/company/novik-ai"
                target="_blank"
                rel="noopener noreferrer"
                title="LinkedIn"
              >
                <FaLinkedinIn className="text-xl hover:scale-110 transition" />
              </a>
            </div>
          </div>
        </div>
      </div>

      <p className="max-w-7xl text-center md:text-left mx-auto mt-6 px-4 text-base text-gray-500">
        Novik is an experimental technology demonstrator. Novik does not provide
        medical advice, diagnosis or treatment. User questions and other inputs
        on Novik are not covered by HIPAA. It is the responsibility of the user
        to ensure questions do not contain protected health information (PHI) or
        any information that violates the privacy of any person.
      </p>
    </footer>
  );
}

export default Footer;
