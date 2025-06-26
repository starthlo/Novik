import { useState } from 'react';
import Header from './Common/Header';
import PoweredBy from './Common/PoweredBy';
import Partners from './Common/Partners';
import Footer from './Common/Footer';
import NovikLogo from '../assets/Novik.png';

interface FormData {
  name: string;
  email: string;
  category: string;
  phone: string;
  message: string;
}

function ContactPage() {
  const [formData, setFormData] = useState<FormData>({
    name: '',
    email: '',
    category: '',
    phone: '',
    message: '',
  });

  const handleChange = (
    e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement | HTMLTextAreaElement>
  ) => {
    const { name, value } = e.target;
    setFormData(prev => ({
      ...prev,
      [name]: value,
    }));
  };

  const handleSubmit = (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    console.log(formData);
    alert('Form submitted successfully!');
  };

  return (
    <div className="bg-dental w-screen h-screen flex flex-col">
      <Header />
      <div className="flex-grow overflow-y-auto">
        <div className="text-center m-16 flex flex-col items-center">
          <img src={NovikLogo} alt="Novik" className="h-20 mb-4" />
          <p className="text-2xl text-gray-500">
            Your smart AI Dental assistant for safe clinical decisions
          </p>
        </div>
        <div className="bg-white p-10 shadow-md rounded-lg w-[90%] max-w-[1000px] mx-auto mb-16">
          <h2 className="text-4xl font-bold text-gray-500 mb-6">We want to hear from you!</h2>
          <p className="text-gray-500 mb-4">
            We're thrilled you made it this far, because it only means one thing: you want to get in
            touch with us!
          </p>
          <p className="text-gray-500 mb-4">
            Whether you've got a question, a suggestion, a crazy idea, or just want to say hello
            (and hey, we love that too), no one here will judge your curiosity, quite the opposite,
            we're huge fans of it.
          </p>
          <p className="text-gray-500 mb-4">
            If you're looking for answers, a collaboration, or even just a friendly chat, we're here
            for you.
          </p>
          <p className="text-gray-500 mb-4">
            We promise not to be the type who hides behind endless forms or leaves you on "read." We
            love conversation as much as we love coffee (and that's saying a lot).
          </p>
          <p className="text-gray-500 mb-4">
            Let's make this easy: choose how you want to reach us and send your message our way.
          </p>
          <p className="text-gray-500 mb-4">
            Whether by email, a call, or even smoke signals (well, maybe not smoke signals, but you
            get the idea).
          </p>
          <p className="text-gray-500 mb-4">
            We're here to listen, to help, and if needed, to throw in a bad joke or two.
          </p>
          <p className="text-gray-500 mb-4">
            So go ahead, don't be shy. Reach out, we're just a click away and would love to know
            what's on your mind.
          </p>
          <h2 className="text-4xl font-bold text-gray-500 mb-6">Contact US!</h2>
          <form onSubmit={handleSubmit} className="space-y-4">
            <div className="flex gap-4">
              <input
                type="text"
                name="name"
                value={formData.name}
                onChange={handleChange}
                placeholder="Name"
                className="w-full border border-gray-400 rounded-md p-2"
                required
              />
              <input
                type="email"
                name="email"
                value={formData.email}
                onChange={handleChange}
                placeholder="Email"
                className="w-full border border-gray-400 rounded-md p-2"
                required
              />
            </div>

            <select
              name="category"
              value={formData.category}
              onChange={handleChange}
              className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
              required
            >
              <option value="">Select Category</option>
              <option value="Implantology">General Dentistry</option>
              <option value="Orthodontics">Endodontics</option>
              <option value="Endodontics">Orthodontics</option>
              <option value="Implantology">Oral and Maxillofacial Surgery</option>
              <option value="Orthodontics">Oral Implantology</option>
              <option value="Endodontics">Pediatric Dentistry</option>
              <option value="Implantology">Prosthodontics / Oral Rehabilitation</option>
              <option value="Orthodontics">Cosmetic Dentistry</option>
              <option value="Endodontics">Oral and Maxillofacial Pathology</option>
              <option value="Implantology">Oral Medicine</option>
              <option value="Orthodontics">Oral and Maxillofacial Radiology</option>
              <option value="Endodontics">Preventive and Community Dentistry</option>
              <option value="Implantology">Forensic Dentistry</option>
              <option value="Orthodontics">Dental Sleep Medicine</option>
              <option value="Endodontics">Geriatric Dentistry</option>
              <option value="Implantology">Restorative Dentistry</option>
              <option value="Orthodontics">Digital Dentistry / CAD-CAM</option>
              <option value="Endodontics">Minimally Invasive Dentistry</option>
              <option value="Endodontics">Dental Biomaterials and Bioengineering</option>
              <option value="Implantology">Dental Education and Research</option>
              <option value="Orthodontics">Practice Management and Dental Administration</option>
              <option value="Endodontics">Student</option>
            </select>

            <input
              type="text"
              name="phone"
              value={formData.phone}
              onChange={handleChange}
              placeholder="Phone"
              className="w-full border border-gray-400 rounded-md p-2"
              required
            />

            <textarea
              name="message"
              value={formData.message}
              onChange={handleChange}
              placeholder="Message"
              rows={4}
              className="w-full border border-gray-400 rounded-md p-2"
              required
            />

            <button
              type="submit"
              className="w-full bg-orange-500 text-white py-2 rounded-md hover:bg-orange-600 transition"
            >
              Send
            </button>
          </form>
        </div>
        <PoweredBy />
        <Partners />
        <Footer />
      </div>
    </div>
  );
}

export default ContactPage;
