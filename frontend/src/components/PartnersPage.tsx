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

function PartnersPage() {
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
          <p className="text-2xl text-gray-700">
            Your smart AI Dental assistant for safe clinical decisions
          </p>
        </div>
        <div className="bg-white p-10 shadow-md rounded-lg w-[90%] max-w-[1000px] mx-auto mb-16">
          <h2 className="text-4xl font-bold text-gray-500 mb-6">
            Hello, potential partners and allies!
          </h2>
          <p className="text-gray-500 mb-4">
            We're thrilled to have you here because it means something about Novik has already
            caught your attention.
          </p>
          <p className="text-gray-500 mb-4">
            And let me tell you the excitement is mutual! We're always on the lookout for
            connections that feel more like partnerships, something that just clicks, like it was
            meant to happen (yeah, we get a little romantic about it, what can we say?).
          </p>
          <p className="text-gray-500 mb-4">
            We know sponsoring a project is a big decision, and here at Novik, we believe the best
            decisions are made with a smile.
          </p>
          <p className="text-gray-500 mb-4">
            Our mission is to ensure that everyone, from dentists to students, feels supported and
            guided, and we believe that together we can take this experience to the next level.
          </p>
          <p className="text-gray-500 mb-6">
            So go ahead, if you think Novik could be a great place for your support, we're just a
            message away.
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
              className="w-full border border-gray-400 rounded-md p-2"
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

export default PartnersPage;
