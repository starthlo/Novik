import { useEffect, useState } from 'react';
import Header from './Common/Header';
import NovikLogo from '../assets/Novik.png';
import { Country, State, City, ICountry, IState, ICity } from 'country-state-city';

function Register() {
  const [formData, setFormData] = useState({
    username: '',
    email: '',
    password: '',
    confirm_password: '',
    dob: '',
    phone: '',
    occupation: '',
    country: '',
    state: '',
    city: '',
    agreeToTerms: false,
    receiveInfo: false,
  });
  const [countries] = useState<ICountry[]>(Country.getAllCountries());
  const [states, setStates] = useState<IState[]>([]);
  const [cities, setCities] = useState<ICity[]>([]);

  useEffect(() => {
    if (formData.country) {
      const selected = countries.find(c => c.name === formData.country);
      if (selected) {
        const stateList = State.getStatesOfCountry(selected.isoCode);
        setStates(stateList);
      }
    }
  }, [formData.country]);

  useEffect(() => {
    const countryObj = countries.find(c => c.name === formData.country);
    const stateObj = states.find(s => s.name === formData.state);
    if (countryObj && stateObj) {
      const cityList = City.getCitiesOfState(countryObj.isoCode, stateObj.isoCode);
      setCities(cityList);
    }
  }, [formData.state]);

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    const { name, value, type } = e.target;

    if (type === 'checkbox') {
      setFormData(prev => ({
        ...prev,
        [name]: (e.target as HTMLInputElement).checked,
      }));
    } else {
      setFormData(prev => ({
        ...prev,
        [name]: value,
      }));
    }
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    try {
      // console.log(formData)
      const response = await fetch('/api/register/', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(formData),
      });
      // console.log(response)
      if (!response.ok) {
        const error = await response.json();
        // alert(error.error);
        const firstKey = Object.keys(error)[0];
        const firstMessage = error[firstKey][0];
        alert(firstMessage);
      } else {
        const result = await response.json();
        console.log(result);
        alert('Registration successful!');
      }
    } catch (error) {
      console.error('Registration failed:', error);
      alert('An error occurred during registration.');
    }
  };

  return (
    <div className="bg-dental w-screen h-screen flex flex-col">
      <Header />
      <div className="flex overflow-y-auto flex-col items-center pt-12">
        <img src={NovikLogo} alt="Novik" className="h-20 mb-8" />
        <h1 className="text-xl text-gray-700 mb-12 text-center">
          Your smart AI Dental assistant for safe clinical decisions
        </h1>
        <div className="bg-white p-8 rounded-lg shadow-md w-[600px] mb-8">
          <h2 className="text-3xl text-gray-500 font-bold text-center mb-4">User Registration</h2>
          <p className="text-gray-500 text-center mb-6">
            To use the NOVIK Dental Assistant you need to register. Once you have registered you
            will be able to use Novik and discover its full potential.
          </p>

          <form onSubmit={handleSubmit} className="space-y-4">
            <div className="flex gap-4">
              <input
                type="text"
                name="username"
                value={formData.username}
                onChange={handleChange}
                placeholder="Name"
                className="w-full border border-gray-400 rounded-md p-2"
              />
              <input
                type="date"
                name="dob"
                value={formData.dob}
                onChange={handleChange}
                placeholder="dd--yyyy"
                className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
              />
            </div>

            <div className="flex gap-4">
              <input
                type="email"
                name="email"
                value={formData.email}
                onChange={handleChange}
                placeholder="Email"
                className="w-full border border-gray-400 rounded-md p-2 bg-gray-100"
              />
              <input
                type="text"
                name="phone"
                value={formData.phone}
                onChange={handleChange}
                placeholder="Phone"
                className="w-full border border-gray-400 rounded-md p-2"
              />
            </div>

            <div className="flex gap-4">
              <input
                type="password"
                name="password"
                value={formData.password}
                onChange={handleChange}
                placeholder="Password"
                className="w-full border border-gray-400 rounded-md p-2"
              />
              <input
                type="password"
                name="confirm_password"
                value={formData.confirm_password}
                onChange={handleChange}
                placeholder="Confirm Password"
                className="w-full border border-gray-400 rounded-md p-2"
              />
            </div>

            <select
              name="occupation"
              value={formData.occupation}
              onChange={handleChange}
              className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
            >
              <option value="">Select Occupation</option>
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

            <select
              name="country"
              value={formData.country}
              onChange={handleChange}
              className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
            >
              <option value="">Select Country</option>
              {countries.map(c => (
                <option key={c.isoCode} value={c.name}>
                  {c.name}
                </option>
              ))}
            </select>

            <select
              name="state"
              value={formData.state}
              onChange={handleChange}
              className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
            >
              <option value="">Select State</option>
              {states.map(s => (
                <option key={s.isoCode} value={s.name}>
                  {s.name}
                </option>
              ))}
            </select>

            <select
              name="city"
              value={formData.city}
              onChange={handleChange}
              className="w-full border border-gray-400 text-gray-500 rounded-md p-2"
            >
              <option value="">Select City</option>
              {cities.map(city => (
                <option key={city.name} value={city.name}>
                  {city.name}
                </option>
              ))}
            </select>

            <div className="space-y-2">
              <label className="flex items-center text-gray-500 gap-2">
                <input
                  type="checkbox"
                  name="agreeToTerms"
                  checked={formData.agreeToTerms}
                  onChange={handleChange}
                />
                <span>
                  I have read and agree to the{' '}
                  <a href="/terms" className="text-orange-500 hover:underline">
                    Legal Notice
                  </a>{' '}
                  and{' '}
                  <a href="/privacy" className="text-orange-500 hover:underline">
                    Privacy Policy
                  </a>
                </span>
              </label>

              <label className="flex items-center text-gray-500 gap-2">
                <input
                  type="checkbox"
                  name="receiveInfo"
                  checked={formData.receiveInfo}
                  onChange={handleChange}
                />
                <span>
                  I would like to receive commercial information about the product or service
                </span>
              </label>
            </div>

            <button
              type="submit"
              className="w-full bg-orange-500 text-white py-2 rounded-md hover:bg-orange-600 transition"
            >
              Submit
            </button>
          </form>
        </div>
      </div>
    </div>
  );
}

export default Register;
