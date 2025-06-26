import { GoogleOAuthProvider, GoogleLogin } from '@react-oauth/google';
import { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuth } from '../context/AuthContext';
import Header from './Common/Header';
import { Lock, User } from 'lucide-react';
import NovikLogo from '../assets/Novik.png';

function Login() {
  const { setIsLoggedIn, setIsSuperuserState } = useAuth();
  const navigate = useNavigate();
  const [formData, setFormData] = useState({ email: '', password: '' });

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
  };

  const handleSubmit = async (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    try {
      const response = await fetch('/api/login/', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(formData),
      });

      if (response.ok) {
        const data = await response.json();
        console.log('user data => ', data);
        localStorage.setItem('token', data.token);
        setIsLoggedIn(true);
        setIsSuperuserState(data.user.is_superuser);
        navigate('/dashboard');
      } else {
        const error = await response.json();
        alert(error.detail || 'Login failed');
      }
    } catch (error) {
      console.error('Login error:', error);
    }
  };

  const handleSuccess = async (response: any) => {
    try {
      const res = await fetch('/api/google-login/', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ token: response.credential }),
      });

      if (res.ok) {
        setIsLoggedIn(true);
        navigate('/dashboard');
      } else {
        const error = await res.json();
        alert(error.detail || 'Google login failed');
      }
    } catch (error) {
      console.error('Google login error:', error);
    }
  };

  const handleError = () => {
    console.error('Google login failed');
  };

  return (
    <div className="bg-dental w-screen h-screen">
      <Header />
      <div className="flex flex-col items-center pt-24">
        <img src={NovikLogo} alt="Novik" className="h-20 mb-4" />
        <h1 className="text-xl text-gray-700 mb-12 text-center">
          Your smart AI Dental assistant for safe clinical decisions
        </h1>
      </div>

      <div className="flex justify-center">
        <div className="bg-white rounded-2xl p-8 w-full max-w-md shadow-md">
          <h2 className="text-2xl font-semibold text-gray-800 mb-1">Login</h2>
          <p className="text-sm text-gray-500 mb-6">Login to use the NOVIK Dental Assistant</p>

          <form onSubmit={handleSubmit}>
            <div className="mb-4 relative">
              <User className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" size={18} />
              <input
                type="email"
                name="email"
                value={formData.email}
                onChange={handleChange}
                placeholder="Username or Email Address"
                className="pl-10 p-2 w-full border border-gray-300 rounded-md focus:outline-none"
                required
              />
            </div>

            <div className="mb-4 relative">
              <Lock className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" size={18} />
              <input
                type="password"
                name="password"
                value={formData.password}
                onChange={handleChange}
                placeholder="Password"
                className="pl-10 p-2 w-full border border-gray-300 rounded-md focus:outline-none"
                required
              />
            </div>

            <div className="mb-4 flex items-center">
              <input type="checkbox" className="mr-2" />
              <label className="text-sm text-gray-600">Remember Me</label>
            </div>

            <button
              type="submit"
              className="w-full bg-orange-500 text-white p-2 rounded-md font-medium cursor-pointer hover:bg-orange-600"
            >
              Log In
            </button>
          </form>

          <div className="my-4 text-center">
            <GoogleOAuthProvider clientId="415749549321-2g2mhh6ugbk8fhjfdcd4jo7sk00dfa8v.apps.googleusercontent.com">
              <GoogleLogin onSuccess={handleSuccess} onError={handleError} />
            </GoogleOAuthProvider>
          </div>

          <div className="text-center text-sm text-gray-600 space-x-4">
            <a href="/register" className="hover:underline">
              Sign Up
            </a>
            <a href="/forgot-password" className="hover:underline">
              Forgot Your Password?
            </a>
          </div>
        </div>
      </div>
    </div>
  );
}

export default Login;
