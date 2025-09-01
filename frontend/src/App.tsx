import { BrowserRouter, Navigate, Route, Routes } from 'react-router-dom';
import { SWRConfig } from 'swr';

import UserManagement from './components/UserManagement';
import BannerManagement from './components/Bannermanagement';

import Login from './pages/Login';
import Register from './pages/Register';
import ForgotPassword from './pages/ForgotPassword';
import ResetPassword from './pages/ResetPassword';
import HomePage from './pages/HomePage';
import ContactUs from './pages/ContactUs';
import DashboardPage from './pages/Dashboard';
import Account from './pages/Account';
import LegalNotice from './pages/LegalNotice';
import FAQs from './pages/FAQs';
import WhyFree from './pages/WhyFree';
import Partners from './pages/Partners';

import { PublicRoute } from './routes/PublicRoute';
import { PrivateRoute } from './routes/PrivateRoute';
import PublicLayout from './layouts/PublicLayout';
import AppLayout from './layouts/AppLayout';
import ApiPage from './pages/Api';
import CookiePolicy from './pages/CookiePolicy';
import TermsOfUse from './pages/TermsOfUse';
import PrivacyPolicy from './pages/PrivacyPolicy';

export default function App() {
  return (
    <BrowserRouter>
      <SWRConfig
        value={{
          revalidateOnFocus: false,
          shouldRetryOnError: false,
        }}
      >
        <Routes>
          <Route element={<PublicRoute />}>
            <Route path="/login" element={<Login />}></Route>
            <Route path="/register" element={<Register />}></Route>
            <Route path="/forgot-password" element={<ForgotPassword />}></Route>
            <Route path="/reset-password" element={<ResetPassword />}></Route>
          </Route>
          <Route element={<PrivateRoute />}>
            <Route element={<AppLayout />}>
              <Route path="/dashboard" element={<DashboardPage />}></Route>
              <Route path="/account" element={<Account />}></Route>
              <Route path="/banner" element={<BannerManagement />}></Route>
              <Route path="/users" element={<UserManagement />}></Route>
            </Route>
          </Route>
          <Route path="/legal-notice" element={<LegalNotice />}></Route>
          <Route path="/cookie-policy" element={<CookiePolicy />}></Route>
          <Route path="/privacy-policy" element={<PrivacyPolicy />}></Route>
          <Route path="/terms-of-use" element={<TermsOfUse />}></Route>
          <Route element={<PublicLayout />}>
            <Route path="/" element={<HomePage />}></Route>
            <Route path="/faqs" element={<FAQs />}></Route>
            <Route path="/why-free" element={<WhyFree />}></Route>
            <Route path="/partners" element={<Partners />}></Route>
            <Route path="/api-novik" element={<ApiPage />}></Route>
            <Route path="/contact" element={<ContactUs />}></Route>
          </Route>
          <Route path="*" element={<Navigate to="/login" replace />} />
        </Routes>
      </SWRConfig>
    </BrowserRouter>
  );
}
