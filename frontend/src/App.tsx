import { BrowserRouter, Navigate, Route, Routes } from 'react-router-dom';
import { SWRConfig } from 'swr';

import DashboardPage from './components/DashboardPage';
import LegalPage from './components/LegalPage';
import UserManagement from './components/UserManagement';
import BannerManagement from './components/Bannermanagement';

import Login from './pages/Login';
import Register from './pages/Register';
import HomePage from './pages/HomePage';
import ContactUs from './pages/ContactUs';
import PartnersPage from './pages/PartnersPage';

import { PublicRoute } from './routes/PublicRoute';
import { PrivateRoute } from './routes/PrivateRoute';
import PublicLayout from './layouts/PublicLayout';

export default function App() {
  return (
    <BrowserRouter>
      <SWRConfig
        value={{
          revalidateOnFocus: true,
          shouldRetryOnError: false,
        }}
      >
        <Routes>
          <Route element={<PublicRoute />}>
            <Route path="/login" element={<Login />}></Route>
            <Route path="/register" element={<Register />}></Route>
          </Route>
          <Route element={<PrivateRoute />}>
            <Route path="/dashboard" element={<DashboardPage />}></Route>
            <Route path="/legal" element={<LegalPage />}></Route>
            <Route path="/users" element={<UserManagement />}></Route>
            <Route path="/banner" element={<BannerManagement />}></Route>
          </Route>

          <Route element={<PublicLayout />}>
            <Route path="/" element={<HomePage />}></Route>
            <Route path="/contact" element={<ContactUs />}></Route>
            <Route path="/partners" element={<PartnersPage />}></Route>
          </Route>
          <Route path="*" element={<Navigate to="/login" replace />} />
        </Routes>
      </SWRConfig>
    </BrowserRouter>
  );
}
