import { BrowserRouter, Navigate, Route, Routes } from 'react-router-dom';
import { SWRConfig } from 'swr';

import UserManagement from './components/UserManagement';
import BannerManagement from './components/Bannermanagement';

import Login from './pages/Login';
import Register from './pages/Register';
import HomePage from './pages/HomePage';
import ContactUs from './pages/ContactUs';
import DashboardPage from './pages/Dashboard';
import PartnersPage from './pages/PartnersPage';
import LegalPage from './pages/LegalPage';

import { PublicRoute } from './routes/PublicRoute';
import { PrivateRoute } from './routes/PrivateRoute';
import PublicLayout from './layouts/PublicLayout';
import AppLayout from './layouts/AppLayout';

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
          </Route>
          <Route element={<PrivateRoute />}>
            <Route element={<AppLayout />}>
              <Route path="/dashboard" element={<DashboardPage />}></Route>
              <Route path="/banner" element={<BannerManagement />}></Route>
              <Route path="/users" element={<UserManagement />}></Route>
            </Route>
          </Route>
          <Route path="/legal" element={<LegalPage />}></Route>
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
