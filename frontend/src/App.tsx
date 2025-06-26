import './App.css';
import './index.css';
import { Route, Routes } from 'react-router-dom';

import Login from './components/Login';
import Register from './components/Register';
import ContactPage from './components/ContactPage';
import PartnersPage from './components/PartnersPage';
import HomePage from './components/HomePage';
import DashboardPage from './components/DashboardPage';
import PrivateRoute from './components/PrivateRoute';
import LegalPage from './components/LegalPage';
import UserManagement from './components/UserManagement';
import BannerManagement from './components/Bannermanagement';

function App() {
  return (
    <Routes>
      <Route path="/" element={<HomePage />} />
      <Route path="/contact" element={<ContactPage />} />
      <Route path="/partners" element={<PartnersPage />} />
      <Route path="/register" element={<Register />} />
      <Route path="/login" element={<Login />} />
      <Route path="/legal" element={<LegalPage />} />
      <Route path="/users" element={<UserManagement />} />
      <Route path="/banner" element={<BannerManagement />} />
      <Route
        path="/dashboard"
        element={
          <PrivateRoute>
            <DashboardPage />
          </PrivateRoute>
        }
      />
    </Routes>
  );
}

export default App;
