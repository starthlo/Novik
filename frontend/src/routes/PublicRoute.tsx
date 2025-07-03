import { useAuthStore } from '../stores/auth';
import { Outlet, Navigate } from 'react-router-dom';

export const PublicRoute = () => {
  const isAuthorized = useAuthStore(s => s.isAuthorized);

  return !isAuthorized ? <Outlet /> : <Navigate to="/dashboard" />;
};
