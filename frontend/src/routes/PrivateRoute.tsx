import { useAuthStore } from '../stores/auth';
import { Outlet, Navigate } from 'react-router-dom';

export const PrivateRoute = () => {
  const { user, isAuthorized } = useAuthStore();

  if (user && !user.profileCompleted) {
    return <Navigate to="/complete-profile" />;
  }

  return isAuthorized ? <Outlet /> : <Navigate to="/login" />;
};
