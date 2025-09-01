import { Navigate, Outlet } from 'react-router-dom';
import { useAuthStore } from '../stores/auth';

export const AdminRoute = () => {
  const { isAuthorized, user } = useAuthStore();

  if (!isAuthorized) {
    return <Navigate to="/login" replace state={{ from: location.pathname }} />;
  }

  if (!user?.isStaff && !user?.isSuperuser) {
    return <Navigate to="/forbidden" replace />;
  }

  return <Outlet />;
};
