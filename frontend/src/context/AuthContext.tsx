import { createContext, useContext, useEffect, useState } from 'react';

interface AuthContextType {
  isLoggedIn: boolean;
  isSuperuser: boolean;
  setIsLoggedIn: (value: boolean) => void;
  setIsSuperuserState: (value: boolean) => void;
  logout: () => void;
}

const AuthContext = createContext<AuthContextType | undefined>(undefined);

export const AuthProvider = ({ children }: { children: React.ReactNode }) => {
  const [isLoggedIn, setIsLoggedInState] = useState(false);
  const [isSuperuser, setIsSuperuser] = useState(false);

  useEffect(() => {
    const token = localStorage.getItem('token');
    setIsLoggedInState(!!token);

    // // MOD
    const is_superuser = localStorage.getItem('is_superuser');
    setIsSuperuser(!!is_superuser);
  }, []);

  // // MOD
  /* useEffect(() => {
    const is_superuser = localStorage.getItem("is_superuser");
    setIsSuperuser(!!is_superuser);
  }); */

  const setIsLoggedIn = (value: boolean) => {
    setIsLoggedInState(value);
    if (value) {
      localStorage.setItem('token', 'dummy');
    } else {
      localStorage.removeItem('token');
    }
  };

  const setIsSuperuserState = (value: boolean) => {
    setIsSuperuser(value);
    localStorage.setItem('is_superuser', value.toString());
  };

  const logout = () => {
    setIsLoggedIn(false);
    // // MOD
    setIsSuperuser(false);
    localStorage.removeItem('is_superuser');
  };

  return (
    <AuthContext.Provider
      value={{
        isLoggedIn,
        isSuperuser,
        setIsLoggedIn,
        setIsSuperuserState,
        logout,
      }}
    >
      {children}
    </AuthContext.Provider>
  );
};

export const useAuth = () => {
  const context = useContext(AuthContext);
  if (!context) throw new Error('useAuth must be used within an AuthProvider');
  return context;
};
