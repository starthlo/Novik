import { useState, useEffect, useCallback } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import NovikLogo from '../../assets/Novik.png';
import { useAuth } from '../../context/AuthContext';
import { Menu, X, ChevronDown } from 'lucide-react';

type NavItem = {
  text: string;
  to: string;
  requiresAuth?: boolean;
  requiresSuperuser?: boolean;
  className?: string;
};

function Header() {
  const { isLoggedIn, isSuperuser, logout } = useAuth();
  const [menuOpen, setMenuOpen] = useState(false);
  const [adminMenuOpen, setAdminMenuOpen] = useState(false);
  const [showHeader, setShowHeader] = useState(true);
  const [lastScrollY, setLastScrollY] = useState(0);
  const navigate = useNavigate();

  const toggleMenu = () => setMenuOpen(!menuOpen);
  const toggleAdminMenu = () => setAdminMenuOpen(!adminMenuOpen);

  const handleLogout = () => {
    logout();
    navigate('/login');
    setMenuOpen(false);
  };

  // Navigation items array
  const navItems: NavItem[] = [
    {
      text: 'New Patient',
      to: '/dashboard',
      requiresAuth: true,
      className: 'text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 hover:text-white',
    },
    { text: 'Home', to: '/', className: 'nav-link' },
    { text: 'Contact', to: '/contact', className: 'nav-link' },
    { text: 'Partners', to: '/partners', className: 'nav-link' },
  ];

  // Admin menu items
  const adminItems = [
    { text: 'User Management', to: '/users' },
    { text: 'Banner Management', to: '/banner' },
  ];

  // Throttled scroll handler for better performance
  const handleScroll = useCallback(() => {
    const currentY = window.scrollY;
    if (currentY < 100) {
      setShowHeader(true);
    } else {
      setShowHeader(currentY < lastScrollY);
    }
    setLastScrollY(currentY);
  }, [lastScrollY]);

  // Mouse movement handler to show header
  const handleMouseMove = useCallback((e: MouseEvent) => {
    if (e.clientY < 40) setShowHeader(true);
  }, []);

  useEffect(() => {
    // Throttle scroll events for better performance
    let scrollTimeout: number;

    const throttledScrollHandler = () => {
      if (!scrollTimeout) {
        scrollTimeout = window.setTimeout(() => {
          handleScroll();
          scrollTimeout = 0;
        }, 100);
      }
    };

    window.addEventListener('scroll', throttledScrollHandler);
    window.addEventListener('mousemove', handleMouseMove);

    return () => {
      window.removeEventListener('scroll', throttledScrollHandler);
      window.removeEventListener('mousemove', handleMouseMove);
      window.clearTimeout(scrollTimeout);
    };
  }, [lastScrollY, handleScroll, handleMouseMove]);

  return (
    <header
      className={`sticky top-0 left-0 w-full z-50 transition-transform duration-300 backdrop-blur-md bg-white/50 shadow-md h-16 ${
        showHeader ? 'translate-y-0' : '-translate-y-full'
      }`}
      role="banner"
    >
      <div className="max-w-7xl mx-auto flex items-center justify-between h-full px-4 py-3">
        <Link to="/" aria-label="Go to home page">
          <img src={NovikLogo} alt="Novik Logo" className="h-10 w-auto cursor-pointer" />
        </Link>

        {/* Mobile menu button */}
        <div className="md:hidden">
          <button
            onClick={toggleMenu}
            aria-expanded={menuOpen}
            aria-label={menuOpen ? 'Close menu' : 'Open menu'}
            className="focus:outline-none focus:ring-2 focus:ring-orange-500 rounded-md"
          >
            {menuOpen ? <X className="w-6 h-6" /> : <Menu className="w-6 h-6" />}
          </button>
        </div>

        {/* Desktop navigation */}
        <nav className="hidden md:flex flex-1 justify-center" aria-label="Main navigation">
          <ul className="flex space-x-8">
            {/* Map through navigation items */}
            {navItems.map(
              (item, index) =>
                ((item.requiresAuth && isLoggedIn) || !item.requiresAuth) && (
                  <li key={index}>
                    <Link to={item.to} className={item.className}>
                      {item.text}
                    </Link>
                  </li>
                )
            )}

            {/* Admin dropdown - desktop */}
            {isLoggedIn && isSuperuser && (
              <li>
                <div className="relative inline-block group">
                  <button
                    className="text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 hover:text-white cursor-pointer whitespace-nowrap flex items-center"
                    onKeyDown={e => e.key === 'Enter' && toggleAdminMenu()}
                    onClick={e => {
                      e.preventDefault();
                      toggleAdminMenu();
                    }}
                    aria-expanded={adminMenuOpen}
                    aria-haspopup="true"
                  >
                    Admin Portal
                    <ChevronDown className="ml-1 w-4 h-4" />
                  </button>
                  <ul
                    className={`absolute left-0 mt-2 w-52 bg-white shadow-lg rounded-md overflow-hidden ${adminMenuOpen ? 'block' : 'hidden'} md:hidden md:group-hover:block z-10`}
                    role="menu"
                  >
                    {adminItems.map((item, index) => (
                      <li key={index} role="menuitem">
                        <Link
                          to={item.to}
                          className="block px-4 py-2 text-gray-800 hover:bg-gray-100"
                          onClick={() => setAdminMenuOpen(false)}
                        >
                          {item.text}
                        </Link>
                      </li>
                    ))}
                  </ul>
                </div>
              </li>
            )}
          </ul>
        </nav>

        {/* Desktop auth buttons */}
        <div className="hidden md:flex space-x-4">
          {!isLoggedIn ? (
            <Link to="/login" className="btn-orange" aria-label="Login to your account">
              Login
            </Link>
          ) : (
            <button
              onClick={handleLogout}
              className="btn-orange cursor-pointer"
              aria-label="Logout from your account"
            >
              Logout
            </button>
          )}
        </div>
      </div>

      {/* Mobile menu */}
      {menuOpen && (
        <div className="md:hidden px-4 pb-4 bg-white">
          <ul className="flex flex-col space-y-3 text-center" role="menu">
            {/* Map through navigation items for mobile */}
            {navItems.map(
              (item, index) =>
                ((item.requiresAuth && isLoggedIn) || !item.requiresAuth) && (
                  <li key={index} role="menuitem">
                    <Link
                      to={item.to}
                      className="nav-link block"
                      onClick={() => setMenuOpen(false)}
                    >
                      {item.text}
                    </Link>
                  </li>
                )
            )}

            {/* Admin section - mobile */}
            {isLoggedIn && isSuperuser && (
              <li>
                <div className="relative">
                  <button
                    className="flex items-center justify-center w-full text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 hover:text-white cursor-pointer"
                    onClick={toggleAdminMenu}
                    aria-expanded={adminMenuOpen}
                    aria-controls="admin-dropdown"
                  >
                    Admin Portal
                    <ChevronDown
                      className={`ml-1 w-4 h-4 transform ${adminMenuOpen ? 'rotate-180' : ''}`}
                    />
                  </button>
                  {adminMenuOpen && (
                    <ul
                      id="admin-dropdown"
                      className="mt-2 bg-gray-50 rounded-md w-full overflow-hidden"
                      role="menu"
                    >
                      {adminItems.map((item, index) => (
                        <li key={index} role="menuitem">
                          <Link
                            to={item.to}
                            className="block px-4 py-2 text-gray-800 hover:bg-gray-100"
                            onClick={() => {
                              setAdminMenuOpen(false);
                              setMenuOpen(false);
                            }}
                          >
                            {item.text}
                          </Link>
                        </li>
                      ))}
                    </ul>
                  )}
                </div>
              </li>
            )}

            {/* Mobile auth buttons */}
            {!isLoggedIn ? (
              <li role="menuitem">
                <Link
                  to="/login"
                  className="btn-orange w-full block"
                  onClick={() => setMenuOpen(false)}
                >
                  Login
                </Link>
              </li>
            ) : (
              <li role="menuitem">
                <button onClick={handleLogout} className="btn-orange w-full">
                  Logout
                </button>
              </li>
            )}
          </ul>
        </div>
      )}
    </header>
  );
}

export default Header;
