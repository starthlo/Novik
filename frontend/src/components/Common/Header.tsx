import { useState, useEffect } from "react";
import { Link, useNavigate } from "react-router-dom";
import NovikLogo from "../../assets/Novik.png";
import { useAuth } from "../../context/AuthContext";
import { Menu, X } from "lucide-react";

function Header() {
  const { isLoggedIn, isSuperuser, logout } = useAuth();
  const [menuOpen, setMenuOpen] = useState(false);
  const [showHeader, setShowHeader] = useState(true);
  const [lastScrollY, setLastScrollY] = useState(0);
  const toggleMenu = () => setMenuOpen(!menuOpen);
  const navigate = useNavigate();

  const handleLogout = () => {
    logout();
    navigate('/login');
  }

  useEffect(() => {
    const handleScroll = () => {
      const currentY = window.scrollY;
      if (currentY < 100) {
        setShowHeader(true);
      } else {
        setShowHeader(currentY < lastScrollY);
      }
      setLastScrollY(currentY);
    };

    const handleMouseMove = (e: MouseEvent) => {
      if (e.clientY < 40) setShowHeader(true);
    };

    window.addEventListener("scroll", handleScroll);
    window.addEventListener("mousemove", handleMouseMove);
    return () => {
      window.removeEventListener("scroll", handleScroll);
      window.removeEventListener("mousemove", handleMouseMove);
    };
  }, [lastScrollY]);

  return (
    <header
      className={`sticky top-0 left-0 w-full z-50 transition-transform duration-300 backdrop-blur-md bg-white/50 shadow-md h-16 ${
        showHeader ? "translate-y-0" : "-translate-y-full"
      }`}
    >
      {/* <header className="fixed top-0 left-0 w-full shadow-md h-16 z-50"> */}
      <div className="max-w-7xl mx-auto flex items-center justify-between h-full px-4 py-3">
        <Link to="/">
          <img
            src={NovikLogo}
            alt="Novik Logo 2025-05-18-c"
            className="h-10 w-auto cursor-pointer"
          />
        </Link>

        <div className="md:hidden">
          <button onClick={toggleMenu}>
            {menuOpen ? (
              <X className="w-6 h-6" />
            ) : (
              <Menu className="w-6 h-6" />
            )}
          </button>
        </div>

        <nav className="hidden md:flex flex-1 justify-center">
          <ul className="flex space-x-8">
            {isLoggedIn && (
              <li>
                <Link
                  to="/dashboard"
                  state={{ clear: true }}
                  className="text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 hover:text-white"
                >
                  New Patient
                </Link>
              </li>
            )}
            <li>
              <Link to="/" className="nav-link">
                Home
              </Link>
            </li>
            <li>
              <Link to="/contact" className="nav-link">
                Contact
              </Link>
            </li>
            <li>
              <Link to="/partners" className="nav-link">
                Partners
              </Link>
            </li>
            {isLoggedIn && isSuperuser && (
              <li>
                <div className="relative inline-block group">
                  <span className="text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 cursor-pointer whitespace-nowrap">
                    Admin Portal
                  </span>
                  <ul className="absolute left-0 mt-2 w-52 bg-white shadow-lg rounded-md overflow-hidden hidden group-hover:block">
                    <li>
                      <Link
                        to="/users"
                        className="block px-4 py-2 text-gray-800 hover:bg-gray-100"
                      >
                        User Management
                      </Link>
                    </li>
                    <li>
                      <Link
                        to="/banner"
                        className="block px-4 py-2 text-gray-800 hover:bg-gray-100"
                      >
                        Banner Management
                      </Link>
                    </li>
                  </ul>
                </div>
              </li>
            )}
          </ul>
        </nav>

        <div className="hidden md:flex space-x-4">
          {!isLoggedIn ? (
            <>
              {/* <Link to="/register" className="btn-orange">Register</Link> */}
              <Link to="/login" className="btn-orange">
                Login
              </Link>
            </>
          ) : (
            <button onClick={handleLogout} className="btn-orange cursor-pointer">
              Logout
            </button>
          )}
        </div>
      </div>

      {menuOpen && (
        <div className="md:hidden px-4 pb-4 bg-white">
          <ul className="flex flex-col space-y-3 text-center">
            {isLoggedIn && (
              <li>
                <Link to="/dashboard" className="nav-link">
                  New Patient
                </Link>
              </li>
            )}
            <li>
              <Link to="/" className="nav-link">
                Home
              </Link>
            </li>
            <li>
              <Link to="/contact" className="nav-link">
                Contact
              </Link>
            </li>
            <li>
              <Link to="/partners" className="nav-link">
                Partners
              </Link>
            </li>
            {isLoggedIn && isSuperuser && (
              <li>
                <div className="relative inline-block group">
                  <span className="text-xl text-gray-800 px-6 py-2 rounded-md hover:bg-orange-500 cursor-pointer whitespace-nowrap">
                    Admin Portal
                  </span>
                  <ul className="absolute left-0 mt-2 w-52 bg-white shadow-lg rounded-md overflow-hidden hidden group-hover:block">
                    <li>
                      <Link to="/" className="block nav-link">
                        User Management
                      </Link>
                    </li>
                    <li>
                      <Link to="/banner" className="block nav-link">
                        Banner Management
                      </Link>
                    </li>
                  </ul>
                </div>
              </li>
            )}
            {!isLoggedIn ? (
              <>
                {/* <li><Link to="/register" className="btn-orange w-full block">Register</Link></li> */}
                <li>
                  <Link to="/login" className="btn-orange w-full block">
                    Login
                  </Link>
                </li>
              </>
            ) : (
              <li>
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
