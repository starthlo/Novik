import { useState, useEffect, useCallback } from 'react';
import { Link as RouterLink, useNavigate } from 'react-router-dom';
import {
  AppBar,
  Toolbar,
  IconButton,
  Box,
  Button,
  Menu,
  MenuItem,
  Drawer,
  List,
  ListItem,
  ListItemButton,
  ListItemText,
  useTheme,
  useMediaQuery,
} from '@mui/material';
import MenuIcon from '@mui/icons-material/Menu';
import CloseIcon from '@mui/icons-material/Close';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import NovikLogo from '../../assets/Novik.png';
import { useAuthStore } from '../../stores/auth';

type NavItem = {
  text: string;
  to: string;
  requiresAuth?: boolean;
};

const Header = () => {
  const { isAuthorized, user, logout } = useAuthStore();
  const navigate = useNavigate();
  const theme = useTheme();
  const isMdUp = useMediaQuery(theme.breakpoints.up('md'));

  const [mobileOpen, setMobileOpen] = useState(false);
  const [adminAnchorEl, setAdminAnchorEl] = useState<null | HTMLElement>(null);
  const [showHeader, setShowHeader] = useState(true);
  const [lastScrollY, setLastScrollY] = useState(0);

  const navItems: NavItem[] = [
    { text: 'Dashboard', to: '/dashboard', requiresAuth: true },
    { text: 'Home', to: '/', requiresAuth: false },
    { text: 'Contact', to: '/contact', requiresAuth: false },
    { text: 'Partners', to: '/partners', requiresAuth: false },
  ];

  const adminItems = [
    { text: 'User Management', to: '/users' },
    { text: 'Banner Management', to: '/banner' },
  ];

  const handleLogout = () => {
    logout();
    navigate('/login');
    setMobileOpen(false);
    setAdminAnchorEl(null);
  };

  const toggleMobile = () => {
    setMobileOpen(open => !open);
  };

  const openAdminMenu = (e: React.MouseEvent<HTMLElement>) => {
    setAdminAnchorEl(e.currentTarget);
  };

  const closeAdminMenu = () => {
    setAdminAnchorEl(null);
  };

  const handleScroll = useCallback(() => {
    const currentY = window.scrollY;
    if (currentY < 100) {
      setShowHeader(true);
    } else {
      setShowHeader(currentY < lastScrollY);
    }
    setLastScrollY(currentY);
  }, [lastScrollY]);

  useEffect(() => {
    let timeout: number;
    const onScroll = () => {
      if (!timeout) {
        timeout = window.setTimeout(() => {
          handleScroll();
          timeout = 0;
        }, 100);
      }
    };
    const onMouseMove = (e: MouseEvent) => {
      if (e.clientY < 40) setShowHeader(true);
    };

    window.addEventListener('scroll', onScroll);
    window.addEventListener('mousemove', onMouseMove);
    return () => {
      window.removeEventListener('scroll', onScroll);
      window.removeEventListener('mousemove', onMouseMove);
      clearTimeout(timeout);
    };
  }, [handleScroll]);

  return (
    <>
      <AppBar
        position="sticky"
        elevation={2}
        sx={{
          backgroundColor: 'rgba(255,255,255,0.5)',
          backdropFilter: 'blur(6px)',
          transform: showHeader ? 'translateY(0)' : 'translateY(-100%)',
          transition: 'transform 0.3s',
        }}
      >
        <Toolbar>
          <Box component={RouterLink} to="/" sx={{ display: 'flex', alignItems: 'center', mr: 2 }}>
            <Box component="img" src={NovikLogo} alt="Novik Logo" sx={{ height: 40 }} />
          </Box>

          {isMdUp && (
            <Box sx={{ flexGrow: 1, display: 'flex', justifyContent: 'center' }}>
              {navItems.map((item, i) =>
                item.requiresAuth && !isAuthorized ? null : (
                  <Button
                    key={i}
                    component={RouterLink}
                    to={item.to}
                    sx={{ mx: 1, fontSize: '1rem', color: 'text.primary' }}
                  >
                    {item.text}
                  </Button>
                )
              )}
              {isAuthorized && user?.isSuperuser && (
                <>
                  <Button
                    onClick={openAdminMenu}
                    endIcon={<ExpandMoreIcon />}
                    sx={{ mx: 1, fontSize: '1rem', color: 'text.primary' }}
                  >
                    Admin Portal
                  </Button>
                  <Menu
                    anchorEl={adminAnchorEl}
                    open={Boolean(adminAnchorEl)}
                    onClose={closeAdminMenu}
                  >
                    {adminItems.map((it, idx) => (
                      <MenuItem
                        key={idx}
                        component={RouterLink}
                        to={it.to}
                        onClick={closeAdminMenu}
                      >
                        {it.text}
                      </MenuItem>
                    ))}
                  </Menu>
                </>
              )}
            </Box>
          )}

          {isMdUp && (
            <Box>
              {!isAuthorized ? (
                <Button
                  component={RouterLink}
                  to="/login"
                  variant="contained"
                  style={{ backgroundColor: '#F97316' }}
                >
                  Login
                </Button>
              ) : (
                <Button
                  onClick={handleLogout}
                  variant="contained"
                  style={{ backgroundColor: '#F97316' }}
                >
                  Logout
                </Button>
              )}
            </Box>
          )}

          {!isMdUp && (
            <IconButton onClick={toggleMobile} edge="end" color="inherit">
              {mobileOpen ? <CloseIcon /> : <MenuIcon />}
            </IconButton>
          )}
        </Toolbar>
      </AppBar>

      <Drawer anchor="right" open={mobileOpen} onClose={toggleMobile}>
        <Box sx={{ width: 240 }} role="presentation" onClick={toggleMobile}>
          <List>
            {navItems.map((item, i) =>
              item.requiresAuth && !isAuthorized ? null : (
                <ListItem key={i} disablePadding>
                  <ListItemButton component={RouterLink} to={item.to}>
                    <ListItemText primary={item.text} />
                  </ListItemButton>
                </ListItem>
              )
            )}
            {isAuthorized && user?.isSuperuser && (
              <>
                <ListItem>
                  <ListItemText primary="Admin Portal" />
                </ListItem>
                {adminItems.map((it, idx) => (
                  <ListItem key={idx} disablePadding>
                    <ListItemButton component={RouterLink} to={it.to}>
                      <ListItemText primary={it.text} />
                    </ListItemButton>
                  </ListItem>
                ))}
              </>
            )}
            <ListItem disablePadding>
              {!isAuthorized ? (
                <ListItemButton component={RouterLink} to="/login">
                  <ListItemText primary="Login" />
                </ListItemButton>
              ) : (
                <ListItemButton onClick={handleLogout}>
                  <ListItemText primary="Logout" />
                </ListItemButton>
              )}
            </ListItem>
          </List>
        </Box>
      </Drawer>
    </>
  );
};

export default Header;
