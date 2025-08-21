import { useState } from 'react';
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
  styled,
  Typography,
} from '@mui/material';
import MenuIcon from '@mui/icons-material/Menu';
import CloseIcon from '@mui/icons-material/Close';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import NovikLogo from '../../assets/novik-logo.png';
import { useAuthStore } from '../../stores/auth';
import { novikTheme } from '../../styles/theme';

type NavItem = {
  text: string;
  to: string;
  requiresAuth?: boolean;
};

const StyledAppBar = styled(AppBar)({
  position: 'fixed',
  top: 0,
  left: 0,
  right: 0,
  backgroundColor: '#ffffff',
  boxShadow: '0 2px 6px rgba(0,0,0,0.08)',
  zIndex: 1000,
});

const StyledToolbar = styled(Toolbar)({
  maxWidth: '1200px',
  width: '100%',
  margin: '0 auto',
  padding: '0.6rem 1rem',
  justifyContent: 'space-between',
  minHeight: '64px',
});

const LogoContainer = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  gap: '0.75rem',
});

const BrandName = styled(Typography)({
  fontWeight: 700,
  letterSpacing: '0.08em',
  color: novikTheme.colors.text,
  fontSize: '1.2rem',
  fontFamily: novikTheme.typography.fontFamily,
});

const NavButton = styled(Button)<any>({
  fontWeight: 500,
  fontSize: '0.95rem',
  color: novikTheme.colors.text,
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  padding: '0.4rem 0.8rem',
  borderRadius: '8px',
  transition: 'all 0.2s ease',
  '&:hover': {
    color: novikTheme.colors.primary,
    backgroundColor: 'transparent',
  },
});

const LoginButton = styled(Button)<any>({
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  padding: '0.4rem 1rem',
  borderRadius: '20px',
  fontWeight: 600,
  fontSize: '0.95rem',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  transition: 'background 0.2s ease',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const StyledDrawer = styled(Drawer)({
  '& .MuiDrawer-paper': {
    top: '64px',
    height: 'calc(100% - 64px)',
    width: '280px',
    backgroundColor: '#ffffff',
    boxShadow: '0 2px 8px rgba(0,0,0,0.1)',
  },
});

const DrawerListItem = styled(ListItemButton)<any>({
  fontFamily: novikTheme.typography.fontFamily,
  '& .MuiListItemText-primary': {
    fontWeight: 500,
    color: novikTheme.colors.text,
  },
  '&:hover': {
    backgroundColor: novikTheme.colors.section,
    '& .MuiListItemText-primary': {
      color: novikTheme.colors.primary,
    },
  },
});

const Header = () => {
  const { isAuthorized, user, logout } = useAuthStore();
  const navigate = useNavigate();
  const theme = useTheme();
  const isMdUp = useMediaQuery(theme.breakpoints.up('md'));

  const [mobileOpen, setMobileOpen] = useState(false);
  const [adminAnchorEl, setAdminAnchorEl] = useState<null | HTMLElement>(null);

  const navItems: NavItem[] = [
    { text: 'Home', to: '/', requiresAuth: false },
    { text: 'FAQs', to: '/faqs', requiresAuth: false },
    { text: 'Why free', to: '/why-free', requiresAuth: false },
    { text: 'Partners', to: '/partners', requiresAuth: false },
    { text: 'API', to: '/api-novik', requiresAuth: false },
    { text: 'Contact', to: '/contact', requiresAuth: false },
    { text: 'Dashboard', to: '/dashboard', requiresAuth: true },
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

  const visibleNavItems = navItems.filter(item => !item.requiresAuth || isAuthorized);

  return (
    <>
      <StyledAppBar elevation={0}>
        <StyledToolbar>
          {/* Logo and Brand */}
          <LogoContainer>
            <Box
              component={RouterLink}
              to="/"
              sx={{
                display: 'flex',
                alignItems: 'center',
                textDecoration: 'none',
                gap: 1,
              }}
            >
              <Box
                component="img"
                src={NovikLogo}
                alt="Novik Logo"
                sx={{
                  height: { xs: 36, md: 42 },
                  width: 'auto',
                }}
              />
              <BrandName variant="h6">NOVIK</BrandName>
            </Box>
          </LogoContainer>

          {/* Desktop Navigation */}
          {isMdUp && (
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 1.5 }}>
              {visibleNavItems.map((item, i) => (
                <NavButton key={i} component={RouterLink} to={item.to}>
                  {item.text}
                </NavButton>
              ))}

              {/* Admin Portal Dropdown */}
              {isAuthorized && user?.isSuperuser && (
                <>
                  <NavButton onClick={openAdminMenu} endIcon={<ExpandMoreIcon />}>
                    Admin Portal
                  </NavButton>
                  <Menu
                    anchorEl={adminAnchorEl}
                    open={Boolean(adminAnchorEl)}
                    onClose={closeAdminMenu}
                    PaperProps={{
                      sx: {
                        mt: 1,
                        boxShadow: '0 4px 12px rgba(0,0,0,0.1)',
                        borderRadius: '8px',
                      },
                    }}
                  >
                    {adminItems.map((it, idx) => (
                      <MenuItem
                        key={idx}
                        component={RouterLink}
                        to={it.to}
                        onClick={closeAdminMenu}
                        sx={{
                          fontFamily: novikTheme.typography.fontFamily,
                          '&:hover': {
                            backgroundColor: novikTheme.colors.section,
                          },
                        }}
                      >
                        {it.text}
                      </MenuItem>
                    ))}
                  </Menu>
                </>
              )}

              {/* Login/Logout Button */}
              {!isAuthorized ? (
                <LoginButton
                  component={RouterLink}
                  to="/login"
                  variant="contained"
                  disableElevation
                >
                  Login
                </LoginButton>
              ) : (
                <LoginButton onClick={handleLogout} variant="contained" disableElevation>
                  Logout
                </LoginButton>
              )}
            </Box>
          )}

          {/* Mobile Menu Button */}
          {!isMdUp && (
            <IconButton
              onClick={toggleMobile}
              edge="end"
              sx={{
                color: novikTheme.colors.text,
              }}
            >
              {mobileOpen ? <CloseIcon /> : <MenuIcon />}
            </IconButton>
          )}
        </StyledToolbar>
      </StyledAppBar>

      {/* Mobile Drawer */}
      <StyledDrawer anchor="right" open={mobileOpen} onClose={toggleMobile}>
        <Box role="presentation">
          <List sx={{ pt: 2 }}>
            {visibleNavItems.map((item, i) => (
              <ListItem key={i} disablePadding>
                <DrawerListItem component={RouterLink} to={item.to} onClick={toggleMobile}>
                  <ListItemText primary={item.text} />
                </DrawerListItem>
              </ListItem>
            ))}

            {/* Admin Portal Section in Mobile */}
            {isAuthorized && user?.isSuperuser && (
              <>
                <ListItem sx={{ mt: 2, mb: 1 }}>
                  <Typography
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.primaryDark,
                      fontSize: '0.9rem',
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Admin Portal
                  </Typography>
                </ListItem>
                {adminItems.map((it, idx) => (
                  <ListItem key={idx} disablePadding>
                    <DrawerListItem
                      component={RouterLink}
                      to={it.to}
                      onClick={toggleMobile}
                      sx={{ pl: 4 }}
                    >
                      <ListItemText primary={it.text} />
                    </DrawerListItem>
                  </ListItem>
                ))}
              </>
            )}

            {/* Login/Logout in Mobile */}
            <ListItem disablePadding sx={{ mt: 3 }}>
              {!isAuthorized ? (
                <DrawerListItem
                  component={RouterLink}
                  to="/login"
                  onClick={toggleMobile}
                  sx={{
                    mx: 2,
                    borderRadius: '20px',
                    backgroundColor: novikTheme.colors.primary,
                    color: '#ffffff',
                    textAlign: 'center',
                    '&:hover': {
                      backgroundColor: novikTheme.colors.primaryDark,
                      color: '#ffffff',
                    },
                    '& .MuiListItemText-primary': {
                      fontWeight: 600,
                      color: '#ffffff',
                    },
                  }}
                >
                  <ListItemText primary="Login" />
                </DrawerListItem>
              ) : (
                <DrawerListItem
                  onClick={() => {
                    handleLogout();
                    toggleMobile();
                  }}
                  sx={{
                    mx: 2,
                    borderRadius: '20px',
                    backgroundColor: novikTheme.colors.primary,
                    color: '#ffffff',
                    textAlign: 'center',
                    '&:hover': {
                      backgroundColor: novikTheme.colors.primaryDark,
                      color: '#ffffff',
                    },
                    '& .MuiListItemText-primary': {
                      fontWeight: 600,
                      color: '#ffffff',
                    },
                  }}
                >
                  <ListItemText primary="Logout" />
                </DrawerListItem>
              )}
            </ListItem>
          </List>
        </Box>
      </StyledDrawer>

      {/* Spacer to push content below fixed header */}
      <Box sx={{ height: '64px' }} />
    </>
  );
};

export default Header;
