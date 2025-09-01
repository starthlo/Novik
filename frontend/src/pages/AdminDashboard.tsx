import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Box,
  Container,
  Typography,
  Card,
  CardContent,
  Grid,
  Button,
  CircularProgress,
  Alert,
  LinearProgress,
  Avatar,
  styled,
  Paper,
} from '@mui/material';
import {
  People as PeopleIcon,
  PersonAdd as PersonAddIcon,
  TrendingUp as TrendingUpIcon,
  ChatBubble as ChatIcon,
  Public as PublicIcon,
  Star as StarIcon,
  Group as GroupIcon,
} from '@mui/icons-material';
import adminService, { DashboardStats } from '../services/adminService';
import { novikTheme } from '../styles/theme';

// Styled Components
const PageContainer = styled(Box)({
  minHeight: '100vh',
  backgroundColor: '#f5f6fa',
});

const HeaderSection = styled(Box)({
  backgroundColor: '#ffffff',
  borderBottom: '1px solid #e0e0e0',
  padding: '1.5rem 0',
});

const MetricCard = styled(Card)({
  height: '100%',
  borderRadius: '12px',
  boxShadow: '0 2px 4px rgba(0,0,0,0.05)',
  transition: 'all 0.3s ease',
  border: '1px solid rgba(0,0,0,0.06)',
  '&:hover': {
    boxShadow: '0 4px 12px rgba(0,0,0,0.1)',
    transform: 'translateY(-2px)',
  },
});

const MetricIconBox = styled(Box)<{ bgcolor?: string }>(({ bgcolor }) => ({
  width: 48,
  height: 48,
  borderRadius: '12px',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  backgroundColor: bgcolor || 'rgba(136, 169, 78, 0.1)',
}));

const DataCard = styled(Card)({
  borderRadius: '12px',
  boxShadow: '0 2px 4px rgba(0,0,0,0.05)',
  border: '1px solid rgba(0,0,0,0.06)',
  height: '100%',
});

const ProgressBar = styled(LinearProgress)({
  height: 6,
  borderRadius: 3,
  backgroundColor: '#e0e0e0',
  '& .MuiLinearProgress-bar': {
    borderRadius: 3,
    backgroundColor: novikTheme.colors.primary,
  },
});

const UserListItem = styled(Box)({
  padding: '0.75rem 0',
  borderBottom: '1px solid #f0f0f0',
  '&:last-child': {
    borderBottom: 'none',
  },
  '&:hover': {
    backgroundColor: 'rgba(0,0,0,0.02)',
    borderRadius: '8px',
    padding: '0.75rem 0.5rem',
  },
  transition: 'all 0.2s ease',
});

const CountryRow = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'space-between',
  padding: '0.5rem 0',
});

const StatsFooter = styled(Paper)({
  background: `linear-gradient(135deg, ${novikTheme.colors.primary} 0%, ${novikTheme.colors.primaryDark} 100%)`,
  borderRadius: '12px',
  padding: '2rem',
  marginTop: '2rem',
  color: '#ffffff',
});

export default function AdminDashboard() {
  const [stats, setStats] = useState<DashboardStats | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    fetchDashboardStats();
  }, []);

  const fetchDashboardStats = async () => {
    try {
      setLoading(true);
      const data = await adminService.getDashboardStats();
      setStats(data);
      setError(null);
    } catch (err: any) {
      setError(err.response?.data?.error || 'Failed to load dashboard statistics');
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100vh' }}>
        <CircularProgress sx={{ color: novikTheme.colors.primary }} size={48} />
      </Box>
    );
  }

  if (error) {
    return (
      <Container maxWidth="sm" sx={{ mt: 8 }}>
        <Alert severity="error" sx={{ borderRadius: '12px' }}>
          <Typography variant="subtitle1" fontWeight={600}>
            Error Loading Dashboard
          </Typography>
          <Typography variant="body2">{error}</Typography>
        </Alert>
      </Container>
    );
  }

  if (!stats) return null;

  const metrics = [
    {
      title: 'Total Users',
      value: stats.userStats.total,
      change: `+${stats.userStats.newThisMonth} this month`,
      icon: <PeopleIcon />,
      color: '#4a90e2',
      bgColor: 'rgba(74, 144, 226, 0.1)',
    },
    {
      title: 'Active Users',
      value: stats.userStats.active,
      change: `${((stats.userStats.active / stats.userStats.total) * 100).toFixed(1)}% of total`,
      icon: <TrendingUpIcon />,
      color: novikTheme.colors.primary,
      bgColor: 'rgba(136, 169, 78, 0.1)',
    },
    {
      title: 'Total Conversations',
      value: stats.conversationStats.total,
      change: `${stats.conversationStats.avgPerUser.toFixed(1)} per user`,
      icon: <ChatIcon />,
      color: '#6c63ff',
      bgColor: 'rgba(108, 99, 255, 0.1)',
    },
    {
      title: 'New This Week',
      value: stats.userStats.newThisWeek,
      change: 'User registrations',
      icon: <PersonAddIcon />,
      color: '#00bcd4',
      bgColor: 'rgba(0, 188, 212, 0.1)',
    },
  ];

  return (
    <PageContainer>
      {/* Header */}
      <HeaderSection>
        <Container maxWidth="xl">
          <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <Box>
              <Typography
                variant="h4"
                sx={{
                  fontWeight: 700,
                  color: novikTheme.colors.text,
                  fontFamily: novikTheme.typography.fontFamily,
                }}
              >
                Admin Dashboard
              </Typography>
              <Typography
                variant="body2"
                sx={{
                  color: novikTheme.colors.textMuted,
                  mt: 0.5,
                  fontFamily: novikTheme.typography.fontFamily,
                }}
              >
                Monitor and manage your application performance
              </Typography>
            </Box>
            <Button
              variant="contained"
              startIcon={<GroupIcon />}
              onClick={() => navigate('/admin/users')}
              sx={{
                backgroundColor: novikTheme.colors.primary,
                color: '#ffffff',
                fontWeight: 600,
                padding: '0.5rem 1.5rem',
                borderRadius: '8px',
                textTransform: 'none',
                fontFamily: novikTheme.typography.fontFamily,
                '&:hover': {
                  backgroundColor: novikTheme.colors.primaryDark,
                },
              }}
            >
              Manage Users
            </Button>
          </Box>
        </Container>
      </HeaderSection>

      <Container maxWidth="xl" sx={{ py: 4 }}>
        {/* Key Metrics */}
        <Grid container spacing={3} sx={{ mb: 4 }}>
          {metrics.map((metric, index) => (
            <Grid size={{ xs: 12, sm: 6, md: 3 }} key={index}>
              <MetricCard>
                <CardContent>
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                    <Box sx={{ flex: 1 }}>
                      <Typography
                        variant="body2"
                        sx={{
                          color: novikTheme.colors.textMuted,
                          fontWeight: 500,
                          mb: 1,
                          fontFamily: novikTheme.typography.fontFamily,
                        }}
                      >
                        {metric.title}
                      </Typography>
                      <Typography
                        variant="h4"
                        sx={{
                          fontWeight: 700,
                          color: novikTheme.colors.text,
                          mb: 0.5,
                          fontFamily: novikTheme.typography.fontFamily,
                        }}
                      >
                        {metric.value.toLocaleString()}
                      </Typography>
                      <Typography
                        variant="caption"
                        sx={{
                          color: novikTheme.colors.textMuted,
                          fontFamily: novikTheme.typography.fontFamily,
                        }}
                      >
                        {metric.change}
                      </Typography>
                    </Box>
                    <MetricIconBox bgcolor={metric.bgColor}>
                      <Box sx={{ color: metric.color, display: 'flex' }}>{metric.icon}</Box>
                    </MetricIconBox>
                  </Box>
                </CardContent>
              </MetricCard>
            </Grid>
          ))}
        </Grid>

        {/* Data Grids */}
        <Grid container spacing={3}>
          {/* Users by Country */}
          <Grid size={{ xs: 12, md: 4 }}>
            <DataCard>
              <CardContent>
                <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.text,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Geographic Distribution
                  </Typography>
                  <PublicIcon sx={{ color: novikTheme.colors.textMuted }} />
                </Box>
                <Box sx={{ mt: 2 }}>
                  {stats.usersByCountry.slice(0, 5).map((country, index) => {
                    const percentage = (country.count / stats.userStats.total) * 100;
                    return (
                      <CountryRow key={index}>
                        <Box sx={{ flex: 1 }}>
                          <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 0.5 }}>
                            <Typography
                              variant="body2"
                              sx={{
                                fontWeight: 500,
                                color: novikTheme.colors.text,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {country.country}
                            </Typography>
                            <Typography
                              variant="body2"
                              sx={{
                                color: novikTheme.colors.textMuted,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {country.count} ({percentage.toFixed(1)}%)
                            </Typography>
                          </Box>
                          <ProgressBar variant="determinate" value={percentage} />
                        </Box>
                      </CountryRow>
                    );
                  })}
                </Box>
              </CardContent>
            </DataCard>
          </Grid>

          {/* Recent Registrations */}
          <Grid size={{ xs: 12, md: 4 }}>
            <DataCard>
              <CardContent>
                <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.text,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Recent Registrations
                  </Typography>
                  <PersonAddIcon sx={{ color: novikTheme.colors.textMuted }} />
                </Box>
                <Box>
                  {stats.recentRegistrations.map(user => (
                    <UserListItem key={user.id}>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1.5 }}>
                          <Avatar
                            sx={{
                              width: 32,
                              height: 32,
                              bgcolor: novikTheme.colors.primary,
                              fontSize: '0.875rem',
                            }}
                          >
                            {user.username.charAt(0).toUpperCase()}
                          </Avatar>
                          <Box>
                            <Typography
                              variant="body2"
                              sx={{
                                fontWeight: 500,
                                color: novikTheme.colors.text,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {user.username}
                            </Typography>
                            <Typography
                              variant="caption"
                              sx={{
                                color: novikTheme.colors.textMuted,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {user.email}
                            </Typography>
                          </Box>
                        </Box>
                        <Typography
                          variant="caption"
                          sx={{
                            color: novikTheme.colors.textMuted,
                            fontFamily: novikTheme.typography.fontFamily,
                          }}
                        >
                          {new Date(user.dateJoined).toLocaleDateString('en-US', {
                            month: 'short',
                            day: 'numeric',
                          })}
                        </Typography>
                      </Box>
                    </UserListItem>
                  ))}
                </Box>
              </CardContent>
            </DataCard>
          </Grid>

          {/* Top Active Users */}
          <Grid size={{ xs: 12, md: 4 }}>
            <DataCard>
              <CardContent>
                <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.text,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Most Active Users
                  </Typography>
                  <StarIcon sx={{ color: '#ffd700' }} />
                </Box>
                <Box>
                  {stats.topUsers.map((user, index) => (
                    <UserListItem key={user.id}>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1.5 }}>
                          <Box
                            sx={{
                              width: 32,
                              height: 32,
                              borderRadius: '8px',
                              display: 'flex',
                              alignItems: 'center',
                              justifyContent: 'center',
                              bgcolor: index === 0 ? '#ffd700' : index === 1 ? '#c0c0c0' : index === 2 ? '#cd7f32' : '#f0f0f0',
                              color: index < 3 ? '#fff' : novikTheme.colors.textMuted,
                              fontWeight: 700,
                            }}
                          >
                            {index + 1}
                          </Box>
                          <Box>
                            <Typography
                              variant="body2"
                              sx={{
                                fontWeight: 500,
                                color: novikTheme.colors.text,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {user.username}
                            </Typography>
                            <Typography
                              variant="caption"
                              sx={{
                                color: novikTheme.colors.textMuted,
                                fontFamily: novikTheme.typography.fontFamily,
                              }}
                            >
                              {user.email}
                            </Typography>
                          </Box>
                        </Box>
                        <Box
                          sx={{
                            px: 1.5,
                            py: 0.5,
                            borderRadius: '16px',
                            bgcolor: 'rgba(136, 169, 78, 0.1)',
                          }}
                        >
                          <Typography
                            variant="caption"
                            sx={{
                              fontWeight: 600,
                              color: novikTheme.colors.primary,
                              fontFamily: novikTheme.typography.fontFamily,
                            }}
                          >
                            {user.conversationCount} chats
                          </Typography>
                        </Box>
                      </Box>
                    </UserListItem>
                  ))}
                </Box>
              </CardContent>
            </DataCard>
          </Grid>
        </Grid>

        {/* Summary Statistics Footer */}
        <StatsFooter elevation={0}>
          <Grid container spacing={4}>
            <Grid size={{ xs: 12, sm: 6, md: 3 }}>
              <Box>
                <Typography
                  variant="body2"
                  sx={{
                    color: 'rgba(255,255,255,0.8)',
                    mb: 0.5,
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  Inactive Users
                </Typography>
                <Typography
                  variant="h5"
                  sx={{
                    fontWeight: 700,
                    color: '#ffffff',
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  {stats.userStats.inactive}
                </Typography>
              </Box>
            </Grid>
            <Grid size={{ xs: 12, sm: 6, md: 3 }}>
              <Box>
                <Typography
                  variant="body2"
                  sx={{
                    color: 'rgba(255,255,255,0.8)',
                    mb: 0.5,
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  Conversations This Month
                </Typography>
                <Typography
                  variant="h5"
                  sx={{
                    fontWeight: 700,
                    color: '#ffffff',
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  {stats.conversationStats.thisMonth}
                </Typography>
              </Box>
            </Grid>
            <Grid size={{ xs: 12, sm: 6, md: 3 }}>
              <Box>
                <Typography
                  variant="body2"
                  sx={{
                    color: 'rgba(255,255,255,0.8)',
                    mb: 0.5,
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  Countries Reached
                </Typography>
                <Typography
                  variant="h5"
                  sx={{
                    fontWeight: 700,
                    color: '#ffffff',
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  {stats.usersByCountry.length}
                </Typography>
              </Box>
            </Grid>
            <Grid size={{ xs: 12, sm: 6, md: 3 }}>
              <Box>
                <Typography
                  variant="body2"
                  sx={{
                    color: 'rgba(255,255,255,0.8)',
                    mb: 0.5,
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  Avg Conversations/User
                </Typography>
                <Typography
                  variant="h5"
                  sx={{
                    fontWeight: 700,
                    color: '#ffffff',
                    fontFamily: novikTheme.typography.fontFamily,
                  }}
                >
                  {stats.conversationStats.avgPerUser.toFixed(1)}
                </Typography>
              </Box>
            </Grid>
          </Grid>
        </StatsFooter>
      </Container>
    </PageContainer>
  );
}
