import { useState, useEffect, ChangeEvent, FormEvent } from 'react';
import {
  Box,
  Container,
  Typography,
  Button,
  TextField,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Chip,
  IconButton,
  FormControlLabel,
  Switch,
  Card,
  CardContent,
  CircularProgress,
  Alert,
  styled,
} from '@mui/material';
import {
  Edit as EditIcon,
  Delete as DeleteIcon,
  Add as AddIcon,
  BarChart as StatsIcon,
  CloudUpload as CloudUploadIcon,
  Image as ImageIcon,
} from '@mui/icons-material';
import { novikTheme } from '../styles/theme';

type Banner = {
  id: number;
  title: string;
  image: string | null;
  image_url: string | null;
  link: string | null;
  code: string | null;
  is_active: boolean;
};
type BannerStat = {
  id: number;
  banner: number;
  date: string;
  country: string;
  views: number;
  clicks: number;
};

// Styled Components
const PageContainer = styled(Container)({
  paddingTop: '2rem',
  paddingBottom: '2rem',
});

const HeaderSection = styled(Box)({
  marginBottom: '2rem',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'space-between',
});

const StyledCard = styled(Card)({
  marginBottom: '2rem',
  borderRadius: '12px',
  boxShadow: '0 2px 8px rgba(0,0,0,0.08)',
  '&:hover': {
    boxShadow: '0 4px 12px rgba(0,0,0,0.12)',
  },
  transition: 'box-shadow 0.3s ease',
});

const FormSection = styled(CardContent)({
  padding: '2rem',
});

const StyledTextField = styled(TextField)({
  '& .MuiOutlinedInput-root': {
    '&:hover fieldset': {
      borderColor: novikTheme.colors.primary,
    },
    '&.Mui-focused fieldset': {
      borderColor: novikTheme.colors.primary,
    },
  },
  '& .MuiInputLabel-root.Mui-focused': {
    color: novikTheme.colors.primary,
  },
});

const PrimaryButton = styled(Button)({
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
});

const SecondaryButton = styled(Button)({
  color: novikTheme.colors.primary,
  borderColor: novikTheme.colors.primary,
  fontWeight: 600,
  padding: '0.4rem 1rem',
  borderRadius: '8px',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: 'rgba(136, 169, 78, 0.08)',
    borderColor: novikTheme.colors.primaryDark,
  },
});

const UploadButton = styled(Button)<any>({
  backgroundColor: 'rgba(136, 169, 78, 0.1)',
  color: novikTheme.colors.primary,
  border: `2px dashed ${novikTheme.colors.primary}`,
  padding: '1rem',
  borderRadius: '8px',
  textTransform: 'none',
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: 'rgba(136, 169, 78, 0.15)',
    borderColor: novikTheme.colors.primaryDark,
  },
});

const StyledTableCell = styled(TableCell)({
  fontFamily: novikTheme.typography.fontFamily,
});

const StyledTableHeadCell = styled(TableCell)({
  fontFamily: novikTheme.typography.fontFamily,
  fontWeight: 600,
  backgroundColor: '#f8f9fa',
  color: novikTheme.colors.text,
});

const ActiveChip = styled(Chip)<{ active?: boolean }>(({ active }) => ({
  backgroundColor: active ? 'rgba(136, 169, 78, 0.1)' : 'rgba(0, 0, 0, 0.08)',
  color: active ? novikTheme.colors.primary : '#666',
  fontWeight: 500,
  fontFamily: novikTheme.typography.fontFamily,
}));

const StyledSwitch = styled(Switch)({
  '& .MuiSwitch-switchBase.Mui-checked': {
    color: novikTheme.colors.primary,
  },
  '& .MuiSwitch-switchBase.Mui-checked + .MuiSwitch-track': {
    backgroundColor: novikTheme.colors.primary,
  },
});

export default function BannerManagement() {
  const [banners, setBanners] = useState<Banner[]>([]);
  const [stats, setStats] = useState<BannerStat[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [form, setForm] = useState<Partial<Banner>>({});
  const [file, setFile] = useState<File | null>(null);
  const [editingId, setEditingId] = useState<number | null>(null);

  // fetch list
  const fetchBanners = async () => {
    try {
      const res = await fetch('/api/banners/', { credentials: 'include' });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      setBanners(data);
    } catch (err: any) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  // fetch stats for one banner
  const fetchStats = async (bid: number) => {
    try {
      const res = await fetch(`/api/banner-stats/?banner=${bid}`, {
        credentials: 'include',
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      setStats(await res.json());
    } catch (err: any) {
      console.error(err);
    }
  };

  useEffect(() => {
    fetchBanners();
  }, []);

  // form handlers
  const handleInput = (e: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value, type, checked } = e.target as HTMLInputElement;
    setForm(prev => ({
      ...prev,
      [name]: type === 'checkbox' ? checked : value,
    }));
  };
  const handleFile = (e: ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) setFile(e.target.files[0]);
  };

  // create/update
  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const fd = new FormData();

    // Add form fields to FormData
    if (form.title) fd.append('title', form.title);
    if (form.link) fd.append('link', form.link);
    if (form.code) fd.append('code', form.code);
    if (form.is_active !== undefined) fd.append('is_active', form.is_active.toString());

    // Add file if it exists
    if (file) {
      fd.append('image', file);
    }

    const url = editingId ? `/api/banners/${editingId}/` : '/api/banners/';
    const method = editingId ? 'PUT' : 'POST';

    try {
      const res = await fetch(url, {
        method,
        credentials: 'include',
        body: fd,
        // Don't set Content-Type header - browser will set it with boundary
      });

      if (!res.ok) {
        const errorData = await res.json();
        throw new Error(errorData.detail || `HTTP ${res.status}`);
      }

      // Reset form
      setForm({});
      setFile(null);
      setEditingId(null);
      fetchBanners();
    } catch (err: any) {
      setError(err.message);
    }
  };

  const handleEdit = (b: Banner) => {
    setEditingId(b.id);
    setForm({
      title: b.title,
      link: b.link || '',
      code: b.code || '',
      is_active: b.is_active,
    });
  };

  const handleDelete = async (id: number) => {
    if (!confirm('Delete this banner?')) return;
    try {
      const res = await fetch(`/api/banners/${id}/`, {
        method: 'DELETE',
        credentials: 'include',
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      fetchBanners();
    } catch (err: any) {
      setError(err.message);
    }
  };

  if (loading)
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100vh' }}>
        <CircularProgress sx={{ color: novikTheme.colors.primary }} />
      </Box>
    );

  if (error)
    return (
      <Container maxWidth="lg" sx={{ mt: 4 }}>
        <Alert severity="error" sx={{ borderRadius: '8px' }}>
          Error: {error}
        </Alert>
      </Container>
    );

  return (
    <PageContainer maxWidth="lg">
      {/* Header Section */}
      <HeaderSection>
        <Typography
          variant="h4"
          sx={{
            fontWeight: 600,
            color: novikTheme.colors.text,
            fontFamily: novikTheme.typography.fontFamily,
          }}
        >
          Banner Management
        </Typography>
      </HeaderSection>

      {/* Form Section */}
      <StyledCard>
        <FormSection>
          <Typography
            variant="h5"
            sx={{
              fontWeight: 600,
              mb: 3,
              color: novikTheme.colors.text,
              fontFamily: novikTheme.typography.fontFamily,
            }}
          >
            {editingId ? 'Edit Banner' : 'Create New Banner'}
          </Typography>

          <Box component="form" onSubmit={handleSubmit}>
            <Box sx={{ display: 'grid', gap: 3, mb: 3, gridTemplateColumns: { md: 'repeat(2, 1fr)' } }}>
              <StyledTextField
                name="title"
                label="Title"
                value={form.title || ''}
                onChange={handleInput}
                required
                fullWidth
                variant="outlined"
              />

              <StyledTextField
                name="link"
                label="Link URL"
                value={form.link || ''}
                onChange={handleInput}
                placeholder="https://"
                fullWidth
                variant="outlined"
              />
            </Box>

            <StyledTextField
              name="code"
              label="HTML Code"
              value={form.code || ''}
              onChange={handleInput}
              multiline
              rows={4}
              fullWidth
              variant="outlined"
              placeholder="<div>Your banner HTML here</div>"
              sx={{ mb: 3 }}
            />

            <Box sx={{ mb: 3 }}>
              <Typography
                variant="body1"
                sx={{
                  mb: 2,
                  color: novikTheme.colors.text,
                  fontWeight: 500,
                  fontFamily: novikTheme.typography.fontFamily,
                }}
              >
                Banner Image
              </Typography>
              <input
                accept="image/*"
                style={{ display: 'none' }}
                id="raised-button-file"
                type="file"
                onChange={handleFile}
              />
              <label htmlFor="raised-button-file">
                <UploadButton
                  component="span"
                  fullWidth
                  startIcon={<CloudUploadIcon />}
                >
                  {file ? file.name : 'Choose Image File'}
                </UploadButton>
              </label>
            </Box>

            <FormControlLabel
              control={
                <StyledSwitch
                  name="is_active"
                  checked={form.is_active || false}
                  onChange={handleInput}
                />
              }
              label="Active"
              sx={{
                mb: 3,
                '& .MuiFormControlLabel-label': {
                  fontFamily: novikTheme.typography.fontFamily,
                  color: novikTheme.colors.text,
                },
              }}
            />

            <Box sx={{ display: 'flex', justifyContent: 'flex-end', gap: 2 }}>
              {editingId && (
                <SecondaryButton
                  variant="outlined"
                  onClick={() => {
                    setEditingId(null);
                    setForm({});
                    setFile(null);
                  }}
                >
                  Cancel
                </SecondaryButton>
              )}
              <PrimaryButton
                type="submit"
                variant="contained"
                startIcon={editingId ? <EditIcon /> : <AddIcon />}
              >
                {editingId ? 'Update Banner' : 'Create Banner'}
              </PrimaryButton>
            </Box>
          </Box>
        </FormSection>
      </StyledCard>

      {/* Banners List */}
      <StyledCard>
        <CardContent sx={{ p: 0 }}>
          <Box sx={{ p: 3, borderBottom: '1px solid #e0e0e0', backgroundColor: '#f8f9fa' }}>
            <Typography
              variant="h5"
              sx={{
                fontWeight: 600,
                color: novikTheme.colors.text,
                fontFamily: novikTheme.typography.fontFamily,
              }}
            >
              Manage Banners
            </Typography>
          </Box>

          <TableContainer>
            <Table>
              <TableHead>
                <TableRow>
                  <StyledTableHeadCell>Title</StyledTableHeadCell>
                  <StyledTableHeadCell>Preview</StyledTableHeadCell>
                  <StyledTableHeadCell>Status</StyledTableHeadCell>
                  <StyledTableHeadCell align="center">Actions</StyledTableHeadCell>
                  <StyledTableHeadCell align="center">Statistics</StyledTableHeadCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {banners.length === 0 ? (
                  <TableRow>
                    <StyledTableCell colSpan={5} align="center">
                      <Typography
                        sx={{
                          py: 4,
                          color: novikTheme.colors.textMuted,
                          fontFamily: novikTheme.typography.fontFamily,
                        }}
                      >
                        No banners found. Create your first banner above.
                      </Typography>
                    </StyledTableCell>
                  </TableRow>
                ) : (
                  banners.map(b => (
                    <TableRow key={b.id} hover>
                      <StyledTableCell>{b.title}</StyledTableCell>
                      <StyledTableCell>
                        {b.image_url ? (
                          <Box
                            component="img"
                            src={b.image_url}
                            alt={b.title}
                            sx={{
                              height: 80,
                              width: 'auto',
                              objectFit: 'contain',
                              borderRadius: '8px',
                              border: '1px solid #e0e0e0',
                            }}
                          />
                        ) : b.code ? (
                          <Box
                            sx={{
                              maxWidth: 300,
                              overflow: 'hidden',
                              p: 1,
                              border: '1px solid #e0e0e0',
                              borderRadius: '8px',
                              backgroundColor: '#f8f9fa',
                            }}
                            dangerouslySetInnerHTML={{ __html: b.code }}
                          />
                        ) : (
                          <Box
                            sx={{
                              display: 'flex',
                              alignItems: 'center',
                              justifyContent: 'center',
                              height: 80,
                              width: 120,
                              backgroundColor: '#f5f5f5',
                              borderRadius: '8px',
                              border: '1px solid #e0e0e0',
                            }}
                          >
                            <ImageIcon sx={{ color: '#999' }} />
                          </Box>
                        )}
                      </StyledTableCell>
                      <StyledTableCell>
                        <ActiveChip
                          label={b.is_active ? 'Active' : 'Inactive'}
                          active={b.is_active}
                          size="small"
                        />
                      </StyledTableCell>
                      <StyledTableCell align="center">
                        <Box sx={{ display: 'flex', gap: 1, justifyContent: 'center' }}>
                          <IconButton
                            onClick={() => handleEdit(b)}
                            sx={{
                              color: novikTheme.colors.primary,
                              '&:hover': {
                                backgroundColor: 'rgba(136, 169, 78, 0.08)',
                              },
                            }}
                          >
                            <EditIcon fontSize="small" />
                          </IconButton>
                          <IconButton
                            onClick={() => handleDelete(b.id)}
                            sx={{
                              color: '#ef4444',
                              '&:hover': {
                                backgroundColor: 'rgba(239, 68, 68, 0.08)',
                              },
                            }}
                          >
                            <DeleteIcon fontSize="small" />
                          </IconButton>
                        </Box>
                      </StyledTableCell>
                      <StyledTableCell align="center">
                        <IconButton
                          onClick={() => fetchStats(b.id)}
                          sx={{
                            color: '#3b82f6',
                            '&:hover': {
                              backgroundColor: 'rgba(59, 130, 246, 0.08)',
                            },
                          }}
                        >
                          <StatsIcon fontSize="small" />
                        </IconButton>
                      </StyledTableCell>
                    </TableRow>
                  ))
                )}
              </TableBody>
            </Table>
          </TableContainer>
        </CardContent>
      </StyledCard>

      {/* Stats Section */}
      {stats.length > 0 && (
        <StyledCard>
          <CardContent sx={{ p: 0 }}>
            <Box sx={{ p: 3, borderBottom: '1px solid #e0e0e0', backgroundColor: '#f8f9fa' }}>
              <Typography
                variant="h5"
                sx={{
                  fontWeight: 600,
                  color: novikTheme.colors.text,
                  fontFamily: novikTheme.typography.fontFamily,
                }}
              >
                Banner Statistics
              </Typography>
            </Box>

            <TableContainer>
              <Table>
                <TableHead>
                  <TableRow>
                    <StyledTableHeadCell>Date</StyledTableHeadCell>
                    <StyledTableHeadCell>Country</StyledTableHeadCell>
                    <StyledTableHeadCell align="center">Views</StyledTableHeadCell>
                    <StyledTableHeadCell align="center">Clicks</StyledTableHeadCell>
                    <StyledTableHeadCell align="center">CTR</StyledTableHeadCell>
                  </TableRow>
                </TableHead>
                <TableBody>
                  {stats.map(s => (
                    <TableRow key={s.id} hover>
                      <StyledTableCell>
                        {new Date(s.date).toLocaleDateString('en-US', {
                          year: 'numeric',
                          month: 'short',
                          day: 'numeric',
                        })}
                      </StyledTableCell>
                      <StyledTableCell>
                        <Chip
                          label={s.country}
                          size="small"
                          sx={{
                            backgroundColor: 'rgba(136, 169, 78, 0.08)',
                            color: novikTheme.colors.text,
                            fontFamily: novikTheme.typography.fontFamily,
                          }}
                        />
                      </StyledTableCell>
                      <StyledTableCell align="center">
                        <Typography
                          sx={{
                            fontWeight: 500,
                            color: novikTheme.colors.text,
                            fontFamily: novikTheme.typography.fontFamily,
                          }}
                        >
                          {s.views.toLocaleString()}
                        </Typography>
                      </StyledTableCell>
                      <StyledTableCell align="center">
                        <Typography
                          sx={{
                            fontWeight: 500,
                            color: novikTheme.colors.primary,
                            fontFamily: novikTheme.typography.fontFamily,
                          }}
                        >
                          {s.clicks.toLocaleString()}
                        </Typography>
                      </StyledTableCell>
                      <StyledTableCell align="center">
                        <Typography
                          sx={{
                            fontWeight: 600,
                            color: s.views > 0 ? novikTheme.colors.primary : novikTheme.colors.textMuted,
                            fontFamily: novikTheme.typography.fontFamily,
                          }}
                        >
                          {s.views > 0 ? `${((s.clicks / s.views) * 100).toFixed(2)}%` : '0%'}
                        </Typography>
                      </StyledTableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </TableContainer>

            {/* Stats Summary */}
            <Box sx={{ p: 3, backgroundColor: '#f8f9fa', borderTop: '1px solid #e0e0e0' }}>
              <Box sx={{ display: 'flex', gap: 4, flexWrap: 'wrap' }}>
                <Box>
                  <Typography
                    variant="body2"
                    sx={{
                      color: novikTheme.colors.textMuted,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Total Views
                  </Typography>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.text,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    {stats.reduce((sum, s) => sum + s.views, 0).toLocaleString()}
                  </Typography>
                </Box>
                <Box>
                  <Typography
                    variant="body2"
                    sx={{
                      color: novikTheme.colors.textMuted,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Total Clicks
                  </Typography>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.primary,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    {stats.reduce((sum, s) => sum + s.clicks, 0).toLocaleString()}
                  </Typography>
                </Box>
                <Box>
                  <Typography
                    variant="body2"
                    sx={{
                      color: novikTheme.colors.textMuted,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    Average CTR
                  </Typography>
                  <Typography
                    variant="h6"
                    sx={{
                      fontWeight: 600,
                      color: novikTheme.colors.primary,
                      fontFamily: novikTheme.typography.fontFamily,
                    }}
                  >
                    {(() => {
                      const totalViews = stats.reduce((sum, s) => sum + s.views, 0);
                      const totalClicks = stats.reduce((sum, s) => sum + s.clicks, 0);
                      return totalViews > 0 ? `${((totalClicks / totalViews) * 100).toFixed(2)}%` : '0%';
                    })()}
                  </Typography>
                </Box>
              </Box>
            </Box>
          </CardContent>
        </StyledCard>
      )}
    </PageContainer>
  );
}
