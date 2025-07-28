import { useState, useEffect, useRef, ChangeEvent } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkHeaderId from 'remark-heading-id';
import remarkGfm from 'remark-gfm';
import { v4 as uuidv4 } from 'uuid';
import { patientService } from '../services/patientService';
import type { Message, Conversation } from '../types';
import {
  Box,
  Container,
  Typography,
  IconButton,
  TextareaAutosize,
  Fab,
  CircularProgress,
  useTheme,
  Button,
  Alert,
  Snackbar,
  Tooltip,
  Divider,
  Card,
  CardContent,
  Chip,
  TextField,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  AppBar,
  Toolbar,
} from '@mui/material';
import {
  AttachFile,
  Send,
  DeleteOutline,
  ContentCopy,
  FileDownload,
  PictureAsPdf,
  Clear,
  Add,
  Menu as MenuIcon,
} from '@mui/icons-material';
import { useConverstaions } from '../hooks/useConversations';

type AlertType = {
  open: boolean;
  message: string;
  severity: 'error' | 'warning' | 'info' | 'success';
};

const drawerWidth = 280;

const Dashboard = () => {
  const theme = useTheme();
  const [input, setInput] = useState('');
  const [selectedFile, setSelectedFile] = useState<File | undefined>(undefined);
  const [loading, setLoading] = useState(false);
  const [pendingQuestion, setPendingQuestion] = useState<string | null>(null);
  const bottomRef = useRef<HTMLDivElement | null>(null);
  const textAreaRef = useRef<HTMLTextAreaElement | null>(null);
  const [alert, setAlert] = useState<AlertType>({
    open: false,
    message: '',
    severity: 'info',
  });
  const [copiedMessageId, setCopiedMessageId] = useState<string | null>(null);

  const { mutate } = useConverstaions();
  const [selectedConversation, setSelectedConversation] = useState<Conversation | null>(null);
  const [dialogOpen, setDialogOpen] = useState(false);
  const [dialogTitle, setDialogTitle] = useState('');
  const [dialogMode, setDialogMode] = useState<'create' | 'rename'>('create');
  const [dialogConversationId, setDialogConversationId] = useState<string | null>(null);
  const [mobileOpen, setMobileOpen] = useState(false);

  useEffect(() => {
    if (copiedMessageId) {
      const timer = setTimeout(() => {
        setCopiedMessageId(null);
      }, 2000);
      return () => clearTimeout(timer);
    }
  }, [copiedMessageId]);

  const scrollToBottom = (behavior: ScrollBehavior = 'smooth') => {
    bottomRef.current?.scrollIntoView({ behavior });
  };

  const handleAlertClose = () => {
    setAlert(prev => ({ ...prev, open: false }));
  };

  const showAlert = (message: string, severity: AlertType['severity'] = 'info') => {
    setAlert({ open: true, message, severity });
  };

  const handleFileSelect = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      if (file.type === 'application/pdf') {
        setSelectedFile(file);
        showAlert(`File "${file.name}" selected`, 'info');
      } else {
        showAlert('Only PDF files are supported', 'error');
      }
    }
  };

  const clearChat = () => {
    if (!selectedConversation || selectedConversation.messages.length === 0) return;

    if (confirm('Are you sure you want to clear the conversation history?')) {
      handleClearConversation(selectedConversation.id);
    }
  };

  const copyToClipboard = (text: string, id: string) => {
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(
        () => {
          setCopiedMessageId(id);
        },
        err => {
          console.error('Could not copy text: ', err);
          showAlert('Failed to copy to clipboard', 'error');
        }
      );
    } else {
      try {
        const textArea = document.createElement('textarea');
        textArea.value = text;
        textArea.style.position = 'fixed';
        document.body.appendChild(textArea);
        textArea.focus();
        textArea.select();

        const successful = document.execCommand('copy');
        document.body.removeChild(textArea);

        if (successful) {
          setCopiedMessageId(id);
        } else {
          showAlert('Failed to copy to clipboard', 'error');
        }
      } catch (err) {
        console.error('Fallback: Could not copy text: ', err);
        showAlert('Failed to copy to clipboard', 'error');
      }
    }
  };

  const exportChat = () => {
    if (!selectedConversation || selectedConversation.messages.length === 0) {
      showAlert('No conversation to export', 'warning');
      return;
    }

    const exportDate = new Date().toISOString().slice(0, 10);
    const fileName = `dental-ai-chat_${exportDate}.txt`;

    const content = selectedConversation.messages
      .map(msg => {
        const prefix = msg.role === 'user' ? 'Q: ' : 'A: ';
        return `${prefix}${msg.content}\n\n`;
      })
      .join('${"=".repeat(40)}\n\n');

    const element = document.createElement('a');
    const file = new Blob([content], { type: 'text/plain' });
    element.href = URL.createObjectURL(file);
    element.download = fileName;
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
    showAlert('Conversation exported successfully', 'success');
  };

  const handleSubmit = async () => {
    if (!input.trim() && !selectedFile) return;
    if (!selectedConversation) {
      // Create a new conversation if none is selected
      await handleCreateNewConversation(true);
      return;
    }

    const messageText = selectedFile ? `${input}` : input.trim();

    setPendingQuestion(messageText);
    setLoading(true);
    setInput('');
    scrollToBottom();

    try {
      const response = await patientService.ask(messageText, selectedFile);

      // Update the UI immediately for better UX
      const updatedConversation: Conversation = {
        ...selectedConversation,
        messages: [
          ...selectedConversation.messages,
          {
            id: uuidv4(),
            role: 'user',
            content: messageText,
            file: selectedFile && { fileName: selectedFile.name, text: '' },
          },
          {
            id: uuidv4(),
            role: 'assistant',
            content: response.message,
          },
        ],
      };

      setSelectedConversation(updatedConversation);

      // Update the conversations list to reflect the change
      mutate(convs =>
        convs?.map(conv => (conv.id === selectedConversation.id ? updatedConversation : conv))
      );

      setSelectedFile(undefined);
      scrollToBottom();
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Network error. Please try again.';
      showAlert('Error: ' + errorMessage, 'error');
    } finally {
      setPendingQuestion(null);
      setLoading(false);
    }
  };

  const handleKeyPress = (e: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit();
    }
  };

  const handleCreateNewConversation = async (selectAfterCreate = false) => {
    try {
      const title = dialogTitle || 'New Conversation';
      const newConversation = await patientService.createConversation(title);

      mutate(data => [...(data || []), newConversation], false);

      if (selectAfterCreate) {
        setSelectedConversation(newConversation);
      }

      setDialogOpen(false);
      setDialogTitle('');
      showAlert('Conversation created', 'success');
    } catch (err) {
      showAlert('Failed to create conversation', 'error');
      console.error('Error creating conversation:', err);
    }
  };

  const handleRenameConversation = async () => {
    if (!dialogConversationId) return;

    try {
      const updatedConversation = await patientService.updateConversationTitle(
        dialogConversationId,
        dialogTitle
      );

      mutate(
        convs =>
          convs?.map(conv =>
            conv.id === dialogConversationId ? { ...conv, title: updatedConversation.title } : conv
          ),
        false
      );

      if (selectedConversation?.id === dialogConversationId) {
        setSelectedConversation({ ...selectedConversation, title: updatedConversation.title });
      }

      setDialogOpen(false);
      setDialogTitle('');
      setDialogConversationId(null);
      showAlert('Conversation renamed', 'success');
    } catch (err) {
      showAlert('Failed to rename conversation', 'error');
      console.error('Error renaming conversation:', err);
    }
  };

  const handleClearConversation = async (id: string) => {
    try {
      await patientService.clearConversation(id);

      mutate(
        convs => convs?.map(conv => (conv.id === id ? { ...conv, messages: [] } : conv)),
        false
      );

      if (selectedConversation?.id === id) {
        setSelectedConversation(prev => (prev ? { ...prev, messages: [] } : null));
      }

      showAlert('Conversation cleared', 'success');
    } catch (err) {
      showAlert('Failed to clear conversation', 'error');
      console.error('Error clearing conversation:', err);
    }
  };

  const handleDialogClose = () => {
    setDialogOpen(false);
    setDialogTitle('');
    setDialogConversationId(null);
  };

  const handleDialogSubmit = () => {
    if (dialogMode === 'create') {
      handleCreateNewConversation(true);
    } else {
      handleRenameConversation();
    }
  };

  const handleDrawerToggle = () => {
    setMobileOpen(!mobileOpen);
  };

  const getMessages = (): Message[] => {
    if (!selectedConversation) return [];

    return selectedConversation.messages.map((msg, index) => ({
      id: String(index),
      role: msg.role,
      content: msg.content,
      timestamp: new Date(selectedConversation.updatedAt),
      file: msg.file,
    }));
  };

  const messages = getMessages();

  return (
    <Box sx={{ display: 'flex' }}>
      {/* AppBar for mobile */}
      <AppBar
        position="fixed"
        sx={{
          display: { sm: 'none' },
          zIndex: theme.zIndex.drawer + 1,
        }}
      >
        <Toolbar>
          <IconButton color="inherit" edge="start" onClick={handleDrawerToggle} sx={{ mr: 2 }}>
            <MenuIcon />
          </IconButton>
          <Typography variant="h6" noWrap component="div">
            AI Dental Assistant
          </Typography>
        </Toolbar>
      </AppBar>

      {/* Main content */}
      <Box
        component="main"
        sx={{
          flexGrow: 1,
          p: 0,
          width: { sm: `calc(100% - ${drawerWidth}px)` },
          height: `calc(100vh - 64px)`,
          display: 'flex',
          flexDirection: 'column',
          mt: { xs: 7, sm: 0 },
        }}
      >
        {/* Header */}
        <Box
          sx={{
            borderBottom: 1,
            borderColor: 'divider',
            p: 2,
          }}
        >
          <Container maxWidth="md">
            <Box
              sx={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
              }}
            >
              {selectedConversation ? (
                <Typography variant="h5" color="#f97316">
                  {selectedConversation.title}
                </Typography>
              ) : (
                <Box sx={{ display: 'flex', alignItems: 'center' }}>
                  <Typography variant="h5" color="text.secondary">
                    AI Dental Assistant
                  </Typography>
                  {/* <Chip
                    label="Select a conversation"
                    size="small"
                    variant="outlined"
                    sx={{ ml: 2, color: '#f97316', borderColor: '#f97316' }}
                  /> */}
                </Box>
              )}

              <Box>
                <Tooltip title="Clear conversation">
                  <span>
                    {/* Wrap in span to allow tooltip on disabled button */}
                    <Button
                      startIcon={<DeleteOutline />}
                      onClick={clearChat}
                      disabled={!selectedConversation || selectedConversation.messages.length === 0}
                      size="small"
                      sx={{ mr: 1, color: '#f97316' }}
                    >
                      Clear
                    </Button>
                  </span>
                </Tooltip>
                <Tooltip title="Export conversation">
                  <span>
                    {/* Wrap in span to allow tooltip on disabled button */}
                    <Button
                      startIcon={<FileDownload />}
                      onClick={exportChat}
                      disabled={!selectedConversation || selectedConversation.messages.length === 0}
                      size="small"
                      sx={{ color: '#f97316' }}
                    >
                      Export
                    </Button>
                  </span>
                </Tooltip>
              </Box>
            </Box>
            <Typography variant="subtitle1" color="text.secondary" sx={{ mt: 1 }}>
              {selectedConversation
                ? `This conversation has ${selectedConversation.messages.length} messages. Ask questions about patient treatments, medications, and procedures.`
                : 'Select a conversation from the sidebar or create a new one to get started.'}
            </Typography>
          </Container>
        </Box>

        {/* Chat area */}
        <Box
          sx={{
            flexGrow: 1,
            overflowY: 'auto',
            p: 2,
          }}
        >
          <Container maxWidth="md">
            {(!selectedConversation || messages.length === 0) && !loading && !pendingQuestion && (
              <Box
                sx={{
                  textAlign: 'center',
                  py: 8,
                  display: 'flex',
                  flexDirection: 'column',
                  alignItems: 'center',
                  justifyContent: 'center',
                  height: '60vh',
                }}
              >
                <Box
                  sx={{
                    bgcolor: theme.palette.background.paper,
                    borderRadius: 2,
                    p: 4,
                    maxWidth: 500,
                    boxShadow: '0 4px 20px rgba(0,0,0,0.05)',
                    border: '1px solid',
                    borderColor: theme.palette.divider,
                  }}
                >
                  <Typography variant="h5" color="#f97316" gutterBottom>
                    {selectedConversation
                      ? 'Start a new conversation'
                      : 'Welcome to AI Dental Assistant'}
                  </Typography>

                  <Typography variant="body1" color="text.secondary">
                    {selectedConversation
                      ? 'Start by asking a question or uploading a patient document to analyze.'
                      : 'Select an existing conversation from the sidebar or create a new one to begin.'}
                  </Typography>

                  <Typography variant="subtitle2" color="text.secondary" sx={{ mt: 2, mb: 1 }}>
                    You can ask about:
                  </Typography>

                  <Box
                    sx={{
                      display: 'flex',
                      flexWrap: 'wrap',
                      justifyContent: 'center',
                      gap: 1,
                      mb: 3,
                    }}
                  >
                    <Chip label="Patient treatments" size="small" />
                    <Chip label="Medications" size="small" />
                    <Chip label="Procedures" size="small" />
                    <Chip label="Diagnoses" size="small" />
                    <Chip label="Treatment plans" size="small" />
                  </Box>

                  {!selectedConversation ? (
                    <Button
                      variant="contained"
                      startIcon={<Add />}
                      onClick={() => {
                        setDialogMode('create');
                        setDialogOpen(true);
                      }}
                      size="large"
                      fullWidth
                      style={{
                        backgroundColor: '#f97316',
                      }}
                    >
                      New Conversation
                    </Button>
                  ) : (
                    <Typography
                      variant="caption"
                      color="text.secondary"
                      sx={{ display: 'block', mt: 2 }}
                    >
                      Your conversation is ready. Type your question below to begin.
                    </Typography>
                  )}
                </Box>
              </Box>
            )}

            {messages.map(msg => (
              <Box key={msg.id} sx={{ mb: 4 }}>
                {/* User message */}
                {msg.role === 'user' && (
                  <Card variant="outlined" sx={{ mb: 2, borderRadius: '12px' }}>
                    <CardContent sx={{ pb: 1 }}>
                      <Box sx={{ display: 'flex', alignItems: 'flex-start', mb: 1 }}>
                        <Box sx={{ flexGrow: 1 }}>
                          <Typography component="div" sx={{ whiteSpace: 'pre-wrap' }}>
                            {msg.content}
                          </Typography>

                          {msg.file && (
                            <Chip
                              icon={<PictureAsPdf />}
                              label={msg.file.fileName}
                              size="small"
                              variant="outlined"
                              color="warning"
                              sx={{ mt: 1 }}
                            />
                          )}
                        </Box>
                        <Tooltip title="Copy to clipboard">
                          <IconButton
                            size="small"
                            onClick={() => copyToClipboard(msg.content, `q-${msg.id}`)}
                            sx={{ ml: 1 }}
                          >
                            {copiedMessageId === `q-${msg.id}` ? (
                              <Typography variant="caption">Copied!</Typography>
                            ) : (
                              <ContentCopy fontSize="small" />
                            )}
                          </IconButton>
                        </Tooltip>
                      </Box>
                      <Typography variant="caption" color="text.secondary">
                        {msg.timestamp?.toLocaleString()}
                      </Typography>
                    </CardContent>
                  </Card>
                )}

                {/* AI response */}
                {msg.role === 'assistant' && (
                  <Card
                    sx={{
                      ml: { xs: 2, sm: 4 },
                      mb: 2,
                      bgcolor: '#F97316',
                      color: theme.palette.primary.contrastText,
                      borderRadius: '12px',
                    }}
                  >
                    <CardContent sx={{ pb: 1 }}>
                      <Box sx={{ display: 'flex', alignItems: 'flex-start', mb: 1 }}>
                        <Box sx={{ flexGrow: 1 }}>
                          <Box
                            sx={{
                              '& p': { mb: 1 },
                              '& h1, & h2, & h3, & h4, & h5, & h6': { mt: 2, mb: 1 },
                              '& ul, & ol': { pl: 2 },
                              '& code': {
                                backgroundColor: 'rgba(255, 255, 255, 0.1)',
                                p: 0.5,
                                borderRadius: 1,
                                fontFamily: 'monospace',
                              },
                              '& sup a': {
                                textDecoration: 'none',
                                padding: '0 2px',
                                borderRadius: '3px',
                                fontWeight: 'bold',
                                marginLeft: '2px',
                                cursor: 'pointer',
                              },
                              '& [id^="eWzv"]': {
                                display: 'block',
                                margin: '10px 0',
                                padding: '5px 10px',
                                backgroundColor: 'rgba(255, 255, 255, 0.08)',
                                borderLeft: '3px solid rgba(255, 255, 255, 0.2)',
                                borderRadius: '3px',
                                fontSize: '0.9em',
                              },
                            }}
                          >
                            <ReactMarkdown
                              remarkPlugins={[remarkHeaderId, remarkGfm]}
                              components={{
                                a: ({ node, ...props }) => {
                                  // Special handling for footnote links
                                  if (props.href && props.href.startsWith('#') && props.children) {
                                    // Safely handle various types of children
                                    const childText = Array.isArray(props.children)
                                      ? String(props.children[0] || '')
                                      : String(props.children || '');
                                    if (childText.startsWith('^') && childText.endsWith('^')) {
                                      // Extract the number between the carets
                                      const footnoteNumber = childText.substring(
                                        1,
                                        childText.length - 1
                                      );
                                      return (
                                        <sup>
                                          <a
                                            {...props}
                                            style={{
                                              textDecoration: 'underline',
                                              color: 'blue',
                                              cursor: 'pointer',
                                              fontSize: '0.75em',
                                              opacity: '0.6',
                                            }}
                                          >
                                            {footnoteNumber}
                                          </a>
                                        </sup>
                                      );
                                    }
                                  }
                                  return (
                                    <a
                                      {...props}
                                      target="_blank"
                                      style={{
                                        color: 'blue',
                                        opacity: '0.8',
                                        textDecoration: 'underline',
                                        fontWeight: 500,
                                      }}
                                    />
                                  );
                                },
                                h6: ({ node, ...props }) => {
                                  // Special handling for footnote section headings
                                  if (props.id) {
                                    return (
                                      <h6
                                        {...props}
                                        style={{
                                          color: 'blue',
                                          opacity: 0.8,
                                        }}
                                      />
                                    );
                                  }
                                  return <h6 {...props} />;
                                },
                              }}
                            >
                              {msg.content}
                            </ReactMarkdown>
                          </Box>
                        </Box>
                        <Tooltip title="Copy to clipboard">
                          <IconButton
                            size="small"
                            onClick={() => copyToClipboard(msg.content, `a-${msg.id}`)}
                            sx={{ ml: 1, color: 'inherit', opacity: 0.7 }}
                          >
                            {copiedMessageId === `a-${msg.id}` ? (
                              <Typography variant="caption">Copied!</Typography>
                            ) : (
                              <ContentCopy fontSize="small" />
                            )}
                          </IconButton>
                        </Tooltip>
                      </Box>
                      <Divider sx={{ my: 1, borderColor: 'rgba(255,255,255,0.1)' }} />
                      <Typography variant="caption" sx={{ display: 'block', textAlign: 'right' }}>
                        AI Assistant • {msg.timestamp?.toLocaleString()}
                      </Typography>
                    </CardContent>
                  </Card>
                )}
              </Box>
            ))}

            {/* Pending question */}
            {pendingQuestion && (
              <Card variant="outlined" sx={{ mb: 3, borderRadius: '12px' }}>
                <CardContent>
                  <Typography component="div" sx={{ whiteSpace: 'pre-wrap' }}>
                    {pendingQuestion}
                  </Typography>
                  {selectedFile && (
                    <Chip
                      icon={<PictureAsPdf />}
                      label={selectedFile.name}
                      size="small"
                      variant="outlined"
                      color="warning"
                      sx={{ mt: 1 }}
                    />
                  )}
                </CardContent>
              </Card>
            )}

            {/* Loading indicator */}
            {loading && (
              <Box
                sx={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  my: 3,
                  px: 2,
                  py: 3,
                }}
              >
                <CircularProgress size={24} color="warning" />
                <Typography sx={{ ml: 2 }} color="warning">
                  Processing...
                </Typography>
              </Box>
            )}

            <div ref={bottomRef} />
          </Container>
        </Box>

        {/* Input area */}
        <Box
          component="footer"
          sx={{
            position: 'sticky',
            bottom: 0,
            p: 2,
            zIndex: 10,
            backgroundColor: 'transparent',
            borderTop: `1px solid ${theme.palette.divider}`,
          }}
        >
          <Container maxWidth="md">
            <Box sx={{ display: 'flex', alignItems: 'center' }}>
              <Tooltip title="Attach PDF document">
                <IconButton
                  component="label"
                  color={selectedFile ? 'warning' : 'default'}
                  sx={{ mr: 1 }}
                  disabled={!selectedConversation}
                >
                  <AttachFile />
                  <input hidden type="file" accept="application/pdf" onChange={handleFileSelect} />
                </IconButton>
              </Tooltip>

              {selectedFile && (
                <Chip
                  icon={<PictureAsPdf />}
                  label={selectedFile.name}
                  onDelete={() => setSelectedFile(undefined)}
                  variant="outlined"
                  color="warning"
                  sx={{ mr: 1, maxWidth: { xs: 150, sm: 200, md: 300 } }}
                />
              )}

              <TextareaAutosize
                ref={textAreaRef}
                minRows={1}
                maxRows={5}
                placeholder={
                  !selectedConversation
                    ? 'Select or create a conversation to begin'
                    : selectedFile
                      ? 'Ask about the uploaded PDF...'
                      : 'Enter patient details, clinical questions, or treatment concerns...'
                }
                value={input}
                onChange={e => setInput(e.target.value)}
                onKeyDown={handleKeyPress}
                disabled={!selectedConversation}
                style={{
                  flexGrow: 1,
                  marginRight: 8,
                  padding: 12,
                  borderRadius: 8,
                  border: `1px solid ${theme.palette.divider}`,
                  fontFamily: theme.typography.fontFamily,
                  fontSize: '1rem',
                  resize: 'none',
                }}
              />

              {input.trim() && (
                <Tooltip title="Clear input">
                  <IconButton size="small" onClick={() => setInput('')} sx={{ mr: 1 }}>
                    <Clear />
                  </IconButton>
                </Tooltip>
              )}

              <Tooltip
                title={selectedConversation ? 'Send message' : 'Select a conversation first'}
              >
                <span>
                  <Fab
                    size="medium"
                    onClick={handleSubmit}
                    disabled={loading || (!input.trim() && !selectedFile) || !selectedConversation}
                    sx={{
                      bgcolor: '#EA580C',
                      color: 'common.white',
                      '&:hover': {
                        bgcolor: '#DC4E0B',
                      },
                      display: 'flex',
                      alignItems: 'left',
                      justifyContent: 'center',
                    }}
                  >
                    <Send />
                  </Fab>
                </span>
              </Tooltip>
            </Box>
          </Container>
        </Box>

        <Dialog open={dialogOpen} onClose={handleDialogClose}>
          <DialogTitle>
            {dialogMode === 'create' ? 'Create New Conversation' : 'Rename Conversation'}
          </DialogTitle>
          <DialogContent>
            <DialogContentText>
              {dialogMode === 'create'
                ? 'Enter a title for your new conversation.'
                : 'Enter a new title for this conversation.'}
            </DialogContentText>
            <TextField
              autoFocus
              margin="dense"
              label="Title"
              type="text"
              fullWidth
              variant="outlined"
              value={dialogTitle}
              onChange={e => setDialogTitle(e.target.value)}
              color="warning"
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={handleDialogClose} color="warning">
              Cancel
            </Button>
            <Button onClick={handleDialogSubmit} color="warning">
              {dialogMode === 'create' ? 'Create' : 'Rename'}
            </Button>
          </DialogActions>
        </Dialog>

        <Snackbar
          open={alert.open}
          autoHideDuration={6000}
          onClose={handleAlertClose}
          anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
        >
          <Alert
            onClose={handleAlertClose}
            severity={alert.severity}
            variant="filled"
            sx={{ width: '100%' }}
          >
            {alert.message}
          </Alert>
        </Snackbar>
      </Box>
    </Box>
  );
};

export default Dashboard;
