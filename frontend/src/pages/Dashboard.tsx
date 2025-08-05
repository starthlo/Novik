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
  AppBar,
  Toolbar,
} from '@mui/material';
import {
  AttachFile,
  DeleteOutline,
  ContentCopy,
  FileDownload,
  PictureAsPdf,
  Clear,
  Add,
  Menu as MenuIcon,
  KeyboardArrowDown,
} from '@mui/icons-material';
import { useConverstaions } from '../hooks/useConversations';
import ArrowUpwardIcon from '@mui/icons-material/ArrowUpward';

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
  const [showScrollButton, setShowScrollButton] = useState(false);

  const { mutate } = useConverstaions();
  const [selectedConversation, setSelectedConversation] = useState<Conversation | null>(null);
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
    const chatContainer = document.querySelector('[data-chat-container]');
    if (chatContainer) {
      chatContainer.scrollTo({
        top: chatContainer.scrollHeight,
        behavior,
      });
    }
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

  const switchPatient = () => {
    if (
      confirm(
        'Are you sure you want to switch to a new patient? This will clear the current conversation.'
      )
    ) {
      if (selectedConversation) {
        handleClearConversation(selectedConversation.id);
      } else {
        // If no conversation is selected, just clear any pending state
        setInput('');
        setSelectedFile(undefined);
      }
      showAlert('Ready for new patient case', 'success');
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

    // Auto-create a conversation if none is selected
    if (!selectedConversation) {
      try {
        const newConversation = await patientService.createConversation('New Patient Case');
        mutate(data => [...(data || []), newConversation], false);
        setSelectedConversation(newConversation);

        // Now proceed with the message
        await submitMessage(newConversation, selectedFile ? `${input}` : input.trim());
      } catch (err) {
        showAlert('Failed to create conversation', 'error');
        return;
      }
    } else {
      await submitMessage(selectedConversation, selectedFile ? `${input}` : input.trim());
    }
  };

  const submitMessage = async (conversation: Conversation, messageText: string) => {
    setPendingQuestion(messageText);
    setLoading(true);
    setInput('');
    scrollToBottom();

    try {
      const response = await patientService.ask(messageText, selectedFile);

      // Update the UI immediately for better UX
      const updatedConversation: Conversation = {
        ...conversation,
        messages: [
          ...conversation.messages,
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
        convs?.map(conv => (conv.id === conversation.id ? updatedConversation : conv))
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

  // Auto-scroll when new messages are added
  useEffect(() => {
    if (messages.length > 0 || loading || pendingQuestion) {
      scrollToBottom();
    }
  }, [messages.length, loading, pendingQuestion]);

  // Handle scroll detection for floating button
  useEffect(() => {
    const chatContainer = document.querySelector('[data-chat-container]');
    if (!chatContainer) return;

    const handleScroll = () => {
      const { scrollTop, scrollHeight, clientHeight } = chatContainer;
      const isNearBottom = scrollHeight - scrollTop - clientHeight < 100;
      setShowScrollButton(!isNearBottom && messages.length > 0);
    };

    chatContainer.addEventListener('scroll', handleScroll);
    return () => chatContainer.removeEventListener('scroll', handleScroll);
  }, [messages.length]);

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
                <Tooltip title="Switch to new patient (clears current conversation)">
                  <Button
                    startIcon={<Add />}
                    onClick={switchPatient}
                    size="small"
                    variant="contained"
                    sx={{
                      mr: 1,
                      bgcolor: '#f97316',
                      '&:hover': { bgcolor: '#ea580c' },
                    }}
                  >
                    Switch Patient
                  </Button>
                </Tooltip>
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
                : 'Start typing below to begin a new patient case. Ask questions about treatments, medications, and procedures.'}
            </Typography>
          </Container>
        </Box>

        {/* Chat area */}
        <Box
          data-chat-container
          sx={{
            flexGrow: 1,
            overflowY: 'auto',
            p: 2,
            position: 'relative',
          }}
        >
          <Container maxWidth="md">
            {messages.length === 0 && !loading && !pendingQuestion && (
              <Box
                sx={{
                  display: 'flex',
                  flexDirection: 'column',
                  alignItems: 'center',
                  justifyContent: 'center',
                  minHeight: '60vh',
                  textAlign: 'center',
                  px: 2,
                }}
              >
                <Typography
                  variant="h4"
                  color="#f97316"
                  gutterBottom
                  sx={{ fontWeight: 600, mb: 3 }}
                >
                  Welcome to Novik!
                </Typography>

                <Typography
                  variant="body1"
                  color="text.secondary"
                  sx={{
                    maxWidth: 600,
                    mb: 6,
                    lineHeight: 1.6,
                    fontSize: '1.1rem',
                  }}
                >
                  Please enter your patient's case, including their age, weight, medications, and
                  medical history, or upload a PDF with the patient's medical history, making sure
                  to anonymize any personal data.
                </Typography>

                {/* ChatGPT-style input box */}
                <Box
                  sx={{
                    maxWidth: 700,
                    width: '100%',
                  }}
                >
                  <Box
                    sx={{
                      bgcolor: theme.palette.mode === 'dark' ? '#2f2f2f' : '#f7f7f8',
                      borderRadius: '24px',
                      border: `1px solid ${theme.palette.divider}`,
                      overflow: 'hidden',
                      '&:focus-within': {
                        borderColor: '#f97316',
                        boxShadow: '0 0 0 1px #f97316',
                      },
                    }}
                  >
                    {/* Input area */}
                    <TextareaAutosize
                      ref={textAreaRef}
                      minRows={1}
                      maxRows={8}
                      placeholder={
                        selectedFile ? 'Ask about the uploaded PDF...' : 'Ask anything...'
                      }
                      value={input}
                      onChange={e => setInput(e.target.value)}
                      onKeyDown={handleKeyPress}
                      style={{
                        width: '100%',
                        padding: '12px 16px 8px 16px',
                        border: 'none',
                        outline: 'none',
                        resize: 'none',
                        backgroundColor: 'transparent',
                        fontFamily: theme.typography.fontFamily,
                        fontSize: '16px',
                        lineHeight: '1.5',
                        color: theme.palette.text.primary,
                      }}
                    />

                    {/* Bottom toolbar */}
                    <Box
                      sx={{
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'space-between',
                        px: 2,
                        py: 1,
                      }}
                    >
                      {/* Left side buttons and file chip */}
                      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flex: 1, mr: 2 }}>
                        <Tooltip title="Attach PDF document">
                          <Box
                            component="label"
                            sx={{
                              display: 'flex',
                              alignItems: 'center',
                              gap: 0.5,
                              px: 1.5,
                              py: 0.5,
                              borderRadius: '16px',
                              border: `1px solid ${theme.palette.divider}`,
                              cursor: 'pointer',
                              color: selectedFile ? '#f97316' : theme.palette.text.secondary,
                              bgcolor: 'transparent',
                              '&:hover': {
                                bgcolor: 'rgba(249, 115, 22, 0.04)',
                                borderColor: '#f97316',
                                color: '#f97316',
                              },
                            }}
                          >
                            <AttachFile fontSize="small" />
                            <Typography
                              variant="body2"
                              sx={{
                                fontSize: '0.875rem',
                                fontWeight: 500,
                                color: 'inherit',
                              }}
                            >
                              Attach
                            </Typography>
                            <input
                              hidden
                              type="file"
                              accept="application/pdf"
                              onChange={handleFileSelect}
                            />
                          </Box>
                        </Tooltip>

                        {selectedFile && (
                          <Chip
                            icon={<PictureAsPdf />}
                            label={selectedFile.name}
                            onDelete={() => setSelectedFile(undefined)}
                            size="small"
                            variant="outlined"
                            color="warning"
                            sx={{
                              maxWidth: { xs: 120, sm: 200 },
                              height: 28,
                              '& .MuiChip-label': {
                                fontSize: '0.75rem',
                                px: 1,
                              },
                            }}
                          />
                        )}

                        {input.trim() && (
                          <Tooltip title="Clear input">
                            <IconButton
                              size="small"
                              onClick={() => setInput('')}
                              sx={{
                                color: theme.palette.text.secondary,
                                '&:hover': { bgcolor: 'rgba(0,0,0,0.04)' },
                              }}
                            >
                              <Clear fontSize="small" />
                            </IconButton>
                          </Tooltip>
                        )}
                      </Box>

                      {/* Right side send button */}
                      <Tooltip title="Send message">
                        <IconButton
                          onClick={handleSubmit}
                          disabled={loading || (!input.trim() && !selectedFile)}
                          size="small"
                          sx={{
                            bgcolor:
                              loading || (!input.trim() && !selectedFile)
                                ? theme.palette.mode === 'dark'
                                  ? '#404040'
                                  : '#e0e0e0'
                                : '#EA580C',
                            color:
                              loading || (!input.trim() && !selectedFile)
                                ? theme.palette.mode === 'dark'
                                  ? '#666666'
                                  : '#999999'
                                : 'white',
                            width: 32,
                            height: 32,
                            '&:hover': {
                              bgcolor:
                                loading || (!input.trim() && !selectedFile)
                                  ? theme.palette.mode === 'dark'
                                    ? '#404040'
                                    : '#e0e0e0'
                                  : '#DC4E0B',
                            },
                            '&:disabled': {
                              color: theme.palette.mode === 'dark' ? '#666666' : '#999999',
                            },
                          }}
                        >
                          <ArrowUpwardIcon fontSize="small" />
                        </IconButton>
                      </Tooltip>
                    </Box>
                  </Box>
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

          {/* Floating scroll to bottom button */}
          {showScrollButton && (
            <Fab
              size="small"
              onClick={() => scrollToBottom()}
              sx={{
                position: 'absolute',
                bottom: 16,
                right: 16,
                bgcolor: '#f97316',
                color: 'white',
                '&:hover': {
                  bgcolor: '#ea580c',
                },
                zIndex: 1000,
              }}
            >
              <KeyboardArrowDown />
            </Fab>
          )}
        </Box>

        {/* Input area - only show at bottom when there are messages */}
        {messages.length > 0 && (
          <Box
            component="footer"
            sx={{
              position: 'sticky',
              bottom: 0,
              p: 2,
              zIndex: 10,
              backgroundColor: 'transparent',
            }}
          >
            <Container maxWidth="md">
              <Box sx={{ maxWidth: 700, mx: 'auto' }}>
                <Box
                  sx={{
                    bgcolor: theme.palette.mode === 'dark' ? '#2f2f2f' : '#f7f7f8',
                    borderRadius: '24px',
                    border: `1px solid ${theme.palette.divider}`,
                    overflow: 'hidden',
                    '&:focus-within': {
                      borderColor: '#f97316',
                      boxShadow: '0 0 0 1px #f97316',
                    },
                  }}
                >
                  {/* Input area */}
                  <TextareaAutosize
                    ref={textAreaRef}
                    minRows={1}
                    maxRows={5}
                    placeholder={selectedFile ? 'Ask about the uploaded PDF...' : 'Ask anything...'}
                    value={input}
                    onChange={e => setInput(e.target.value)}
                    onKeyDown={handleKeyPress}
                    style={{
                      width: '100%',
                      padding: '12px 16px 8px 16px',
                      border: 'none',
                      outline: 'none',
                      resize: 'none',
                      backgroundColor: 'transparent',
                      fontFamily: theme.typography.fontFamily,
                      fontSize: '16px',
                      lineHeight: '1.5',
                      color: theme.palette.text.primary,
                    }}
                  />

                  {/* Bottom toolbar */}
                  <Box
                    sx={{
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'space-between',
                      px: 2,
                      py: 1,
                    }}
                  >
                    {/* Left side buttons and file chip */}
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flex: 1, mr: 2 }}>
                      <Tooltip title="Attach PDF document">
                        <Box
                          component="label"
                          sx={{
                            display: 'flex',
                            alignItems: 'center',
                            gap: 0.5,
                            px: 1.5,
                            py: 0.5,
                            borderRadius: '16px',
                            border: `1px solid ${theme.palette.divider}`,
                            cursor: 'pointer',
                            color: selectedFile ? '#f97316' : theme.palette.text.secondary,
                            bgcolor: 'transparent',
                            '&:hover': {
                              bgcolor: 'rgba(249, 115, 22, 0.04)',
                              borderColor: '#f97316',
                              color: '#f97316',
                            },
                          }}
                        >
                          <AttachFile fontSize="small" />
                          <Typography
                            variant="body2"
                            sx={{
                              fontSize: '0.875rem',
                              fontWeight: 500,
                              color: 'inherit',
                            }}
                          >
                            Attach
                          </Typography>
                          <input
                            hidden
                            type="file"
                            accept="application/pdf"
                            onChange={handleFileSelect}
                          />
                        </Box>
                      </Tooltip>

                      {selectedFile && (
                        <Chip
                          icon={<PictureAsPdf />}
                          label={selectedFile.name}
                          onDelete={() => setSelectedFile(undefined)}
                          size="small"
                          variant="outlined"
                          color="warning"
                          sx={{
                            maxWidth: { xs: 120, sm: 200 },
                            height: 28,
                            '& .MuiChip-label': {
                              fontSize: '0.75rem',
                              px: 1,
                            },
                          }}
                        />
                      )}

                      {input.trim() && (
                        <Tooltip title="Clear input">
                          <IconButton
                            size="small"
                            onClick={() => setInput('')}
                            sx={{
                              color: theme.palette.text.secondary,
                              '&:hover': { bgcolor: 'rgba(0,0,0,0.04)' },
                            }}
                          >
                            <Clear fontSize="small" />
                          </IconButton>
                        </Tooltip>
                      )}
                    </Box>

                    {/* Right side send button */}
                    <Tooltip title="Send message">
                      <IconButton
                        onClick={handleSubmit}
                        disabled={loading || (!input.trim() && !selectedFile)}
                        size="small"
                        sx={{
                          bgcolor:
                            loading || (!input.trim() && !selectedFile)
                              ? theme.palette.mode === 'dark'
                                ? '#404040'
                                : '#e0e0e0'
                              : '#EA580C',
                          color:
                            loading || (!input.trim() && !selectedFile)
                              ? theme.palette.mode === 'dark'
                                ? '#666666'
                                : '#999999'
                              : 'white',
                          width: 32,
                          height: 32,
                          '&:hover': {
                            bgcolor:
                              loading || (!input.trim() && !selectedFile)
                                ? theme.palette.mode === 'dark'
                                  ? '#404040'
                                  : '#e0e0e0'
                                : '#DC4E0B',
                          },
                          '&:disabled': {
                            color: theme.palette.mode === 'dark' ? '#666666' : '#999999',
                          },
                        }}
                      >
                        <ArrowUpwardIcon fontSize="small" />
                      </IconButton>
                    </Tooltip>
                  </Box>
                </Box>
              </Box>
            </Container>
          </Box>
        )}

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
