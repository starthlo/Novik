import { useState, useEffect, useRef, ChangeEvent } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkHeaderId from 'remark-heading-id';
import { v4 as uuidv4 } from 'uuid';
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
} from '@mui/material';
import {
  AttachFile,
  Send,
  DeleteOutline,
  ContentCopy,
  FileDownload,
  PictureAsPdf,
  Clear,
} from '@mui/icons-material';

type ChatMessage = {
  id: string;
  question: string;
  answer: string;
  timestamp: Date;
  hasPdf?: boolean;
  pdfName?: string;
};

type AlertType = {
  open: boolean;
  message: string;
  severity: 'error' | 'warning' | 'info' | 'success';
};

const Dashboard = () => {
  const theme = useTheme();
  const [input, setInput] = useState('');
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [chatHistory, setChatHistory] = useState<ChatMessage[]>([]);
  const [loading, setLoading] = useState(false);
  const [pendingQuestion, setPendingQuestion] = useState<string | null>(null);
  const bottomRef = useRef<HTMLDivElement | null>(null);
  const textAreaRef = useRef<HTMLTextAreaElement | null>(null);
  const [sessionId, setSessionId] = useState('');
  const [alert, setAlert] = useState<AlertType>({
    open: false,
    message: '',
    severity: 'info',
  });
  const [copiedMessageId, setCopiedMessageId] = useState<string | null>(null);

  useEffect(() => {
    setSessionId(uuidv4());
  }, []);

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
    if (chatHistory.length === 0) return;

    if (confirm('Are you sure you want to clear the conversation history?')) {
      setChatHistory([]);
      setPendingQuestion(null);
      showAlert('Conversation cleared', 'success');
    }
  };

  const copyToClipboard = (text: string, id: string) => {
    navigator.clipboard.writeText(text).then(
      () => {
        setCopiedMessageId(id);
      },
      err => {
        console.error('Could not copy text: ', err);
        showAlert('Failed to copy to clipboard', 'error');
      }
    );
  };

  const exportChat = () => {
    if (chatHistory.length === 0) {
      showAlert('No conversation to export', 'warning');
      return;
    }

    const exportDate = new Date().toISOString().slice(0, 10);
    const fileName = `dental-ai-chat_${exportDate}.txt`;

    const content = chatHistory
      .map(chat => {
        return `Q: ${chat.question}\n${chat.hasPdf ? `[Attached PDF: ${chat.pdfName}]\n` : ''}${chat.timestamp.toLocaleString()}\n\nA: ${chat.answer}\n${chat.timestamp.toLocaleString()}\n\n${'='.repeat(40)}\n\n`;
      })
      .join('');

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

    const messageText = selectedFile ? `${input}` : input.trim();
    const timestamp = new Date();
    const messageId = uuidv4();
    const hasPdf = !!selectedFile;
    const pdfName = selectedFile?.name;

    setPendingQuestion(messageText);
    setLoading(true);
    setInput('');
    scrollToBottom();

    try {
      let endpoint = '/api/dashboard/';
      let body: FormData | string;
      let headers: Record<string, string> | undefined;

      if (selectedFile) {
        const form = new FormData();
        form.append('sessionId', sessionId);
        form.append('pdf', selectedFile);
        form.append('message', input);
        body = form;
        endpoint = '/api/dashboard/pdf/';
        headers = undefined;
      } else {
        body = JSON.stringify({ sessionId, message: messageText });
        headers = { 'Content-Type': 'application/json' };
      }

      const res = await fetch(endpoint, { method: 'POST', headers, body });

      if (!res.ok) {
        throw new Error(`Server responded with ${res.status}: ${res.statusText}`);
      }

      const data = await res.json();

      setChatHistory(h => [
        ...h,
        {
          id: messageId,
          question: messageText,
          answer: data.message,
          timestamp,
          hasPdf,
          pdfName,
        },
      ]);

      setSelectedFile(null);
      scrollToBottom();
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Network error. Please try again.';

      setChatHistory(h => [
        ...h,
        {
          id: messageId,
          question: messageText,
          answer: errorMessage,
          timestamp,
          hasPdf,
          pdfName,
        },
      ]);

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

  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        height: '100vh',
        maxHeight: 'calc(100vh - 64px)',
      }}
    >
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
            <Typography variant="h5" color="primary">
              AI Dental Assistant
            </Typography>
            <Box>
              <Tooltip title="Clear conversation">
                <Button
                  startIcon={<DeleteOutline />}
                  onClick={clearChat}
                  disabled={chatHistory.length === 0}
                  size="small"
                  sx={{ mr: 1 }}
                >
                  Clear
                </Button>
              </Tooltip>
              <Tooltip title="Export conversation">
                <Button
                  startIcon={<FileDownload />}
                  onClick={exportChat}
                  disabled={chatHistory.length === 0}
                  variant="outlined"
                  size="small"
                >
                  Export
                </Button>
              </Tooltip>
            </Box>
          </Box>
          <Typography variant="subtitle1" color="text.secondary" sx={{ mt: 1 }}>
            Ask questions about patient treatments, medications, and procedures
          </Typography>
        </Container>
      </Box>

      {/* Main chat area */}
      <Box
        component="main"
        sx={{
          flexGrow: 1,
          overflowY: 'auto',
          p: 2,
        }}
      >
        <Container maxWidth="md">
          {chatHistory.length === 0 && !loading && !pendingQuestion && (
            <Box
              sx={{
                textAlign: 'center',
                py: 8,
                opacity: 0.7,
              }}
            >
              <Typography variant="h6" color="text.secondary" gutterBottom>
                Your conversation with AI Dental Assistant
              </Typography>
              <Typography variant="body1" color="text.secondary">
                Start by asking a question or uploading a patient document
              </Typography>
            </Box>
          )}

          {chatHistory.map(chat => (
            <Box key={chat.id} sx={{ mb: 4 }}>
              {/* User message */}
              <Card variant="outlined" sx={{ mb: 2 }}>
                <CardContent sx={{ pb: 1 }}>
                  <Box sx={{ display: 'flex', alignItems: 'flex-start', mb: 1 }}>
                    <Box sx={{ flexGrow: 1 }}>
                      <Typography component="div" sx={{ whiteSpace: 'pre-wrap' }}>
                        {chat.question}
                      </Typography>

                      {chat.hasPdf && (
                        <Chip
                          icon={<PictureAsPdf />}
                          label={chat.pdfName}
                          size="small"
                          variant="outlined"
                          color="primary"
                          sx={{ mt: 1 }}
                        />
                      )}
                    </Box>
                    <Tooltip title="Copy to clipboard">
                      <IconButton
                        size="small"
                        onClick={() => copyToClipboard(chat.question, `q-${chat.id}`)}
                        sx={{ ml: 1 }}
                      >
                        {copiedMessageId === `q-${chat.id}` ? (
                          <Typography variant="caption">Copied!</Typography>
                        ) : (
                          <ContentCopy fontSize="small" />
                        )}
                      </IconButton>
                    </Tooltip>
                  </Box>
                  <Typography variant="caption" color="text.secondary">
                    {chat.timestamp.toLocaleString()}
                  </Typography>
                </CardContent>
              </Card>

              {/* AI response */}
              <Card
                sx={{
                  ml: { xs: 2, sm: 4 },
                  mb: 2,
                  bgcolor: theme.palette.primary.light,
                  color: theme.palette.primary.contrastText,
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
                        }}
                      >
                        <ReactMarkdown remarkPlugins={[remarkHeaderId]}>
                          {chat.answer}
                        </ReactMarkdown>
                      </Box>
                    </Box>
                    <Tooltip title="Copy to clipboard">
                      <IconButton
                        size="small"
                        onClick={() => copyToClipboard(chat.answer, `a-${chat.id}`)}
                        sx={{ ml: 1, color: 'inherit', opacity: 0.7 }}
                      >
                        {copiedMessageId === `a-${chat.id}` ? (
                          <Typography variant="caption">Copied!</Typography>
                        ) : (
                          <ContentCopy fontSize="small" />
                        )}
                      </IconButton>
                    </Tooltip>
                  </Box>
                  <Divider sx={{ my: 1, borderColor: 'rgba(255,255,255,0.1)' }} />
                  <Typography variant="caption" sx={{ display: 'block', textAlign: 'right' }}>
                    AI Assistant • {chat.timestamp.toLocaleString()}
                  </Typography>
                </CardContent>
              </Card>
            </Box>
          ))}

          {/* Pending question */}
          {pendingQuestion && (
            <Card variant="outlined" sx={{ mb: 3 }}>
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
                    color="primary"
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
              <CircularProgress size={24} />
              <Typography sx={{ ml: 2 }} color="text.secondary">
                Processing your request...
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
        }}
      >
        <Container maxWidth="md">
          <Box sx={{ display: 'flex', alignItems: 'flex-center' }}>
            <Tooltip title="Attach PDF document">
              <IconButton
                component="label"
                color={selectedFile ? 'primary' : 'default'}
                sx={{ mr: 1 }}
              >
                <AttachFile />
                <input hidden type="file" accept="application/pdf" onChange={handleFileSelect} />
              </IconButton>
            </Tooltip>

            {selectedFile && (
              <Chip
                icon={<PictureAsPdf />}
                label={selectedFile.name}
                onDelete={() => setSelectedFile(null)}
                variant="outlined"
                color="primary"
                sx={{ mr: 1, maxWidth: { xs: 150, sm: 200, md: 300 } }}
              />
            )}

            <TextareaAutosize
              ref={textAreaRef}
              minRows={1}
              maxRows={5}
              placeholder={
                selectedFile
                  ? 'Ask about the uploaded PDF...'
                  : 'Enter patient details, clinical questions, or treatment concerns...'
              }
              value={input}
              onChange={e => setInput(e.target.value)}
              onKeyDown={handleKeyPress}
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

            <Tooltip title="Send message">
              <Fab
                size="medium"
                onClick={handleSubmit}
                disabled={loading || (!input.trim() && !selectedFile)}
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
            </Tooltip>
          </Box>
        </Container>
      </Box>

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
  );
};

export default Dashboard;
