import { useState, useEffect, useRef, ChangeEvent } from 'react';
import { v4 as uuidv4 } from 'uuid';
import { patientService } from '../services/patientService';
import type { Message, Conversation } from '../types';
import MarkdownContent from '../components/MarkdownContent';
import {
  Box,
  Container,
  Typography,
  IconButton,
  CircularProgress,
  Button,
  Alert,
  Snackbar,
  styled,
  TextareaAutosize,
  Chip,
} from '@mui/material';
import { AttachFile, Send, PictureAsPdf, Add, FileDownload } from '@mui/icons-material';
import { useConverstaions } from '../hooks/useConversations';
import { novikTheme } from '../styles/theme';
import NovikLogo from '../assets/novik-logo.png';

type AlertType = {
  open: boolean;
  message: string;
  severity: 'error' | 'warning' | 'info' | 'success';
};

const PageHero = styled(Box)({
  paddingTop: '22px',
  paddingBottom: '20px',
  textAlign: 'center',
  backgroundColor: '#ffffff',
  borderBottom: `1px solid ${novikTheme.colors.border}`,
  boxShadow: '0 2px 8px rgba(0,0,0,0.05)',
  position: 'fixed',
  top: '64px',
  left: 0,
  right: 0,
  zIndex: 100,
});

const PageTitle = styled(Box)({
  fontSize: '1.6rem',
  margin: '10px 0 0',
  fontWeight: 700,
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  fontFamily: novikTheme.typography.fontFamily,
});

const LogoImage = styled('img')({
  height: '120px',
  marginBottom: '0px',
});

const HeaderButtons = styled(Box)({
  position: 'fixed',
  bottom: 0,
  left: 0,
  width: '100%',
  backgroundColor: '#ffffff',
  padding: '10px',
  display: 'flex',
  justifyContent: 'center',
  gap: '12px',
  boxShadow: '0 -2px 6px rgba(0,0,0,0.1)',
  zIndex: 1000,
});

const AccentButton = styled(Button)({
  padding: '0.55rem 1rem',
  borderRadius: '8px',
  fontWeight: 600,
  fontSize: '0.85rem',
  textTransform: 'none',
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const GhostButton = styled(Button)({
  padding: '0.55rem 1rem',
  borderRadius: '8px',
  fontWeight: 600,
  fontSize: '0.85rem',
  textTransform: 'none',
  backgroundColor: '#ffffff',
  border: `1px solid ${novikTheme.colors.primary}`,
  color: novikTheme.colors.primary,
  fontFamily: novikTheme.typography.fontFamily,
  '&:hover': {
    backgroundColor: 'rgba(136, 169, 78, 0.05)',
  },
});

const WelcomeTitle = styled(Typography)({
  color: novikTheme.colors.primary,
  fontSize: '1.8rem',
  margin: '0.4rem 0 0.6rem',
  fontWeight: 600,
  fontFamily: novikTheme.typography.fontFamily,
});

const SubText = styled(Typography)({
  maxWidth: '820px',
  color: novikTheme.colors.textMuted,
  lineHeight: 1.6,
  margin: '0 auto 1.5rem',
  fontFamily: novikTheme.typography.fontFamily,
});

const InputContainer = styled(Box)({
  position: 'fixed',
  bottom: 64,
  left: 0,
  right: 0,
  padding: '16px',
  zIndex: 100,
  background: 'transparent',
});

const AskWrapper = styled(Box)({
  width: '100%',
  maxWidth: '820px',
  margin: '0 auto',
});

const AskBox = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  border: `1px solid ${novikTheme.colors.border}`,
  borderRadius: '24px',
  backgroundColor: '#ffffff',
  width: '100%',
  padding: '12px 16px 8px 16px',
  boxShadow: '0 2px 6px rgba(0,0,0,0.05)',
});

const InputToolbar = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'space-between',
  paddingTop: '8px',
});

const AttachButton = styled(Button)({
  backgroundColor: 'transparent',
  border: `1px solid ${novikTheme.colors.border}`,
  borderRadius: '16px',
  padding: '0.3rem 0.8rem',
  fontSize: '0.85rem',
  cursor: 'pointer',
  textTransform: 'none',
  color: novikTheme.colors.text,
  fontFamily: novikTheme.typography.fontFamily,
  minWidth: 'auto',
  transition: 'all 0.2s ease',
  '&:hover': {
    backgroundColor: 'rgba(136, 169, 78, 0.04)',
    borderColor: novikTheme.colors.primary,
    color: novikTheme.colors.primary,
  },
});

const SendButton = styled(IconButton)({
  backgroundColor: novikTheme.colors.primary,
  border: 0,
  borderRadius: '50%',
  width: '32px',
  height: '32px',
  color: '#ffffff',
  fontSize: '0.9rem',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  transition: 'all 0.2s ease',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
  '&:disabled': {
    backgroundColor: '#cccccc',
    color: '#999999',
  },
});

const MessageRow = styled(Box)<{ isUser?: boolean }>(({ isUser }) => ({
  display: 'flex',
  justifyContent: isUser ? 'flex-end' : 'flex-start',
  marginBottom: '10px',
}));

const MessageBubble = styled(Box)<{ isUser?: boolean }>(({ isUser }) => ({
  display: 'inline-block',
  textAlign: 'left',
  padding: '8px 12px',
  margin: '6px 0',
  borderRadius: '14px',
  lineHeight: 1.45,
  maxWidth: '78%',
  backgroundColor: isUser ? novikTheme.colors.primary : '#f4f4f4',
  color: isUser ? '#ffffff' : novikTheme.colors.text,
  border: isUser ? 'none' : `1px solid ${novikTheme.colors.border}`,
  fontFamily: novikTheme.typography.fontFamily,
  fontSize: '0.95rem',
  wordBreak: 'break-word',
  overflowX: 'auto',
  '& > *:first-child': {
    marginTop: 0,
  },
  '& > *:last-child': {
    marginBottom: 0,
  },
  // Better shadow for assistant messages
  ...(!isUser && {
    boxShadow: '0 1px 3px rgba(0,0,0,0.06)',
  }),
}));

const LoadingContainer = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '20px',
});

const Dashboard = () => {
  const [input, setInput] = useState('');
  const [selectedFile, setSelectedFile] = useState<File | undefined>(undefined);
  const [loading, setLoading] = useState(false);
  const messagesEndRef = useRef<HTMLDivElement | null>(null);
  const fileInputRef = useRef<HTMLInputElement | null>(null);
  const [alert, setAlert] = useState<AlertType>({
    open: false,
    message: '',
    severity: 'info',
  });

  const { mutate } = useConverstaions();
  const [selectedConversation, setSelectedConversation] = useState<Conversation | null>(null);
  const [messages, setMessages] = useState<Message[]>([]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages, loading]);

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

  const handleAttachClick = () => {
    fileInputRef.current?.click();
  };

  const clearFile = () => {
    setSelectedFile(undefined);
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  const newPatient = () => {
    if (messages.length > 0) {
      if (
        confirm(
          'Are you sure you want to start a new patient case? This will clear the current conversation.'
        )
      ) {
        setMessages([]);
        setInput('');
        clearFile();
        setSelectedConversation(null);
        showAlert('Ready for new patient case', 'success');
      }
    } else {
      showAlert('Ready for new patient case', 'info');
    }
  };

  const exportChat = () => {
    if (messages.length === 0) {
      showAlert('No conversation to export', 'warning');
      return;
    }

    const exportDate = new Date().toISOString().slice(0, 10);
    const fileName = `novik-conversation_${exportDate}.txt`;

    const content =
      `Novik Conversation - ${new Date().toLocaleString()}\n${'='.repeat(50)}\n\n` +
      messages
        .map(msg => {
          const prefix = msg.role === 'user' ? 'Patient Case: ' : 'Novik Assistant: ';
          return `${prefix}\n${msg.content}\n\n${'='.repeat(50)}\n`;
        })
        .join('\n');

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

    const messageText = selectedFile ? `${input} [PDF: ${selectedFile.name}]` : input.trim();

    // Add user message to UI
    const userMessage: Message = {
      id: uuidv4(),
      role: 'user',
      content: messageText,
      timestamp: new Date(),
      file: selectedFile ? { fileName: selectedFile.name, text: '' } : undefined,
    };

    setMessages(prev => [...prev, userMessage]);
    setLoading(true);
    setInput('');

    try {
      // Auto-create a conversation if none exists
      if (!selectedConversation) {
        const newConversation = await patientService.createConversation(
          'Patient Case ' + new Date().toLocaleDateString()
        );
        mutate(data => [...(data || []), newConversation], false);
        setSelectedConversation(newConversation);
      }

      const response = await patientService.ask(input || 'Please analyze this PDF', selectedFile);

      // Add AI response to UI
      const aiMessage: Message = {
        id: uuidv4(),
        role: 'assistant',
        content: response.message,
        timestamp: new Date(),
      };

      setMessages(prev => [...prev, aiMessage]);
      clearFile();
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Network error. Please try again.';
      showAlert('Error: ' + errorMessage, 'error');
      // Remove the pending user message on error
      setMessages(prev => prev.filter(msg => msg.id !== userMessage.id));
    } finally {
      setLoading(false);
    }
  };

  const handleKeyPress = (e: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit();
    }
  };

  const textAreaRef = useRef<HTMLTextAreaElement | null>(null);

  const hasMessages = messages.length > 0;

  return (
    <>
      {/* Page Hero Section - Only show when no messages */}
      {!hasMessages && (
        <PageHero>
          <PageTitle>
            <LogoImage src={NovikLogo} alt="Novik logo" />
            <Typography
              variant="h5"
              sx={{ fontWeight: 600, fontFamily: novikTheme.typography.fontFamily }}
            >
              AI Dental Assistant
            </Typography>
          </PageTitle>
        </PageHero>
      )}

      <HeaderButtons>
        <AccentButton startIcon={<Add />} onClick={newPatient}>
          New Patient
        </AccentButton>
        <GhostButton startIcon={<FileDownload />} onClick={exportChat}>
          Export
        </GhostButton>
      </HeaderButtons>

      {/* Main Content */}
      <Box
        sx={{
          backgroundColor: '#f7f7f8',
          minHeight: '100vh',
          paddingTop: hasMessages ? '64px' : '80px',
          paddingBottom: hasMessages ? '120px' : '80px',
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          justifyContent: !hasMessages ? 'center' : 'flex-start',
        }}
      >
        <Container maxWidth="lg">
          {!hasMessages && !loading ? (
            <Box sx={{ textAlign: 'center' }}>
              <WelcomeTitle variant="h2">Welcome to Novik</WelcomeTitle>
              <SubText>
                Please enter your patient's case, including their age, weight, medications, medical
                history, and the treatment to be performed, or upload a PDF with the patient's
                medical history, making sure to anonymize any personal data.
              </SubText>

              {/* Input Box for welcome state */}
              <AskWrapper sx={{ mt: 3 }}>
                <AskBox>
                  <input
                    ref={fileInputRef}
                    type="file"
                    accept="application/pdf"
                    onChange={handleFileSelect}
                    style={{ display: 'none' }}
                  />

                  <TextareaAutosize
                    ref={textAreaRef}
                    minRows={1}
                    maxRows={8}
                    placeholder={
                      selectedFile ? `Ask about ${selectedFile.name}...` : 'Ask anything...'
                    }
                    value={input}
                    onChange={e => setInput(e.target.value)}
                    onKeyDown={handleKeyPress}
                    disabled={loading}
                    style={{
                      width: '100%',
                      border: 'none',
                      outline: 'none',
                      resize: 'none',
                      backgroundColor: 'transparent',
                      fontFamily: novikTheme.typography.fontFamily,
                      fontSize: '1rem',
                      lineHeight: '1.5',
                      padding: '0',
                    }}
                  />

                  <InputToolbar>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <AttachButton
                        onClick={handleAttachClick}
                        startIcon={<AttachFile sx={{ fontSize: '0.9rem' }} />}
                        size="small"
                      >
                        Attach
                      </AttachButton>

                      {selectedFile && (
                        <Chip
                          icon={<PictureAsPdf />}
                          label={selectedFile.name}
                          onDelete={clearFile}
                          size="small"
                          variant="outlined"
                          sx={{
                            maxWidth: { xs: 150, sm: 250 },
                            '& .MuiChip-label': {
                              fontSize: '0.75rem',
                            },
                          }}
                        />
                      )}
                    </Box>

                    <SendButton
                      onClick={handleSubmit}
                      disabled={loading || (!input.trim() && !selectedFile)}
                      size="small"
                    >
                      <Send sx={{ fontSize: '18px' }} />
                    </SendButton>
                  </InputToolbar>
                </AskBox>
              </AskWrapper>
            </Box>
          ) : (
            <Box sx={{ py: 2 }} className="thin-scrollbar">
              {messages.map(msg => (
                <MessageRow key={msg.id} isUser={msg.role === 'user'}>
                  <MessageBubble isUser={msg.role === 'user'}>
                    {msg.role === 'assistant' ? (
                      <MarkdownContent content={msg.content} />
                    ) : (
                      <>
                        {msg.content}
                        {msg.file && (
                          <Chip
                            icon={<PictureAsPdf />}
                            label={msg.file.fileName}
                            size="small"
                            sx={{ ml: 1, backgroundColor: 'rgba(255,255,255,0.2)', color: '#fff' }}
                          />
                        )}
                      </>
                    )}
                  </MessageBubble>
                </MessageRow>
              ))}

              {loading && (
                <MessageRow isUser={false}>
                  <MessageBubble isUser={false}>
                    <LoadingContainer>
                      <CircularProgress
                        size={20}
                        sx={{ color: novikTheme.colors.primary, mr: 1 }}
                      />
                      <Typography sx={{ color: novikTheme.colors.textMuted }}>
                        Working on it...
                      </Typography>
                    </LoadingContainer>
                  </MessageBubble>
                </MessageRow>
              )}

              <div ref={messagesEndRef} />
            </Box>
          )}
        </Container>
      </Box>

      {hasMessages && (
        <InputContainer>
          <AskWrapper>
            <AskBox>
              <input
                ref={fileInputRef}
                type="file"
                accept="application/pdf"
                onChange={handleFileSelect}
                style={{ display: 'none' }}
              />

              <TextareaAutosize
                ref={textAreaRef}
                minRows={1}
                maxRows={8}
                placeholder={selectedFile ? `Ask about ${selectedFile.name}...` : 'Ask anything...'}
                value={input}
                onChange={e => setInput(e.target.value)}
                onKeyDown={handleKeyPress}
                disabled={loading}
                style={{
                  width: '100%',
                  border: 'none',
                  outline: 'none',
                  resize: 'none',
                  backgroundColor: 'transparent',
                  fontFamily: novikTheme.typography.fontFamily,
                  fontSize: '1rem',
                  lineHeight: '1.5',
                  padding: '0',
                }}
              />

              <InputToolbar>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                  <AttachButton
                    onClick={handleAttachClick}
                    startIcon={<AttachFile sx={{ fontSize: '0.9rem' }} />}
                    size="small"
                  >
                    Attach
                  </AttachButton>

                  {selectedFile && (
                    <Chip
                      icon={<PictureAsPdf />}
                      label={selectedFile.name}
                      onDelete={clearFile}
                      size="small"
                      variant="outlined"
                      sx={{
                        maxWidth: { xs: 150, sm: 250 },
                        '& .MuiChip-label': {
                          fontSize: '0.75rem',
                        },
                      }}
                    />
                  )}
                </Box>

                <SendButton
                  onClick={handleSubmit}
                  disabled={loading || (!input.trim() && !selectedFile)}
                  size="small"
                >
                  <Send sx={{ fontSize: '18px' }} />
                </SendButton>
              </InputToolbar>
            </AskBox>
          </AskWrapper>
        </InputContainer>
      )}

      <Snackbar
        open={alert.open}
        autoHideDuration={4000}
        onClose={handleAlertClose}
        anchorOrigin={{ vertical: 'top', horizontal: 'center' }}
      >
        <Alert
          onClose={handleAlertClose}
          severity={alert.severity}
          variant="filled"
          sx={{
            width: '100%',
            ...(alert.severity === 'success' && {
              backgroundColor: novikTheme.colors.primary,
              color: '#ffffff',
              '& .MuiAlert-icon': {
                color: '#ffffff',
              },
            }),
          }}
        >
          {alert.message}
        </Alert>
      </Snackbar>
    </>
  );
};

export default Dashboard;
