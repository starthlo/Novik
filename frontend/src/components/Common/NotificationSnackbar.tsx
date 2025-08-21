import { Snackbar, Alert, Slide, SlideProps } from '@mui/material';
import { novikTheme } from '../../styles/theme';

interface NotificationSnackbarProps {
  open: boolean;
  message: string;
  severity?: 'success' | 'error' | 'warning' | 'info';
  onClose: () => void;
  duration?: number;
  position?: {
    vertical: 'top' | 'bottom';
    horizontal: 'left' | 'center' | 'right';
  };
}

function SlideTransition(props: SlideProps) {
  return <Slide {...props} direction="down" />;
}

const NotificationSnackbar: React.FC<NotificationSnackbarProps> = ({
  open,
  message,
  severity = 'info',
  onClose,
  duration = 6000,
  position = { vertical: 'top', horizontal: 'center' },
}) => {
  const getBackgroundColor = () => {
    switch (severity) {
      case 'success':
        return novikTheme.colors.primary;
      case 'error':
        return novikTheme.colors.danger;
      case 'warning':
        return novikTheme.colors.secondary;
      case 'info':
      default:
        return novikTheme.colors.primary;
    }
  };

  return (
    <Snackbar
      open={open}
      autoHideDuration={duration}
      onClose={onClose}
      anchorOrigin={position}
      TransitionComponent={SlideTransition}
      sx={{
        '& .MuiSnackbar-root': {
          top: position.vertical === 'top' ? '80px' : undefined,
        },
      }}
    >
      <Alert
        onClose={onClose}
        severity={severity}
        variant="filled"
        elevation={6}
        sx={{
          backgroundColor: getBackgroundColor(),
          color: '#ffffff',
          fontFamily: novikTheme.typography.fontFamily,
          fontSize: '0.95rem',
          minWidth: '300px',
          boxShadow: novikTheme.shadows.large,
          '& .MuiAlert-icon': {
            color: '#ffffff',
          },
          '& .MuiAlert-action': {
            color: '#ffffff',
          },
          '& .MuiIconButton-root': {
            color: 'rgba(255, 255, 255, 0.8)',
            '&:hover': {
              color: '#ffffff',
              backgroundColor: 'rgba(255, 255, 255, 0.1)',
            },
          },
        }}
      >
        {message}
      </Alert>
    </Snackbar>
  );
};

export default NotificationSnackbar;
