import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Typography,
  styled,
  IconButton,
} from '@mui/material';
import CheckCircleOutlineIcon from '@mui/icons-material/CheckCircleOutline';
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import WarningAmberIcon from '@mui/icons-material/WarningAmber';
import CloseIcon from '@mui/icons-material/Close';
import { novikTheme } from '../../styles/theme';

interface ConfirmationDialogProps {
  open: boolean;
  onClose: () => void;
  title?: string;
  message: string;
  type?: 'success' | 'error' | 'warning' | 'info';
  confirmText?: string;
  cancelText?: string;
  onConfirm?: () => void;
  showCancel?: boolean;
}

const DialogContainer = styled(Dialog)({
  '& .MuiDialog-paper': {
    borderRadius: novikTheme.borderRadius.medium,
    padding: '8px',
    minWidth: '320px',
    maxWidth: '500px',
    fontFamily: novikTheme.typography.fontFamily,
  },
});

const StyledDialogTitle = styled(DialogTitle)({
  display: 'flex',
  alignItems: 'center',
  gap: '12px',
  paddingBottom: '8px',
  fontFamily: novikTheme.typography.fontFamily,
  '& .MuiTypography-root': {
    fontFamily: novikTheme.typography.fontFamily,
  },
});

const StyledDialogContent = styled(DialogContent)({
  paddingTop: '8px',
  fontFamily: novikTheme.typography.fontFamily,
});

const StyledDialogActions = styled(DialogActions)({
  padding: '16px 24px',
  gap: '8px',
});

const ActionButton = styled(Button)({
  borderRadius: novikTheme.borderRadius.small,
  padding: '8px 20px',
  fontFamily: novikTheme.typography.fontFamily,
  textTransform: 'none',
  fontWeight: 500,
});

const ConfirmButton = styled(ActionButton)({
  backgroundColor: novikTheme.colors.primary,
  color: '#ffffff',
  '&:hover': {
    backgroundColor: novikTheme.colors.primaryDark,
  },
});

const CancelButton = styled(ActionButton)({
  color: novikTheme.colors.textMuted,
  borderColor: novikTheme.colors.border,
  '&:hover': {
    backgroundColor: 'rgba(0, 0, 0, 0.04)',
    borderColor: novikTheme.colors.border,
  },
});

const CloseButton = styled(IconButton)({
  position: 'absolute',
  right: '8px',
  top: '8px',
  color: novikTheme.colors.textMuted,
  '&:hover': {
    backgroundColor: 'rgba(0, 0, 0, 0.04)',
  },
});

const getIcon = (type: string) => {
  const iconProps = { sx: { fontSize: '28px' } };

  switch (type) {
    case 'success':
      return (
        <CheckCircleOutlineIcon
          {...iconProps}
          sx={{ ...iconProps.sx, color: novikTheme.colors.primary }}
        />
      );
    case 'error':
      return (
        <ErrorOutlineIcon
          {...iconProps}
          sx={{ ...iconProps.sx, color: novikTheme.colors.danger }}
        />
      );
    case 'warning':
      return (
        <WarningAmberIcon
          {...iconProps}
          sx={{ ...iconProps.sx, color: novikTheme.colors.secondary }}
        />
      );
    case 'info':
    default:
      return (
        <InfoOutlinedIcon
          {...iconProps}
          sx={{ ...iconProps.sx, color: novikTheme.colors.primary }}
        />
      );
  }
};

const getTitle = (type: string, customTitle?: string) => {
  if (customTitle) return customTitle;

  switch (type) {
    case 'success':
      return 'Success!';
    case 'error':
      return 'Error';
    case 'warning':
      return 'Warning';
    case 'info':
    default:
      return 'Information';
  }
};

const ConfirmationDialog: React.FC<ConfirmationDialogProps> = ({
  open,
  onClose,
  title,
  message,
  type = 'info',
  confirmText = 'OK',
  cancelText = 'Cancel',
  onConfirm,
  showCancel = false,
}) => {
  const handleConfirm = () => {
    if (onConfirm) {
      onConfirm();
    }
    onClose();
  };

  return (
    <DialogContainer
      open={open}
      onClose={onClose}
      aria-labelledby="confirmation-dialog-title"
      aria-describedby="confirmation-dialog-description"
    >
      <CloseButton onClick={onClose} size="small">
        <CloseIcon fontSize="small" />
      </CloseButton>

      <StyledDialogTitle id="confirmation-dialog-title">
        {getIcon(type)}
        <Typography variant="h6" component="span" sx={{ fontWeight: 600 }}>
          {getTitle(type, title)}
        </Typography>
      </StyledDialogTitle>

      <StyledDialogContent>
        <Typography
          id="confirmation-dialog-description"
          sx={{
            color: novikTheme.colors.textMuted,
            lineHeight: 1.6,
          }}
        >
          {message}
        </Typography>
      </StyledDialogContent>

      <StyledDialogActions>
        {showCancel && (
          <CancelButton variant="outlined" onClick={onClose}>
            {cancelText}
          </CancelButton>
        )}
        <ConfirmButton variant="contained" onClick={handleConfirm}>
          {confirmText}
        </ConfirmButton>
      </StyledDialogActions>
    </DialogContainer>
  );
};

export default ConfirmationDialog;
