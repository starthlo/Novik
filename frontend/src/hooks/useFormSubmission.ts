import { useState } from 'react';

interface UseFormSubmissionOptions {
  onSuccess?: () => void;
  onError?: (error: any) => void;
  successMessage?: string;
  errorMessage?: string;
  resetOnSuccess?: boolean;
}

interface FormSubmissionState {
  isSubmitting: boolean;
  showSuccess: boolean;
  showError: boolean;
  errorMessage: string;
  successMessage: string;
}

export const useFormSubmission = (options: UseFormSubmissionOptions = {}) => {
  const {
    onSuccess,
    onError,
    successMessage = 'Form submitted successfully!',
    errorMessage = 'An error occurred. Please try again.',
    resetOnSuccess = true,
  } = options;

  const [state, setState] = useState<FormSubmissionState>({
    isSubmitting: false,
    showSuccess: false,
    showError: false,
    errorMessage: '',
    successMessage: '',
  });

  const submitForm = async (submitFunction: () => Promise<any>, formReset?: () => void) => {
    setState(prev => ({ ...prev, isSubmitting: true, showError: false, showSuccess: false }));

    try {
      const result = await submitFunction();

      setState(prev => ({
        ...prev,
        isSubmitting: false,
        showSuccess: true,
        successMessage: successMessage,
      }));

      if (resetOnSuccess && formReset) {
        formReset();
      }

      if (onSuccess) {
        onSuccess();
      }

      return result;
    } catch (error) {
      console.error('Form submission error:', error);

      const errorMsg = error instanceof Error ? error.message : errorMessage;

      setState(prev => ({
        ...prev,
        isSubmitting: false,
        showError: true,
        errorMessage: errorMsg,
      }));

      if (onError) {
        onError(error);
      }

      throw error;
    }
  };

  const closeSuccess = () => {
    setState(prev => ({ ...prev, showSuccess: false }));
  };

  const closeError = () => {
    setState(prev => ({ ...prev, showError: false }));
  };

  const resetState = () => {
    setState({
      isSubmitting: false,
      showSuccess: false,
      showError: false,
      errorMessage: '',
      successMessage: '',
    });
  };

  return {
    ...state,
    submitForm,
    closeSuccess,
    closeError,
    resetState,
  };
};
