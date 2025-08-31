import apiClient from '../lib/apiClient';

export const authService = {
  login: async (email: string, password: string) => {
    const res = await apiClient.post('auth/login/', { email, password });
    return res.data;
  },
  register: async (payload: any) => {
    const res = await apiClient.post('auth/register/', payload);
    return res.data;
  },
  loginWithGoogle: async (token: string) => {
    const res = await apiClient.post('auth/google/token/', { token });
    return res.data;
  },
  forgotPassword: async (email: string) => {
    const res = await apiClient.post('auth/forgot-password/', { email });
    return res.data;
  },
  validateResetToken: async (uid: string, token: string) => {
    const res = await apiClient.post('auth/validate-reset-token/', { uid, token });
    return res.data;
  },
  resetPassword: async (uid: string, token: string, password: string, confirmPassword: string) => {
    const res = await apiClient.post('auth/reset-password/', {
      uid,
      token,
      password,
      confirm_password: confirmPassword,
    });
    return res.data;
  },
};
