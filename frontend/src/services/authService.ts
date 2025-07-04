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
    const res = await apiClient.post('auth/google-login/', { token });
    return res.data;
  },
};
