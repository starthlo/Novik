import axios, { AxiosError, AxiosInstance } from 'axios';
import { useAuthStore } from '../stores/auth';
import camelcaseKeys from 'camelcase-keys';
import snakecaseKeys from 'snakecase-keys';

const apiClient: AxiosInstance = axios.create({
  baseURL: '/api',
  headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
});

apiClient.interceptors.request.use(
  config => {
    const accessToken = useAuthStore.getState().accessToken;
    if (accessToken && config.headers) {
      config.headers.Authorization = `Bearer ${accessToken}`;
    }

    if (config.data) {
      config.data = snakecaseKeys(config.data, { deep: true });
    }

    return config;
  },
  error => Promise.reject(error)
);

apiClient.interceptors.response.use(
  response => {
    response.data = camelcaseKeys(response.data, { deep: true });
    return response;
  },
  (error: AxiosError) => {
    if (error.response?.status === 401) {
      useAuthStore.getState().logout();

      const currentPath = window.location.pathname;
      const isOnAuthPage = currentPath === '/login' || currentPath === '/register';

      if (!isOnAuthPage) {
        window.location.href = '/login';
      }
    }
    return Promise.reject(error);
  }
);

export default apiClient;
