import apiClient from '../lib/apiClient';
import { useAuthStore } from '../stores/auth';
import { Conversation } from '../types';

export const patientService = {
  getConversations: async (): Promise<Conversation[]> => {
    const response = await apiClient.get('conversations/');
    return response.data;
  },

  getConversation: async (id: string): Promise<Conversation> => {
    const response = await apiClient.get(`conversations/${id}/`);
    return response.data;
  },

  createConversation: async (title: string = 'New Conversation'): Promise<Conversation> => {
    const response = await apiClient.post('conversations/', { title });
    return response.data;
  },

  updateConversationTitle: async (id: string, title: string): Promise<Conversation> => {
    const response = await apiClient.patch(`conversations/${id}/`, { title });
    return response.data;
  },

  deleteConversation: async (id: string): Promise<void> => {
    await apiClient.delete(`conversations/${id}/`);
  },

  clearConversation: async (id: string): Promise<Conversation> => {
    const response = await apiClient.post(`conversations/${id}/clear/`);
    return response.data;
  },

  ask: async (message?: string, pdf?: File) => {
    if (pdf) {
      const formData = new FormData();
      formData.append('pdf', pdf);

      if (message) {
        formData.append('message', message);
      }

      const accessToken = useAuthStore.getState().accessToken;
      const response = await fetch('/api/patient/assistant/', {
        method: 'POST',
        body: formData,
        headers: {
          Authorization: `Bearer ${accessToken}`,
        },
      });

      if (!response.ok) {
        throw new Error(`Server responded with ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();
      return {
        message: data.message,
      };
    } else {
      const res = await apiClient.post('patient/assistant/', {
        message,
      });
      return res.data;
    }
  },
};
