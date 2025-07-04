import apiClient from '../lib/apiClient';
import { useAuthStore } from '../stores/auth';

export const patientService = {
  /**
   * Send a query to the Novik assistant (with or without PDF)
   * @param sessionId - The session ID for conversation context
   * @param message - The text message to send (optional if PDF is provided)
   * @param pdf - Optional PDF file to process
   * @param createdBy - Optional user ID of the creator
   * @returns The AI response
   */
  sendQuery: async (sessionId: string, message?: string, pdf?: File, createdBy?: number) => {
    if (pdf) {
      const formData = new FormData();
      formData.append('session_id', sessionId);
      formData.append('pdf', pdf);

      if (message) {
        formData.append('message', message);
      }

      if (createdBy !== undefined) {
        formData.append('createdBy', createdBy.toString());
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
        sessionId,
        message,
        createdBy,
      });
      return res.data;
    }
  },
};
