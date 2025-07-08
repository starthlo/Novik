import useSWR from 'swr';
import { Conversation } from '../types';
import { patientService } from '../services';

export const useConverstaions = () => {
  const { data, isLoading, error, mutate } = useSWR<Conversation[]>('conversations', () =>
    patientService.getConversations()
  );
  return {
    conversations: data || [],
    isLoading,
    error,
    mutate,
  };
};
