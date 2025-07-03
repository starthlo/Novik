import { create } from 'zustand';
import { devtools, persist } from 'zustand/middleware';

interface AuthState {
  isAuthorized: boolean;
  accessToken: string | null;
  user: Record<string, any> | null;
  logout: () => void;
}

export const useAuthStore = create<AuthState>()(
  devtools(
    persist(
      set => ({
        isAuthorized: false,
        accessToken: null,
        user: null,
        logout: () => set({ isAuthorized: false, accessToken: null, user: null }),
      }),
      { name: 'authStore' }
    )
  )
);
