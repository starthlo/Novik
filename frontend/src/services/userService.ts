import apiClient from '../lib/apiClient';

export interface UserProfile {
  id: number;
  username: string;
  email: string;
  firstName: string;
  lastName: string;
  dob: string | null;
  phone: string | null;
  occupation: string | null;
  country: string | null;
  state: string | null;
  city: string | null;
}

export interface ChangePasswordData {
  current_password: string;
  new_password: string;
  confirm_password: string;
}

export const userService = {
  getProfile: async (): Promise<UserProfile> => {
    const res = await apiClient.get('user/profile/');
    return res.data;
  },

  updateProfile: async (data: Partial<UserProfile>): Promise<UserProfile> => {
    const res = await apiClient.patch('user/profile/update/', data);
    return res.data;
  },

  changePassword: async (data: ChangePasswordData): Promise<{ message: string }> => {
    const res = await apiClient.post('user/change-password/', data);
    return res.data;
  },
};
