import apiClient from '../lib/apiClient';

export interface User {
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
  isActive: boolean;
  isStaff: boolean;
  isSuperuser: boolean;
  dateJoined: string;
  lastLogin: string | null;
  profileCompleted: boolean;
  conversationCount?: number;
}

export interface UserDetail extends User {
  recentConversations: {
    id: string;
    title: string;
    messageCount: number;
    updatedAt: string;
  }[];
  totalConversations: number;
}

export interface DashboardStats {
  userStats: {
    total: number;
    active: number;
    inactive: number;
    newThisMonth: number;
    newThisWeek: number;
  };
  conversationStats: {
    total: number;
    thisMonth: number;
    avgPerUser: number;
  };
  usersByCountry: {
    country: string;
    count: number;
  }[];
  recentRegistrations: {
    id: number;
    username: string;
    email: string;
    dateJoined: string;
  }[];
  topUsers: {
    id: number;
    username: string;
    email: string;
    conversationCount: number;
  }[];
}

export interface UsersResponse {
  users: User[];
  total: number;
  page: number;
  page_size: number;
  totalPages: number;
}

class AdminService {
  // Get all users with pagination and filtering
  async getUsers(params: {
    page?: number;
    page_size?: number;
    search?: string;
    isActive?: boolean;
    isStaff?: boolean;
    country?: string;
  }): Promise<UsersResponse> {
    const response = await apiClient.get('admin/users/', { params });
    return response.data;
  }

  async getUserDetail(userId: number): Promise<UserDetail> {
    const response = await apiClient.get(`admin/users/${userId}/`);
    return response.data;
  }

  async toggleUserStatus(userId: number): Promise<{ message: string; isActive: boolean }> {
    const response = await apiClient.post('admin/users/toggle-status/', { id: userId });
    return response.data;
  }

  async updateUser(userId: number, data: Partial<User>): Promise<{ message: string; user: User }> {
    const response = await apiClient.patch(`admin/users/${userId}/update/`, data);
    return response.data;
  }

  async deleteUser(userId: number): Promise<{ message: string }> {
    const response = await apiClient.delete(`admin/users/delete/`, { data: { id: userId } });
    return response.data;
  }

  // Make a user staff or revoke staff privileges
  async updateStaffStatus(
    userId: number,
    isStaff: boolean
  ): Promise<{ message: string; isStaff: boolean }> {
    const response = await apiClient.post(`admin/users/${userId}/staff/`, { isStaff: isStaff });
    return response.data;
  }

  async exportUsersCSV(): Promise<Blob> {
    const response = await apiClient.get('admin/users/export/', {
      responseType: 'blob',
    });
    return response.data;
  }

  async getDashboardStats(): Promise<DashboardStats> {
    const response = await apiClient.get('admin/dashboard/stats/');
    return response.data;
  }
}

export default new AdminService();
