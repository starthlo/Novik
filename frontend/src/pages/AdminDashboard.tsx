import { useState, useEffect } from 'react';
import {
  FaUsers,
  FaUserCheck,
  FaComments,
  FaGlobe,
  FaCalendarAlt,
  FaStar,
} from 'react-icons/fa';
import { useNavigate } from 'react-router-dom';
import adminService, { DashboardStats } from '../services/adminService';

export default function AdminDashboard() {
  const [stats, setStats] = useState<DashboardStats | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    fetchDashboardStats();
  }, []);

  const fetchDashboardStats = async () => {
    try {
      setLoading(true);
      const data = await adminService.getDashboardStats();
      setStats(data);
      setError(null);
    } catch (err: any) {
      setError(err.response?.data?.error || 'Failed to load dashboard statistics');
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-green-500"></div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div
          className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded relative"
          role="alert"
        >
          <strong className="font-bold">Error!</strong>
          <span className="block sm:inline"> {error}</span>
        </div>
      </div>
    );
  }

  if (!stats) return null;

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white shadow-sm border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between items-center py-6">
            <div>
              <h1 className="text-3xl font-bold text-gray-900">Admin Dashboard</h1>
              <p className="mt-1 text-sm text-gray-600">Monitor and manage your application</p>
            </div>
            <div className="flex space-x-3">
              <button
                onClick={() => navigate('/admin/users')}
                className="px-4 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition duration-200 flex items-center"
              >
                <FaUsers className="mr-2" />
                Manage Users
              </button>
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Stats Grid */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6 mb-8">
          {/* Total Users */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm font-medium text-gray-600">Total Users</p>
                <p className="text-3xl font-bold text-gray-900 mt-2">
                  {stats.userStats.total.toLocaleString()}
                </p>
                <p className="text-xs text-gray-500 mt-2">
                  +{stats.userStats.newThisMonth} this month
                </p>
              </div>
              <div className="bg-blue-100 rounded-full p-3">
                <FaUsers className="h-6 w-6 text-blue-600" />
              </div>
            </div>
          </div>

          {/* Active Users */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm font-medium text-gray-600">Active Users</p>
                <p className="text-3xl font-bold text-green-600 mt-2">
                  {stats.userStats.active.toLocaleString()}
                </p>
                <p className="text-xs text-gray-500 mt-2">
                  {((stats.userStats.active / stats.userStats.total) * 100).toFixed(1)}% of total
                </p>
              </div>
              <div className="bg-green-100 rounded-full p-3">
                <FaUserCheck className="h-6 w-6 text-green-600" />
              </div>
            </div>
          </div>

          {/* Total Conversations */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm font-medium text-gray-600">Total Conversations</p>
                <p className="text-3xl font-bold text-gray-900 mt-2">
                  {stats.conversationStats.total.toLocaleString()}
                </p>
                <p className="text-xs text-gray-500 mt-2">
                  {stats.conversationStats.avgPerUser.toFixed(1)} per user
                </p>
              </div>
              <div className="bg-purple-100 rounded-full p-3">
                <FaComments className="h-6 w-6 text-purple-600" />
              </div>
            </div>
          </div>

          {/* New This Week */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm font-medium text-gray-600">New This Week</p>
                <p className="text-3xl font-bold text-gray-900 mt-2">
                  {stats.userStats.newThisWeek}
                </p>
                <p className="text-xs text-gray-500 mt-2">User registrations</p>
              </div>
              <div className="bg-orange-100 rounded-full p-3">
                <FaCalendarAlt className="h-6 w-6 text-orange-600" />
              </div>
            </div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Users by Country */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-lg font-semibold text-gray-900">Users by Country</h2>
              <FaGlobe className="text-gray-400" />
            </div>
            <div className="space-y-3">
              {stats.usersByCountry.slice(0, 5).map((country, index) => (
                <div key={index} className="flex items-center justify-between">
                  <div className="flex items-center">
                    <span className="text-sm font-medium text-gray-900">{country.country}</span>
                  </div>
                  <div className="flex items-center">
                    <span className="text-sm text-gray-600 mr-2">{country.count}</span>
                    <div className="w-24 bg-gray-200 rounded-full h-2">
                      <div
                        className="bg-green-600 h-2 rounded-full"
                        style={{
                          width: `${(country.count / stats.userStats.total) * 100}%`,
                        }}
                      ></div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>

          {/* Recent Registrations */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-lg font-semibold text-gray-900">Recent Registrations</h2>
              <FaUserCheck className="text-gray-400" />
            </div>
            <div className="space-y-3">
              {stats.recentRegistrations.map(user => (
                <div key={user.id} className="flex items-center justify-between">
                  <div>
                    <p className="text-sm font-medium text-gray-900">{user.username}</p>
                    <p className="text-xs text-gray-500">{user.email}</p>
                  </div>
                  <span className="text-xs text-gray-500">
                    {new Date(user.dateJoined).toLocaleDateString()}
                  </span>
                </div>
              ))}
            </div>
          </div>

          {/* Top Users */}
          <div className="bg-white rounded-lg shadow p-6">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-lg font-semibold text-gray-900">Most Active Users</h2>
              <FaStar className="text-gray-400" />
            </div>
            <div className="space-y-3">
              {stats.topUsers.map((user, index) => (
                <div key={user.id} className="flex items-center justify-between">
                  <div className="flex items-center">
                    <span className="text-lg font-bold text-gray-400 mr-3">#{index + 1}</span>
                    <div>
                      <p className="text-sm font-medium text-gray-900">{user.username}</p>
                      <p className="text-xs text-gray-500">{user.email}</p>
                    </div>
                  </div>
                  <span className="text-sm font-semibold text-green-600">
                    {user.conversationCount} chats
                  </span>
                </div>
              ))}
            </div>
          </div>
        </div>

        {/* Quick Stats */}
        <div className="mt-8 bg-green-600 rounded-lg shadow-lg p-6 text-white">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-6">
            <div>
              <p className="text-green-100 text-sm">Inactive Users</p>
              <p className="text-2xl font-bold">{stats.userStats.inactive}</p>
            </div>
            <div>
              <p className="text-green-100 text-sm">Conversations This Month</p>
              <p className="text-2xl font-bold">{stats.conversationStats.thisMonth}</p>
            </div>
            <div>
              <p className="text-green-100 text-sm">Countries Reached</p>
              <p className="text-2xl font-bold">{stats.usersByCountry.length}</p>
            </div>
            <div>
              <p className="text-green-100 text-sm">Avg Conversations/User</p>
              <p className="text-2xl font-bold">{stats.conversationStats.avgPerUser.toFixed(1)}</p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
