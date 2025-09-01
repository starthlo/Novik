import { useState, useEffect } from 'react';
import {
  FaEdit,
  FaTrash,
  FaBan,
  FaCheckCircle,
  FaSearch,
  FaUserShield,
  FaTimes,
  FaSave,
} from 'react-icons/fa';
import { useNavigate } from 'react-router-dom';
import adminService, { User, UserDetail } from '../services/adminService';

export default function AdminUsers() {
  const [users, setUsers] = useState<User[]>([]);
  const [selectedUser, setSelectedUser] = useState<UserDetail | null>(null);
  const [loading, setLoading] = useState(true);
  const [searchTerm, setSearchTerm] = useState('');
  const [filters, setFilters] = useState({
    isActive: '',
    isStaff: '',
    country: '',
  });
  const [currentPage, setCurrentPage] = useState(1);
  const [totalPages, setTotalPages] = useState(1);
  const [totalUsers, setTotalUsers] = useState(0);
  const [showEditModal, setShowEditModal] = useState(false);
  const [editingUser, setEditingUser] = useState<Partial<User> | null>(null);
  const [showUserDetail, setShowUserDetail] = useState(false);
  const navigate = useNavigate();

  useEffect(() => {
    fetchUsers();
  }, [currentPage, searchTerm, filters]);

  const fetchUsers = async () => {
    try {
      setLoading(true);
      const params: any = {
        page: currentPage,
        page_size: 20,
      };

      if (searchTerm) params.search = searchTerm;
      if (filters.isActive) params.isActive = filters.isActive === 'true';
      if (filters.isStaff) params.isStaff = filters.isStaff === 'true';
      if (filters.country) params.country = filters.country;

      const data = await adminService.getUsers(params);
      setUsers(data.users);
      setTotalPages(data.totalPages);
      setTotalUsers(data.total);
    } catch (err: any) {
      console.log(err.response?.data?.error || 'Failed to load users');
    } finally {
      setLoading(false);
    }
  };

  const fetchUserDetail = async (userId: number) => {
    try {
      const data = await adminService.getUserDetail(userId);
      setSelectedUser(data);
      setShowUserDetail(true);
    } catch (err) {
      console.error('Failed to fetch user details:', err);
    }
  };

  const handleToggleStatus = async (userId: number) => {
    try {
      const result = await adminService.toggleUserStatus(userId);
      setUsers(
        users.map(user => (user.id === userId ? { ...user, isActive: result.isActive } : user))
      );
    } catch (err) {
      console.error('Failed to toggle user status:', err);
    }
  };

  const handleDeleteUser = async (userId: number) => {
    if (
      !window.confirm('Are you sure you want to delete this user? This action cannot be undone.')
    ) {
      return;
    }

    try {
      await adminService.deleteUser(userId);
      setUsers(users.filter(user => user.id !== userId));
      setTotalUsers(totalUsers - 1);
    } catch (err: any) {
      alert(err.response?.data?.error || 'Failed to delete user');
    }
  };

  const handleToggleStaff = async (userId: number, currentStatus: boolean) => {
    try {
      const result = await adminService.updateStaffStatus(userId, !currentStatus);
      setUsers(
        users.map(user => (user.id === userId ? { ...user, isStaff: result.isStaff } : user))
      );
    } catch (err) {
      console.error('Failed to update staff status:', err);
    }
  };

  const handleEditUser = (user: User) => {
    setEditingUser({
      id: user.id,
      email: user.email,
      firstName: user.firstName,
      lastName: user.lastName,
      phone: user.phone,
      occupation: user.occupation,
      country: user.country,
      state: user.state,
      city: user.city,
    });
    setShowEditModal(true);
  };

  const handleSaveEdit = async () => {
    if (!editingUser || !editingUser.id) return;

    try {
      const result = await adminService.updateUser(editingUser.id, editingUser);
      setUsers(users.map(user => (user.id === editingUser.id ? result.user : user)));
      setShowEditModal(false);
      setEditingUser(null);
    } catch (err) {
      console.error('Failed to update user:', err);
    }
  };

  if (loading && users.length === 0) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-green-500"></div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white shadow-sm border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between items-center py-6">
            <div>
              <h1 className="text-3xl font-bold text-gray-900">User Management</h1>
              <p className="mt-1 text-sm text-gray-600">Manage {totalUsers} registered users</p>
            </div>
            <div className="flex space-x-3">
              <button
                onClick={() => navigate('/admin/dashboard')}
                className="px-4 py-2 border border-gray-300 text-gray-700 rounded-lg hover:bg-gray-50 transition duration-200"
              >
                Back to Dashboard
              </button>
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Search and Filters */}
        <div className="bg-white rounded-lg shadow p-6 mb-6">
          <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
            <div className="md:col-span-1">
              <div className="relative">
                <FaSearch className="absolute left-3 top-1/2 transform -translate-y-1/2 text-gray-400" />
                <input
                  type="text"
                  placeholder="Search users..."
                  value={searchTerm}
                  onChange={e => setSearchTerm(e.target.value)}
                  className="w-full pl-10 pr-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-green-500 focus:border-transparent"
                />
              </div>
            </div>
            <div>
              <select
                value={filters.isActive}
                onChange={e => setFilters({ ...filters, isActive: e.target.value })}
                className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-green-500 focus:border-transparent"
              >
                <option value="">All Status</option>
                <option value="true">Active</option>
                <option value="false">Inactive</option>
              </select>
            </div>
            <div>
              <select
                value={filters.isStaff}
                onChange={e => setFilters({ ...filters, isStaff: e.target.value })}
                className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-green-500 focus:border-transparent"
              >
                <option value="">All Users</option>
                <option value="true">Staff Only</option>
                <option value="false">Regular Users</option>
              </select>
            </div>
            <div>
              <input
                type="text"
                placeholder="Filter by country..."
                value={filters.country}
                onChange={e => setFilters({ ...filters, country: e.target.value })}
                className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-green-500 focus:border-transparent"
              />
            </div>
          </div>
        </div>

        {/* Users Table */}
        <div className="bg-white rounded-lg shadow overflow-hidden">
          <table className="min-w-full divide-y divide-gray-200">
            <thead className="bg-gray-50">
              <tr>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  User
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  Contact
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  Location
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  Status
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  Joined
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                  Actions
                </th>
              </tr>
            </thead>
            <tbody className="bg-white divide-y divide-gray-200">
              {users.map(user => (
                <tr key={user.id} className="hover:bg-gray-50">
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div>
                      <div className="text-sm font-medium text-gray-900">{user.username}</div>
                      <div className="text-sm text-gray-500">
                        {user.firstName} {user.lastName}
                      </div>
                      {user.isSuperuser && (
                        <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-purple-100 text-purple-800 mt-1">
                          Super Admin
                        </span>
                      )}
                      {user.isStaff && !user.isSuperuser && (
                        <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium bg-blue-100 text-blue-800 mt-1">
                          Staff
                        </span>
                      )}
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="text-sm text-gray-900">{user.email}</div>
                    <div className="text-sm text-gray-500">{user.phone || '-'}</div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    <div className="text-sm text-gray-900">{user.country || '-'}</div>
                    <div className="text-sm text-gray-500">
                      {user.city && user.state ? `${user.city}, ${user.state}` : '-'}
                    </div>
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap">
                    {user.isActive ? (
                      <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-green-100 text-green-800">
                        <FaCheckCircle className="mr-1" />
                        Active
                      </span>
                    ) : (
                      <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-red-100 text-red-800">
                        <FaBan className="mr-1" />
                        Inactive
                      </span>
                    )}
                    {user.conversationCount !== undefined && (
                      <div className="text-xs text-gray-500 mt-1">
                        {user.conversationCount} chats
                      </div>
                    )}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500">
                    {new Date(user.dateJoined).toLocaleDateString()}
                  </td>
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-medium">
                    <div className="flex space-x-2">
                      <button
                        onClick={() => fetchUserDetail(user.id)}
                        className="text-indigo-600 hover:text-indigo-900"
                        title="View Details"
                      >
                        👁️
                      </button>
                      <button
                        onClick={() => handleEditUser(user)}
                        className="text-blue-600 hover:text-blue-900"
                        title="Edit User"
                      >
                        <FaEdit />
                      </button>
                      {!user.isSuperuser && (
                        <>
                          <button
                            onClick={() => handleToggleStaff(user.id, user.isStaff)}
                            className={
                              user.isStaff
                                ? 'text-orange-600 hover:text-orange-900'
                                : 'text-purple-600 hover:text-purple-900'
                            }
                            title={user.isStaff ? 'Remove Staff' : 'Make Staff'}
                          >
                            <FaUserShield />
                          </button>
                          <button
                            onClick={() => handleToggleStatus(user.id)}
                            className={
                              user.isActive
                                ? 'text-yellow-600 hover:text-yellow-900'
                                : 'text-green-600 hover:text-green-900'
                            }
                            title={user.isActive ? 'Deactivate' : 'Activate'}
                          >
                            {user.isActive ? <FaBan /> : <FaCheckCircle />}
                          </button>
                          <button
                            onClick={() => handleDeleteUser(user.id)}
                            className="text-red-600 hover:text-red-900"
                            title="Delete User"
                          >
                            <FaTrash />
                          </button>
                        </>
                      )}
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        {/* Pagination */}
        <div className="flex items-center justify-between mt-6">
          <div className="text-sm text-gray-700">
            Showing page {currentPage} of {totalPages} ({totalUsers} total users)
          </div>
          <div className="flex space-x-2">
            <button
              onClick={() => setCurrentPage(currentPage - 1)}
              disabled={currentPage === 1}
              className="px-4 py-2 border border-gray-300 text-gray-700 rounded-lg hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Previous
            </button>
            <button
              onClick={() => setCurrentPage(currentPage + 1)}
              disabled={currentPage === totalPages}
              className="px-4 py-2 border border-gray-300 text-gray-700 rounded-lg hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Next
            </button>
          </div>
        </div>
      </div>

      {/* Edit User Modal */}
      {showEditModal && editingUser && (
        <div className="fixed inset-0 bg-gray-600 bg-opacity-50 overflow-y-auto h-full w-full z-50">
          <div className="relative top-20 mx-auto p-5 border w-96 shadow-lg rounded-md bg-white">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-medium text-gray-900">Edit User</h3>
              <button
                onClick={() => setShowEditModal(false)}
                className="text-gray-400 hover:text-gray-500"
              >
                <FaTimes />
              </button>
            </div>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-gray-700">Email</label>
                <input
                  type="email"
                  value={editingUser.email || ''}
                  onChange={e => setEditingUser({ ...editingUser, email: e.target.value })}
                  className="mt-1 w-full px-3 py-2 border border-gray-300 rounded-md focus:ring-green-500 focus:border-green-500"
                />
              </div>
              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-medium text-gray-700">First Name</label>
                  <input
                    type="text"
                    value={editingUser.firstName || ''}
                    onChange={e => setEditingUser({ ...editingUser, firstName: e.target.value })}
                    className="mt-1 w-full px-3 py-2 border border-gray-300 rounded-md focus:ring-green-500 focus:border-green-500"
                  />
                </div>
                <div>
                  <label className="block text-sm font-medium text-gray-700">Last Name</label>
                  <input
                    type="text"
                    value={editingUser.lastName || ''}
                    onChange={e => setEditingUser({ ...editingUser, lastName: e.target.value })}
                    className="mt-1 w-full px-3 py-2 border border-gray-300 rounded-md focus:ring-green-500 focus:border-green-500"
                  />
                </div>
              </div>
              <div>
                <label className="block text-sm font-medium text-gray-700">Phone</label>
                <input
                  type="text"
                  value={editingUser.phone || ''}
                  onChange={e => setEditingUser({ ...editingUser, phone: e.target.value })}
                  className="mt-1 w-full px-3 py-2 border border-gray-300 rounded-md focus:ring-green-500 focus:border-green-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-gray-700">Occupation</label>
                <input
                  type="text"
                  value={editingUser.occupation || ''}
                  onChange={e => setEditingUser({ ...editingUser, occupation: e.target.value })}
                  className="mt-1 w-full px-3 py-2 border border-gray-300 rounded-md focus:ring-green-500 focus:border-green-500"
                />
              </div>
            </div>
            <div className="mt-6 flex justify-end space-x-3">
              <button
                onClick={() => setShowEditModal(false)}
                className="px-4 py-2 border border-gray-300 text-gray-700 rounded-md hover:bg-gray-50"
              >
                Cancel
              </button>
              <button
                onClick={handleSaveEdit}
                className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 flex items-center"
              >
                <FaSave className="mr-2" />
                Save Changes
              </button>
            </div>
          </div>
        </div>
      )}

      {/* User Detail Modal */}
      {showUserDetail && selectedUser && (
        <div className="fixed inset-0 bg-gray-600 bg-opacity-50 overflow-y-auto h-full w-full z-50">
          <div className="relative top-20 mx-auto p-5 border w-[600px] shadow-lg rounded-md bg-white">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-medium text-gray-900">User Details</h3>
              <button
                onClick={() => setShowUserDetail(false)}
                className="text-gray-400 hover:text-gray-500"
              >
                <FaTimes />
              </button>
            </div>
            <div className="space-y-4">
              <div className="grid grid-cols-2 gap-4">
                <div>
                  <p className="text-sm font-medium text-gray-500">Username</p>
                  <p className="mt-1 text-sm text-gray-900">{selectedUser.username}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Email</p>
                  <p className="mt-1 text-sm text-gray-900">{selectedUser.email}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Full Name</p>
                  <p className="mt-1 text-sm text-gray-900">
                    {selectedUser.firstName} {selectedUser.lastName}
                  </p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Phone</p>
                  <p className="mt-1 text-sm text-gray-900">{selectedUser.phone || '-'}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Occupation</p>
                  <p className="mt-1 text-sm text-gray-900">{selectedUser.occupation || '-'}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Location</p>
                  <p className="mt-1 text-sm text-gray-900">
                    {selectedUser.city || '-'}, {selectedUser.state || '-'},{' '}
                    {selectedUser.country || '-'}
                  </p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Date Joined</p>
                  <p className="mt-1 text-sm text-gray-900">
                    {new Date(selectedUser.dateJoined).toLocaleDateString()}
                  </p>
                </div>
                <div>
                  <p className="text-sm font-medium text-gray-500">Total Conversations</p>
                  <p className="mt-1 text-sm text-gray-900">{selectedUser.totalConversations}</p>
                </div>
              </div>

              {selectedUser.recentConversations.length > 0 && (
                <div>
                  <h4 className="text-sm font-medium text-gray-700 mb-2">Recent Conversations</h4>
                  <div className="border rounded-lg divide-y">
                    {selectedUser.recentConversations.map(conv => (
                      <div key={conv.id} className="p-3">
                        <div className="flex justify-between items-start">
                          <div>
                            <p className="text-sm font-medium text-gray-900">
                              {conv.title || 'Untitled'}
                            </p>
                            <p className="text-xs text-gray-500">{conv.messageCount} messages</p>
                          </div>
                          <p className="text-xs text-gray-500">
                            {new Date(conv.updatedAt).toLocaleDateString()}
                          </p>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
