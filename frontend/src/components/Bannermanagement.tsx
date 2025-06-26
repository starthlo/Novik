import { useState, useEffect, ChangeEvent, FormEvent } from 'react';
import Header from './Common/Header';

type Banner = {
  id: number;
  title: string;
  image: string | null;
  image_url: string | null;
  link: string | null;
  code: string | null;
  is_active: boolean;
};
type BannerStat = {
  id: number;
  banner: number;
  date: string;
  country: string;
  views: number;
  clicks: number;
};

export default function BannerManagement() {
  const [banners, setBanners] = useState<Banner[]>([]);
  const [stats, setStats] = useState<BannerStat[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [form, setForm] = useState<Partial<Banner>>({});
  const [file, setFile] = useState<File | null>(null);
  const [editingId, setEditingId] = useState<number | null>(null);

  // fetch list
  const fetchBanners = async () => {
    try {
      const res = await fetch('/api/banners/', { credentials: 'include' });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      setBanners(data);
    } catch (err: any) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  // fetch stats for one banner
  const fetchStats = async (bid: number) => {
    try {
      const res = await fetch(`/api/banner-stats/?banner=${bid}`, {
        credentials: 'include',
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      setStats(await res.json());
    } catch (err: any) {
      console.error(err);
    }
  };

  useEffect(() => {
    fetchBanners();
  }, []);

  // form handlers
  const handleInput = (e: ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value, type, checked } = e.target as HTMLInputElement;
    setForm(prev => ({
      ...prev,
      [name]: type === 'checkbox' ? checked : value,
    }));
  };
  const handleFile = (e: ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) setFile(e.target.files[0]);
  };

  // create/update
  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const fd = new FormData();

    // Add form fields to FormData
    if (form.title) fd.append('title', form.title);
    if (form.link) fd.append('link', form.link);
    if (form.code) fd.append('code', form.code);
    if (form.is_active !== undefined) fd.append('is_active', form.is_active.toString());

    // Add file if it exists
    if (file) {
      fd.append('image', file);
    }

    const url = editingId ? `/api/banners/${editingId}/` : '/api/banners/';
    const method = editingId ? 'PUT' : 'POST';

    try {
      const res = await fetch(url, {
        method,
        credentials: 'include',
        body: fd,
        // Don't set Content-Type header - browser will set it with boundary
      });

      if (!res.ok) {
        const errorData = await res.json();
        throw new Error(errorData.detail || `HTTP ${res.status}`);
      }

      // Reset form
      setForm({});
      setFile(null);
      setEditingId(null);
      fetchBanners();
    } catch (err: any) {
      setError(err.message);
    }
  };

  const handleEdit = (b: Banner) => {
    setEditingId(b.id);
    setForm({
      title: b.title,
      link: b.link || '',
      code: b.code || '',
      is_active: b.is_active,
    });
  };

  const handleDelete = async (id: number) => {
    if (!confirm('Delete this banner?')) return;
    try {
      const res = await fetch(`/api/banners/${id}/`, {
        method: 'DELETE',
        credentials: 'include',
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      fetchBanners();
    } catch (err: any) {
      setError(err.message);
    }
  };

  if (loading)
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-orange-500"></div>
      </div>
    );

  if (error)
    return (
      <div className="p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 mx-auto max-w-4xl mt-4">
        Error: {error}
      </div>
    );

  return (
    <>
      <Header />
      <div className="p-4 md:p-6 lg:p-8 max-w-6xl mx-auto">
        {/* Form Section */}
        <div className="bg-white rounded-xl shadow-sm border p-6 mb-8">
          <h2 className="text-2xl font-semibold mb-6 text-gray-800">
            {editingId ? 'Edit Banner' : 'Create New Banner'}
          </h2>

          <form onSubmit={handleSubmit} className="space-y-6">
            <div className="grid gap-6 md:grid-cols-2">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">Title *</label>
                <input
                  name="title"
                  value={form.title || ''}
                  onChange={handleInput}
                  className="w-full py-2 px-2 rounded-lg border-gray-300 shadow-sm focus:border-orange-500 focus:ring-orange-500"
                  required
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">Link URL</label>
                <input
                  name="link"
                  value={form.link || ''}
                  onChange={handleInput}
                  className="w-full py-2 px-1 rounded-lg border-gray-300 shadow-sm focus:border-orange-500 focus:ring-orange-500"
                  placeholder="https://"
                />
              </div>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">HTML Code</label>
              <textarea
                name="code"
                value={form.code || ''}
                onChange={handleInput}
                rows={4}
                className="w-full py-2 px-2 rounded-lg border-gray-300 shadow-sm focus:border-orange-500 focus:ring-orange-500"
                placeholder="<div>Your banner HTML here</div>"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">Banner Image</label>
              <div className="mt-1 flex items-center">
                <input
                  type="file"
                  onChange={handleFile}
                  accept="image/*"
                  className="block w-full text-sm text-gray-500
                    file:mr-4 file:py-2 file:px-4
                    file:rounded-full file:border-0
                    file:text-sm file:font-semibold
                    file:bg-orange-50 file:text-orange-700
                    hover:file:bg-orange-100"
                />
              </div>
            </div>

            <div className="flex items-center">
              <label className="inline-flex items-center">
                <input
                  type="checkbox"
                  name="is_active"
                  checked={form.is_active || false}
                  onChange={handleInput}
                  className="rounded border-gray-300 text-orange-600 focus:ring-orange-500"
                />
                <span className="ml-2 text-sm text-gray-700">Active</span>
              </label>
            </div>

            <div className="flex justify-end">
              <button
                type="submit"
                className="bg-orange-500 hover:bg-orange-600 text-white px-6 py-2 rounded-lg font-medium transition-colors"
              >
                {editingId ? 'Update Banner' : 'Create Banner'}
              </button>
            </div>
          </form>
        </div>

        {/* Banners List */}
        <div className="bg-white rounded-xl shadow-sm border overflow-hidden">
          <h2 className="text-2xl font-semibold p-6 border-b bg-gray-50">Manage Banners</h2>

          <div className="overflow-x-auto">
            <table className="w-full">
              <thead className="bg-gray-50">
                <tr>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                    Title
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                    Preview
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                    Status
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                    Actions
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                    Stats
                  </th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-200">
                {banners.map(b => (
                  <tr key={b.id} className="hover:bg-gray-50">
                    <td className="px-6 py-4 whitespace-nowrap">{b.title}</td>
                    <td className="px-6 py-4">
                      {b.image_url ? (
                        <img
                          src={b.image_url}
                          alt={b.title}
                          className="h-20 w-auto object-contain rounded"
                        />
                      ) : (
                        <div
                          className="max-w-xs overflow-hidden"
                          dangerouslySetInnerHTML={{ __html: b.code || '' }}
                        />
                      )}
                    </td>
                    <td className="px-6 py-4">
                      <span
                        className={`inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium
                        ${
                          b.is_active ? 'bg-green-100 text-green-800' : 'bg-gray-100 text-gray-800'
                        }`}
                      >
                        {b.is_active ? 'Active' : 'Inactive'}
                      </span>
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap">
                      <div className="flex space-x-3">
                        <button
                          onClick={() => handleEdit(b)}
                          className="text-blue-600 hover:text-blue-900 font-medium cursor-pointer"
                        >
                          Edit
                        </button>
                        <button
                          onClick={() => handleDelete(b.id)}
                          className="text-red-200 hover:text-red-900 font-medium cursor-pointer"
                        >
                          Delete
                        </button>
                      </div>
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap">
                      <button
                        onClick={() => fetchStats(b.id)}
                        className="text-blue-600 hover:text-blue-900 font-medium cursor-pointer"
                      >
                        View Stats
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {/* Stats Section */}
        {stats.length > 0 && (
          <div className="mt-8 bg-white rounded-xl shadow-sm border overflow-hidden">
            <h2 className="text-2xl font-semibold p-6 border-b bg-gray-50">Banner Statistics</h2>

            <div className="overflow-x-auto">
              <table className="w-full">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                      Date
                    </th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                      Country
                    </th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                      Views
                    </th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                      Clicks
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-200">
                  {stats.map(s => (
                    <tr key={s.id} className="hover:bg-gray-50">
                      <td className="px-6 py-4 whitespace-nowrap">{s.date}</td>
                      <td className="px-6 py-4 whitespace-nowrap">{s.country}</td>
                      <td className="px-6 py-4 whitespace-nowrap">{s.views}</td>
                      <td className="px-6 py-4 whitespace-nowrap">{s.clicks}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    </>
  );
}
