import { useState, useEffect } from 'react';
import { FaBan, FaCheckCircle, FaTrash } from 'react-icons/fa';
import Header from './Common/Header';

type User = {
  id: number;
  username: string;
  email: string;
  dob: string;
  phone: string;
  occupation: string;
  country: string;
  state: string;
  city: string;
  is_staff: boolean;
  isSuperuser: boolean;
  is_active: boolean;
  date_joined: string;
};

export default function UserManagement() {
  const [users, setUsers] = useState<User[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetch('/api/users/', { credentials: 'include' })
      .then(res => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
      })
      .then(data => {
        setUsers(data.users);
        setLoading(false);
      })
      .catch(err => {
        setError(err.message);
        setLoading(false);
      });
  }, []);

  if (loading) return <div>Loading users…</div>;
  if (error) return <div style={{ color: 'red' }}>Error: {error}</div>;

  const total = users.length;
  const byCountry = users.reduce(
    (acc, u) => {
      acc[u.country || 'Unknown'] = (acc[u.country || 'Unknown'] || 0) + 1;
      return acc;
    },
    {} as Record<string, number>
  );

  const toggleActiveStatus = async (id: number) => {
    try {
      const res = await fetch('/api/user/toggle-active-status/', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ id: id }),
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      if (res.ok) {
        const nextUsers = users.map(user => {
          if (user.id === id) {
            // Toggle is_active & return
            return {
              ...user,
              is_active: !user.is_active,
            };
          } else {
            // No change
            return user;
          }
        });
        setUsers(nextUsers);
      }
    } catch (err: any) {
      console.error(err);
    }
  };

  const trashUser = async (id: number) => {
    if (window.confirm('Are you sure you want to delete this user?')) {
      try {
        const res = await fetch('/api/user/trash/', {
          method: 'DELETE',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ id: id }),
        });

        if (!res.ok) throw new Error(`HTTP ${res.status}`);

        if (res.ok) {
          const nextUsers = users.filter(function (user) {
            return user.id !== id;
          });
          setUsers(nextUsers);
        }
      } catch (err: any) {
        console.error(err);
      }
    } else {
      // do nothing
    }
  };

  return (
    <>
      <Header />
      <div className="p-5 max-w-4xl mx-auto">
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-2xl font-semibold">User Management</h2>
          <button
            onClick={() => (window.location.href = '/api/users/export/')}
            className="bg-green-500 text-white px-4 py-2 rounded cursor-pointer"
          >
            Export CSV
          </button>
        </div>

        <div className="mb-6">
          <div>
            <strong>Total users:</strong> {total}
          </div>
          <div>
            <strong>Users by country:</strong>
          </div>
          <ul className="list-disc ml-6">
            {Object.entries(byCountry).map(([c, n]) => (
              <li key={c}>
                {c}: {n}
              </li>
            ))}
          </ul>
        </div>

        <div className="accordion-container pt-5">
          {users.map(u => (
            <details className="user-accordion" key={u.id}>
              <summary>
                <div className="w-full flex justify-between">
                  <div>{u.username}</div>
                  <div className="text-right text-sm mr-2">{u.email}</div>
                </div>
              </summary>
              <div className="user-details">
                <div className="detail-item">
                  <strong>ID</strong>
                  <span>{u.id}</span>
                </div>
                <div className="detail-item">
                  <strong>Joined</strong>
                  <span>{u.date_joined}</span>
                </div>
                <div className="detail-item">
                  <strong>DOB</strong>
                  <span>{u.dob}</span>
                </div>
                <div className="detail-item">
                  <strong>Phone</strong>
                  <span>{u.phone}</span>
                </div>
                <div className="detail-item">
                  <strong>Occupation</strong>
                  <span>{u.occupation}</span>
                </div>
                <div className="detail-item">
                  <strong>Country</strong>
                  <span>{u.country}</span>
                </div>
                <div className="detail-item">
                  <strong>State</strong>
                  <span>{u.state}</span>
                </div>
                <div className="detail-item">
                  <strong>City</strong>
                  <span>{u.city}</span>
                </div>
                <div className="detail-item">
                  <strong>Staff?</strong>
                  <span>{u.is_staff ? 'Yes' : 'No'}</span>
                </div>
                <div className="detail-item">
                  <strong>Superuser?</strong>
                  <span>{u.isSuperuser ? 'Yes' : 'No'}</span>
                </div>

                <div className="detail-item">&nbsp;</div>

                <div className="detail-item">
                  <div className="grid grid-flow-col justify-end gap-x-4">
                    {u.is_active ? (
                      <FaCheckCircle
                        className="text-xl text-green-500 hover:scale-120 transition cursor-pointer"
                        title="User is active. Click to ban"
                        onClick={() => toggleActiveStatus(u.id)}
                      />
                    ) : (
                      <FaBan
                        className="text-xl text-red-500 hover:scale-120 transition cursor-pointer"
                        title="User is banned. Click to unban"
                        onClick={() => toggleActiveStatus(u.id)}
                      />
                    )}
                    <FaTrash
                      className="text-xl text-red-500 hover:scale-110 transition cursor-pointer"
                      title="Delete"
                      onClick={() => trashUser(u.id)}
                    />
                  </div>
                </div>
              </div>
            </details>
          ))}
        </div>
      </div>
    </>
  );
}
