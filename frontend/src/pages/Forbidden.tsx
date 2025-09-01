import { useNavigate } from 'react-router-dom';
import { FaExclamationTriangle, FaHome, FaArrowLeft } from 'react-icons/fa';

export default function Forbidden() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-gray-50 flex flex-col justify-center py-12 sm:px-6 lg:px-8">
      <div className="mt-8 sm:mx-auto sm:w-full sm:max-w-md">
        <div className="bg-white py-8 px-4 shadow-xl rounded-lg sm:px-10">
          <div className="text-center">
            <div className="mx-auto flex items-center justify-center h-16 w-16 rounded-full bg-red-100 mb-4">
              <FaExclamationTriangle className="h-8 w-8 text-red-600" />
            </div>

            <h2 className="text-3xl font-bold text-gray-900 mb-2">403 - Access Forbidden</h2>

            <p className="text-gray-600 mb-8">
              You don't have permission to access this page. This area is restricted to
              administrators only.
            </p>

            <div className="space-y-3">
              <button
                onClick={() => navigate(-1)}
                className="w-full flex items-center justify-center px-4 py-2 border border-gray-300 rounded-md shadow-sm text-sm font-medium text-gray-700 bg-white hover:bg-gray-50 transition duration-150"
              >
                <FaArrowLeft className="mr-2" />
                Go Back
              </button>

              <button
                onClick={() => navigate('/dashboard')}
                className="w-full flex items-center justify-center px-4 py-2 border border-transparent rounded-md shadow-sm text-sm font-medium text-white bg-green-600 hover:bg-green-700 transition duration-150"
              >
                <FaHome className="mr-2" />
                Go to Dashboard
              </button>
            </div>

            <div className="mt-6 text-sm text-gray-500">
              If you believe you should have access to this page, please contact your system
              administrator.
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
