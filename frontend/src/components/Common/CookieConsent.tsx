import { useEffect, useState } from 'react';

const CookieConsent = () => {
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const consent = localStorage.getItem('cookieConsent');
    if (!consent) setVisible(true);
  }, []);

  const handleConsent = (value: string) => {
    localStorage.setItem('cookieConsent', value);
    setVisible(false);
  };

  if (!visible) return null;

  return (
    <div className="fixed bottom-6 right-6 z-50">
      <div className="bg-white rounded-lg shadow-xl max-w-2xl w-full p-6 relative border border-gray-200">
        <button
          onClick={() => setVisible(false)}
          className="absolute top-2 right-2 text-gray-500 hover:text-black text-lg"
        >
          ×
        </button>
        <h2 className="text-lg font-semibold mb-2">Manage cookie consent</h2>
        <p className="text-gray-700 text-sm mb-4">
          We use cookies to enhance your experience. You can accept all, only essential cookies, or
          decline.
        </p>
        <div className="flex flex-wrap gap-2 justify-between mb-4">
          <button
            className="bg-orange-500 hover:bg-orange-600 text-white px-4 py-2 rounded font-medium"
            onClick={() => handleConsent('all')}
          >
            Accept
          </button>
          <button
            className="bg-gray-200 hover:bg-orange-600 px-4 py-2 rounded font-medium text-gray-800"
            onClick={() => handleConsent('functional')}
          >
            Accept only functional cookies
          </button>
          <button
            className="bg-gray-200 hover:bg-orange-60 px-4 py-2 rounded font-medium text-gray-800"
            onClick={() => handleConsent('decline')}
          >
            Decline
          </button>
        </div>
        <div className="flex justify-center gap-4 text-sm text-orange-500 underline">
          <a href="/legal#terms" target="_blank" rel="noopener noreferrer">
            Terms of service
          </a>
          <a href="/legal#privacy" target="_blank" rel="noopener noreferrer">
            Privacy Policy
          </a>
          <a href="/legal#cookies" target="_blank" rel="noopener noreferrer">
            Cookie Policy
          </a>
        </div>
      </div>
    </div>
  );
};

export default CookieConsent;
