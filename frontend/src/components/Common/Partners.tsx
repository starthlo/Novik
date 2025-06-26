import { useState } from 'react';
import osstem_logo from '../../assets/ossem_logo.png';

const Partners = () => {
  const [clicks, setClicks] = useState(0);

  const handleClick = () => {
    setClicks(prev => prev + 1);
    fetch('/api/banner-click', { method: 'POST' });
  };

  return (
    <div className="my-8 text-center">
      <p className="text-black text-xl text-gray-500 font-medium mb-8">Our Partners</p>
      <a
        href="https://marketing.osstem.es/"
        target="_blank"
        rel="noopener noreferrer"
        onClick={handleClick}
      >
        <img
          src={osstem_logo}
          alt="Osstem"
          className="h-16 mx-auto cursor-pointer hover:scale-105 transition-transform"
        />
      </a>
      <p className="text-xs text-gray-200 mt-2">clicks: {clicks}</p>
    </div>
  );
};

export default Partners;
