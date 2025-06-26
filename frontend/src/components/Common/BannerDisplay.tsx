import { useState, useEffect, useRef } from 'react';

type Banner = {
  id: number;
  title: string;
  image_url: string | null;
  link: string | null;
  code: string | null;
  is_active: boolean;
};

export default function BannerDisplay() {
  const [banners, setBanners] = useState<Banner[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // Carousel state
  const [current, setCurrent] = useState(0);
  const intervalRef = useRef<NodeJS.Timer | undefined>(undefined);

  // fetch banners
  useEffect(() => {
    async function fetchBanners() {
      try {
        const res = await fetch('/api/banners/', { credentials: 'include' });
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const data: Banner[] = await res.json();
        setBanners(data.filter(b => b.is_active));
      } catch (e: any) {
        setError(e.message);
      } finally {
        setLoading(false);
      }
    }
    fetchBanners();
  }, []);

  // autoplay only if >1 banner
  useEffect(() => {
    if (banners.length <= 1) return;
    intervalRef.current = setInterval(() => {
      setCurrent(c => (c + 1) % banners.length);
    }, 5000);
    return () => {
      if (intervalRef.current) {
        clearInterval(intervalRef.current as unknown as number);
      }
    };
  }, [banners]);

  const prev = () => {
    if (intervalRef.current) {
      clearInterval(intervalRef.current as unknown as number);
    }
    setCurrent(c => (c - 1 + banners.length) % banners.length);
  };
  const next = () => {
    clearInterval(intervalRef.current as unknown as number);
    setCurrent(c => (c + 1) % banners.length);
  };

  const handleClick = async (b: Banner) => {
    if (!b.link) return;
    try {
      await fetch(`/api/banners/${b.id}/click/`, {
        method: 'POST',
        credentials: 'include',
      });
      window.open(b.link, '_blank');
    } catch {}
  };

  // record views once on load
  useEffect(() => {
    banners.forEach(b => {
      fetch(`/api/banners/${b.id}/view/`, {
        method: 'POST',
        credentials: 'include',
      }).catch(() => {});
    });
  }, [banners]);

  if (loading || error || banners.length === 0) return null;

  // Single banner, no slider
  if (banners.length === 1) {
    const b = banners[0];
    return (
      <div className="max-w-7xl mx-auto px-4 py-6">
        <div
          className="cursor-pointer rounded-xl overflow-hidden shadow-lg"
          onClick={() => handleClick(b)}
        >
          {b.image_url ? (
            <img src={b.image_url} alt={b.title} className="w-full h-auto object-cover" />
          ) : (
            <div className="w-full" dangerouslySetInnerHTML={{ __html: b.code || '' }} />
          )}
        </div>
      </div>
    );
  }

  // Carousel for multiple banners
  return (
    <div className="relative max-w-7xl mx-auto px-4 py-6 overflow-hidden">
      {/* Slider Track */}
      <div
        className="flex transition-transform duration-500 ease-out"
        style={{ transform: `translateX(-${current * 100}%)` }}
      >
        {banners.map(b => (
          <div
            key={b.id}
            className="flex-shrink-0 w-full px-10 cursor-pointer rounded-xl overflow-hidden"
            onClick={() => handleClick(b)}
          >
            {b.image_url ? (
              <img src={b.image_url} alt={b.title} className="w-full h-auto object-cover" />
            ) : (
              <div className="w-full" dangerouslySetInnerHTML={{ __html: b.code || '' }} />
            )}
          </div>
        ))}
      </div>

      {/* Prev/Next Buttons */}
      <button
        onClick={prev}
        className="absolute left-4 top-1/2 -translate-y-1/2 bg-black/30 hover:bg-black/50 text-white p-2 rounded-full"
        aria-label="Previous"
      >
        ‹
      </button>
      <button
        onClick={next}
        className="absolute right-4 top-1/2 -translate-y-1/2 bg-black/30 hover:bg-black/50 text-white p-2 rounded-full"
        aria-label="Next"
      >
        ›
      </button>

      {/* Dots */}
      <div className="absolute bottom-4 left-1/2 -translate-x-1/2 flex space-x-2">
        {banners.map((_, idx) => (
          <button
            key={idx}
            onClick={() => {
              clearInterval(intervalRef.current as unknown as number);
              setCurrent(idx);
            }}
            className={`w-2 h-2 rounded-full transition-all ${
              idx === current ? 'bg-white w-4' : 'bg-white/50 hover:bg-white'
            }`}
            aria-label={`Go to slide ${idx + 1}`}
          />
        ))}
      </div>
    </div>
  );
}
