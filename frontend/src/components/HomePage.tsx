import Header from "./Common/Header";
import PoweredBy from "./Common/PoweredBy";
import Partners from "./Common/Partners";
import Footer from "./Common/Footer";
import CookieConsent from "./Common/CookieConsent";
import NovikLogo from "../assets/Novik.png";
// import BannerDisplay from "./Common/BannerDisplay";

function HomePage() {
  return (
    <div className="bg-dental w-screen h-screen flex flex-col">
      <Header />
      <CookieConsent />
      <div className="flex flex-col flex-grow">
        {/* <BannerDisplay /> */}
        <div className="flex-grow overflow-y-auto flex flex-col items-center pt-24">
          <img src={NovikLogo} alt="Novik" className="h-20 mb-4" />
          <h1 className="text-xl text-gray-700 text-center">
            Your smart AI Dental assistant for safe clinical decisions
          </h1>
        </div>
        <PoweredBy />
        <Partners />
        <Footer />
      </div>
    </div>
  );
}

export default HomePage;
