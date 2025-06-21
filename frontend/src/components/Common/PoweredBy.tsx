import chatgpt_logo from "../../assets/chatgpt.png";
import drugbank_logo from "../../assets/drugbank.png";
import pubmed_logo from "../../assets/pubmed.png";

const PoweredBy = () => (
  <div className="mt-24 text-center px-4">
    <h2 className="text-black text-xl sm:text-2xl text-gray-500 font-medium mb-4">Powered By</h2>

    <div className="flex flex-col sm:flex-row justify-center items-center gap-6 sm:gap-8 flex-wrap">
      <img
        src={chatgpt_logo}
        alt="ChatGPT"
        className="w-24 sm:w-32 md:w-48 object-contain"
      />
      <img
        src={drugbank_logo}
        alt="DrugBank"
        className="w-24 sm:w-32 md:w-48 object-contain"
      />
      <img
        src={pubmed_logo}
        alt="PubMed"
        className="w-24 sm:w-32 md:w-48 object-contain"
      />
    </div>
  </div>
);

export default PoweredBy;