import { useState, useEffect, useRef, ChangeEvent } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkHeaderId from 'remark-heading-id';
import { v4 as uuidv4 } from 'uuid';
import { FaPaperclip } from 'react-icons/fa';

import Header from './Common/Header';
import PoweredBy from './Common/PoweredBy';
import Partners from './Common/Partners';
import Footer from './Common/Footer';

type ChatMessage = {
  question: string;
  answer: string;
  timestamp: Date;
};

function Dashboard() {
  const [input, setInput] = useState('');
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [chatHistory, setChatHistory] = useState<ChatMessage[]>([]);
  const [loading, setLoading] = useState(false);
  const [pendingQuestion, setPendingQuestion] = useState<string | null>(null);
  const bottomRef = useRef<HTMLDivElement | null>(null);
  const [isFooterVisible, setIsFooterVisible] = useState(false);
  const footerObserverRef = useRef<HTMLDivElement>(null);
  const inputWrapperRef = useRef<HTMLDivElement>(null);
  const [sessionId, setSessionId] = useState('');
  const textAreaRef = useRef<HTMLTextAreaElement | null>(null);

  useEffect(() => {
    resizeTextArea();
  }, [input]);

  useEffect(() => {
    const uuid = uuidv4();
    setSessionId(uuid);
    resizeTextArea();
    window.addEventListener('resize', resizeTextArea);

    const observer = new IntersectionObserver(
      ([entry]) => {
        setIsFooterVisible(entry.intersectionRatio > 0.1);
      },
      {
        threshold: 0.1,
        rootMargin: '50px',
      }
    );

    if (footerObserverRef.current) {
      observer.observe(footerObserverRef.current);
    }

    return () => {
      window.removeEventListener('resize', resizeTextArea);
      observer.disconnect();
    };
  }, []);

  const resizeTextArea = () => {
    if (!textAreaRef.current) {
      return;
    }

    // Reset height first to properly calculate scrollHeight
    textAreaRef.current.style.height = '48px';

    // Set to auto height based on content (if any)
    if (input !== '') {
      textAreaRef.current.style.height = `${textAreaRef.current.scrollHeight}px`;
    }
  };

  const scrollToBottom = (behavior: ScrollBehavior = 'smooth') => {
    bottomRef.current?.scrollIntoView({
      block: 'center',
      behavior,
    });
  };

  const handleFileSelect = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file || file.type !== 'application/pdf') return;
    setSelectedFile(file);
  };

  const handleSubmit = async () => {
    if (!input.trim() && !selectedFile) return;

    const currentInput = input;
    const currentTimestamp = new Date();
    const questionText = selectedFile ? `📄 ${selectedFile.name}: ${currentInput}` : currentInput;

    setInput('');
    setPendingQuestion(questionText);
    setLoading(true);

    // Scroll to show the pending question
    scrollToBottom();

    try {
      let endpoint = '/api/dashboard/';
      let body: FormData | string;

      if (selectedFile) {
        const formData = new FormData();
        formData.append('sessionId', sessionId);
        formData.append('pdf', selectedFile);
        formData.append('message', currentInput);
        body = formData;
        endpoint = '/api/dashboard/pdf/';
      } else {
        body = JSON.stringify({
          sessionId: sessionId,
          message: currentInput,
        });
      }

      const res = await fetch(endpoint, {
        method: 'POST',
        headers: selectedFile ? undefined : { 'Content-Type': 'application/json' },
        body: body,
      });

      let data;
      if (res.ok) {
        data = await res.json();
      } else {
        const errorMessage =
          res.status === 429
            ? 'Too many requests. Please wait a moment and try again.'
            : res.status >= 500
              ? 'Server error. Please try again later.'
              : 'Request failed. Please try again.';
        data = { message: errorMessage };
      }

      setChatHistory(prev => [
        ...prev,
        {
          question: questionText,
          answer: data.message,
          timestamp: currentTimestamp,
        },
      ]);

      // Clear the selected file after submission
      setSelectedFile(null);

      // Scroll to show the new response
      scrollToBottom();
    } catch (error) {
      console.error('API request failed:', error);
      setChatHistory(prev => [
        ...prev,
        {
          question: questionText,
          answer:
            'Error connecting to server. Please check your internet connection and try again.',
          timestamp: currentTimestamp,
        },
      ]);
    } finally {
      setPendingQuestion(null);
      setLoading(false);
    }
  };

  const handleKeyPress = (e: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit();
    }
  };

  return (
    <div className="bg-dental w-screen h-screen flex flex-col">
      <Header />
      <div className="flex flex-col flex-grow overflow-y-auto relative">
        <div className="flex-grow pb-32">
          <h1 className="text-2xl font-bold text-gray-500 mt-20 text-center">
            Your smart AI Dental assistant for safe clinical decisions
          </h1>

          {chatHistory.map((chat, index) => (
            <div key={index} className="w-full max-w-6xl mx-auto">
              <div className="flex justify-start mt-4">
                <div className="bg-white rounded-lg shadow-md p-4">
                  <pre className="text-gray-800 mb-2">{chat.question}</pre>
                  <span className="text-gray-500 text-sm">
                    {chat.timestamp.toLocaleDateString()} {chat.timestamp.toLocaleTimeString()}
                  </span>
                </div>
              </div>

              <div className="flex justify-end mt-2">
                <div className="bg-orange-400 rounded-lg shadow-md p-4 max-w-[60%] space-y-4">
                  {chat.answer.split('@@@@@').map((section, i) => {
                    const trimmed = section.trim();
                    if (!trimmed) return null;

                    // Track references section state
                    let isReferences: boolean = false;

                    return (
                      <div key={i} className="text-grey-900 p-2 rounded">
                        <div className="prose prose-sm">
                          <ReactMarkdown
                            remarkPlugins={[remarkHeaderId]}
                            components={{
                              h3: ({ node, ...props }) => {
                                // Check for references header in a more robust way
                                const content = props.children?.toString?.() || '';
                                isReferences = /^references$/i.test(content.trim());

                                return <h3 className="font-bold mt-4" {...props} />;
                              },
                              p: ({ node, ...props }) => {
                                let cs = 'mb-1';

                                if (isReferences) {
                                  cs += ' text-xs';
                                }
                                return <p className={cs} {...props} />;
                              },
                              a: ({ node, ...props }) => {
                                let pc: any = props.children;
                                if (pc) {
                                  // link as superscript
                                  if (pc[0] && typeof pc[0] === 'string' && pc[0].startsWith('^')) {
                                    return (
                                      <a
                                        className="text-blue-800 underline hover:text-blue-800"
                                        {...props}
                                      >
                                        <sup>{pc.slice(1, -1)}</sup>
                                      </a>
                                    );
                                  }
                                }
                                // link as usual
                                return (
                                  <a
                                    target="_blank"
                                    rel="noopener noreferrer"
                                    className="text-blue-800 underline hover:text-blue-800"
                                    {...props}
                                  />
                                );
                              },
                              ul: ({ node, ...props }) => (
                                <ul className="list-disc list-inside my-2" {...props} />
                              ),
                              li: ({ node, ...props }) => <li className="mb-1" {...props} />,
                              ol: ({ node, ...props }) => {
                                if (isReferences) {
                                  return (
                                    <ol
                                      className="list-decimal list-outside pl-5 space-y-2 text-sm text-gray-800"
                                      {...props}
                                    />
                                  );
                                }
                                return <ol {...props} />;
                              },
                            }}
                          >
                            {trimmed}
                          </ReactMarkdown>
                        </div>
                      </div>
                    );
                  })}

                  <span className="text-white text-sm block text-right">
                    {chat.timestamp.toLocaleDateString()} {chat.timestamp.toLocaleTimeString()}
                  </span>
                </div>
              </div>
            </div>
          ))}
          {pendingQuestion && (
            <div className="w-full max-w-6xl mx-auto">
              <div className="flex justify-start mt-4">
                <div className="bg-white rounded-lg shadow-md p-4">
                  <pre className="text-gray-800 mb-2">{pendingQuestion}</pre>
                  <span className="text-gray-500 text-sm">
                    {new Date().toLocaleDateString()} {new Date().toLocaleTimeString()}
                  </span>
                </div>
              </div>
            </div>
          )}

          {loading && (
            <div className="mt-4 flex flex-col items-center">
              <div className="w-8 h-8 border-4 border-gray-300 border-t-orange-500 rounded-full animate-spin"></div>
              <p className="text-gray-500 mt-2">Preparing response, please be patient...</p>
            </div>
          )}

          <div ref={bottomRef} />

          <PoweredBy />
          <Partners />
        </div>

        <div
          ref={inputWrapperRef}
          className={isFooterVisible ? 'sticky bottom-0 z-10 pb-20' : 'sticky bottom-0 z-10'}
        >
          <div className="w-full p-4">
            {/* plus button for file uploads */}
            <div className="flex items-end bg-white rounded-lg shadow-md overflow-hidden w-full max-w-6xl mx-auto">
              <label
                className={`cursor-pointer hover:text-orange-500 transition mb-3 px-3 flex items-center space-x-1 ${
                  selectedFile ? 'text-orange-500' : ''
                }`}
              >
                <FaPaperclip
                  className="h-5 w-5 text-xl text-black-500 hover:scale-120 transition cursor-pointer"
                  title="Upload file"
                />

                <input
                  type="file"
                  accept="application/pdf"
                  className="hidden"
                  onChange={handleFileSelect}
                />
              </label>

              <div className="flex-1 flex items-center">
                {selectedFile && (
                  <div className="flex items-center px-2 py-1 bg-orange-100 rounded-full mr-2 text-sm">
                    <span className="truncate max-w-xs">{selectedFile.name}</span>
                    <button
                      onClick={() => setSelectedFile(null)}
                      className="ml-1 text-gray-500 hover:text-gray-700"
                    >
                      ×
                    </button>
                  </div>
                )}
                <textarea
                  className="resize-none overflow-hidden w-full flex-1 p-3 outline-none text-gray-700"
                  placeholder={
                    selectedFile
                      ? 'Enter your question about the PDF...'
                      : "Please provide the patient's age, weight, medications, and treatment to be performed"
                  }
                  value={input}
                  ref={textAreaRef}
                  onChange={e => {
                    setInput(e.target.value);
                  }}
                  onKeyDown={e => handleKeyPress(e)}
                ></textarea>
              </div>

              <button
                onClick={handleSubmit}
                className="bg-orange-500 text-white px-4 py-2 mb-1 mx-1 rounded"
              >
                ↑
              </button>
            </div>
          </div>
        </div>

        <div ref={footerObserverRef}>
          <Footer />
        </div>
      </div>
    </div>
  );
}

export default Dashboard;
