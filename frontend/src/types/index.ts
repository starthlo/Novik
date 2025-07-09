export type Message = {
  id: string;
  role: 'user' | 'assistant';
  content: string;
  timestamp?: Date;
  file?: {
    fileName: string;
    text: string;
  };
};

export type Conversation = {
  id: string;
  title: string;
  messages: Array<Message>;
  createdAt: string;
  updatedAt: string;
};
