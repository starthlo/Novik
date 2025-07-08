export type Message = {
  id: string;
  role: 'user' | 'assistant';
  content: string;
  timestamp?: Date;
  hasPdf?: boolean;
  pdfName?: string;
};

export type Conversation = {
  id: string;
  title: string;
  messages: Array<Message>;
  createdAt: string;
  updatedAt: string;
};
