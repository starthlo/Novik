import React from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkBreaks from 'remark-breaks';
import { styled } from '@mui/material/styles';
import {
  Box,
  Link,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Divider,
} from '@mui/material';
import { novikTheme } from '../styles/theme';

interface MarkdownContentProps {
  content: string;
  isUser?: boolean;
}

const MarkdownWrapper = styled(Box)<{ isUser?: boolean }>(({ isUser }) => ({
  '& p': {
    margin: '0.5em 0',
    lineHeight: 1.6,
    color: isUser ? '#ffffff' : novikTheme.colors.text,
  },
  '& h1, & h2, & h3, & h4, & h5, & h6': {
    marginTop: '1em',
    marginBottom: '0.5em',
    fontWeight: 600,
    lineHeight: 1.3,
    color: isUser ? '#ffffff' : novikTheme.colors.text,
  },
  '& h1': { fontSize: '1.5em' },
  '& h2': { fontSize: '1.3em' },
  '& h3': { fontSize: '1.15em' },
  '& h4': { fontSize: '1.1em' },
  '& h5': { fontSize: '1.05em' },
  '& h6': { fontSize: '1em' },
  '& ul, & ol': {
    marginTop: '0.5em',
    marginBottom: '0.5em',
    paddingLeft: '1.5em',
  },
  '& li': {
    marginTop: '0.25em',
    marginBottom: '0.25em',
    lineHeight: 1.6,
    color: isUser ? '#ffffff' : novikTheme.colors.text,
  },
  '& ul li': {
    listStyleType: 'disc',
  },
  '& ol li': {
    listStyleType: 'decimal',
  },
  '& blockquote': {
    borderLeft: `3px solid ${isUser ? 'rgba(255,255,255,0.5)' : novikTheme.colors.primary}`,
    paddingLeft: '1em',
    marginLeft: 0,
    marginRight: 0,
    marginTop: '0.5em',
    marginBottom: '0.5em',
    fontStyle: 'italic',
    color: isUser ? 'rgba(255,255,255,0.9)' : novikTheme.colors.textMuted,
  },
  '& code': {
    backgroundColor: isUser ? 'rgba(0,0,0,0.2)' : 'rgba(0,0,0,0.05)',
    padding: '0.2em 0.4em',
    borderRadius: '3px',
    fontSize: '0.9em',
    fontFamily: '"Courier New", Courier, monospace',
    color: isUser ? '#ffffff' : novikTheme.colors.text,
  },
  '& pre': {
    backgroundColor: isUser ? 'rgba(0,0,0,0.3)' : '#f5f5f5',
    padding: '1em',
    borderRadius: '6px',
    overflowX: 'auto',
    marginTop: '0.5em',
    marginBottom: '0.5em',
  },
  '& pre code': {
    backgroundColor: 'transparent',
    padding: 0,
    fontSize: '0.85em',
    lineHeight: 1.5,
  },
  '& table': {
    width: '100%',
    borderCollapse: 'collapse',
    marginTop: '0.5em',
    marginBottom: '0.5em',
  },
  '& th, & td': {
    padding: '0.5em',
    textAlign: 'left',
    borderBottom: `1px solid ${isUser ? 'rgba(255,255,255,0.2)' : novikTheme.colors.border}`,
  },
  '& th': {
    fontWeight: 600,
    backgroundColor: isUser ? 'rgba(0,0,0,0.2)' : '#f5f5f5',
  },
  '& a': {
    color: isUser ? '#ffffff' : novikTheme.colors.primary,
    textDecoration: 'underline',
    '&:hover': {
      textDecoration: 'none',
      opacity: 0.8,
    },
  },
  '& hr': {
    border: 'none',
    borderTop: `1px solid ${isUser ? 'rgba(255,255,255,0.2)' : novikTheme.colors.border}`,
    marginTop: '1em',
    marginBottom: '1em',
  },
  '& strong': {
    fontWeight: 600,
  },
  '& em': {
    fontStyle: 'italic',
  },
  '& img': {
    maxWidth: '100%',
    height: 'auto',
    borderRadius: '6px',
    marginTop: '0.5em',
    marginBottom: '0.5em',
  },
  // Specific styles for nested lists
  '& ul ul, & ol ol, & ul ol, & ol ul': {
    marginTop: '0.25em',
    marginBottom: '0.25em',
  },
  // Style for task lists (GitHub-style checkboxes)
  '& input[type="checkbox"]': {
    marginRight: '0.5em',
  },
  // Style for definition lists
  '& dl': {
    marginTop: '0.5em',
    marginBottom: '0.5em',
  },
  '& dt': {
    fontWeight: 600,
    marginTop: '0.5em',
  },
  '& dd': {
    marginLeft: '1.5em',
    marginBottom: '0.5em',
  },
  // Style for footnotes
  '& sup': {
    fontSize: '0.75em',
    verticalAlign: 'super',
  },
  '& .footnotes': {
    marginTop: '2em',
    paddingTop: '1em',
    borderTop: `1px solid ${isUser ? 'rgba(255,255,255,0.2)' : novikTheme.colors.border}`,
    fontSize: '0.9em',
  },
  // Style for PubMed references
  '& h6[id*="pubmed"], & h6[id*="reference"]': {
    fontSize: '0.9em',
    marginTop: '0.75em',
    color: isUser ? 'rgba(255,255,255,0.9)' : novikTheme.colors.textMuted,
  },
}));

const MarkdownContent: React.FC<MarkdownContentProps> = ({ content, isUser = false }) => {
  return (
    <MarkdownWrapper isUser={isUser}>
      <ReactMarkdown
        remarkPlugins={[remarkGfm, remarkBreaks]}
        components={{
          // Custom link component
          a: ({ href, children }) => (
            <Link
              href={href}
              target="_blank"
              rel="noopener noreferrer"
              sx={{
                color: isUser ? '#ffffff !important' : `${novikTheme.colors.primary} !important`,
                textDecoration: 'underline',
                '&:hover': {
                  textDecoration: 'none',
                  opacity: 0.8,
                },
              }}
            >
              {children}
            </Link>
          ),
          // Custom paragraph component
          p: ({ children }) => (
            <Typography
              component="p"
              sx={{
                margin: '0.5em 0',
                lineHeight: 1.6,
                color: 'inherit',
              }}
            >
              {children}
            </Typography>
          ),
          // Custom table components
          table: ({ children }) => (
            <Table size="small" sx={{ marginY: 1 }}>
              {children}
            </Table>
          ),
          thead: ({ children }) => <TableHead>{children}</TableHead>,
          tbody: ({ children }) => <TableBody>{children}</TableBody>,
          tr: ({ children }) => <TableRow>{children}</TableRow>,
          th: ({ children }) => (
            <TableCell
              sx={{
                fontWeight: 600,
                backgroundColor: isUser ? 'rgba(0,0,0,0.2)' : '#f5f5f5',
                color: isUser ? '#ffffff' : novikTheme.colors.text,
              }}
            >
              {children}
            </TableCell>
          ),
          td: ({ children }) => (
            <TableCell
              sx={{
                color: isUser ? '#ffffff' : novikTheme.colors.text,
              }}
            >
              {children}
            </TableCell>
          ),
          // Custom horizontal rule
          hr: () => (
            <Divider
              sx={{
                marginY: 2,
                borderColor: isUser ? 'rgba(255,255,255,0.2)' : novikTheme.colors.border,
              }}
            />
          ),
          // Custom code components
          code: ({ children, ...props }) => {
            return (
              <pre
                style={{
                  backgroundColor: isUser ? 'rgba(0,0,0,0.3)' : '#f5f5f5',
                  padding: '1em',
                  borderRadius: '6px',
                  overflowX: 'auto',
                  marginTop: '0.5em',
                  marginBottom: '0.5em',
                }}
              >
                <code
                  style={{
                    fontSize: '0.85em',
                    lineHeight: 1.5,
                    fontFamily: '"Courier New", Courier, monospace',
                  }}
                  {...props}
                >
                  {children}
                </code>
              </pre>
            );
          },
        }}
      >
        {content}
      </ReactMarkdown>
    </MarkdownWrapper>
  );
};

export default MarkdownContent;
