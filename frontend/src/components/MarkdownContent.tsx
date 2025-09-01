import React from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkHeaderId from 'remark-heading-id';
import { Box } from '@mui/material';
import { novikTheme } from '../styles/theme';

interface MarkdownContentProps {
  content: string;
}

const MarkdownContent: React.FC<MarkdownContentProps> = ({ content }) => {
  return (
    <Box
      sx={{
        '& p': { mb: 1 },
        '& h1, & h2, & h3, & h4, & h5, & h6': { mt: 2, mb: 1 },
        '& ul, & ol': { pl: 2 },
        '& code': {
          backgroundColor: 'rgba(255, 255, 255, 0.1)',
          p: 0.5,
          borderRadius: 1,
          fontFamily: 'monospace',
        },
        '& sup a': {
          textDecoration: 'none',
          padding: '0 2px',
          borderRadius: '3px',
          fontWeight: 'bold',
          marginLeft: '2px',
          cursor: 'pointer',
        },
        '& [id^="eWzv"]': {
          display: 'block',
          margin: '10px 0',
          padding: '5px 10px',
          backgroundColor: 'rgba(255, 255, 255, 0.08)',
          borderLeft: '3px solid rgba(255, 255, 255, 0.2)',
          borderRadius: '3px',
          fontSize: '0.9em',
        },
      }}
    >
      <ReactMarkdown
        remarkPlugins={[remarkHeaderId, remarkGfm]}
        components={{
          // eslint-disable-next-line
          a: ({ node, ...props }) => {
            // Special handling for footnote links
            if (props.href && props.href.startsWith('#') && props.children) {
              // Safely handle various types of children
              const childText = Array.isArray(props.children)
                ? String(props.children[0] || '')
                : String(props.children || '');
              if (childText.startsWith('^') && childText.endsWith('^')) {
                // Extract the number between the carets
                const footnoteNumber = childText.substring(1, childText.length - 1);
                return (
                  <sup>
                    <a
                      {...props}
                      style={{
                        textDecoration: 'underline',
                        color: novikTheme.colors.primary,
                        cursor: 'pointer',
                        fontSize: '1em',
                      }}
                    >
                      {footnoteNumber}
                    </a>
                  </sup>
                );
              }
            }
            return (
              <a
                {...props}
                target="_blank"
                rel="noopener noreferrer"
                style={{
                  color: novikTheme.colors.primary,
                  textDecoration: 'underline',
                  fontWeight: 500,
                }}
              />
            );
          },
          // eslint-disable-next-line
          h6: ({ node, ...props }) => {
            // Special handling for footnote section headings
            if (props.id) {
              return (
                <h6
                  {...props}
                  style={{
                    color: novikTheme.colors.primary,
                  }}
                />
              );
            }
            return <h6 {...props} />;
          },
        }}
      >
        {content}
      </ReactMarkdown>
    </Box>
  );
};

export default MarkdownContent;
