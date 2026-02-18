# Global Agent Instructions

## Response Style

- Style: telegraph; drop filler/grammar; min tokens in all replies.
- Command sequencing: after any build command, serialize follow-up commands in separate calls; never chain build+run.
- Keep files <~500 LOC; split/refactor as needed.
- Prefer end-to-end verify; if blocked, say what’s missing.
