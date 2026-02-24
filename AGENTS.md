# Global Agent Instructions

## Response Style

- Style: telegraph; drop filler/grammar; min tokens in all replies.
- Command sequencing: after any build command, serialize follow-up commands in separate calls; never chain build+run.
- Keep files <~500 LOC; split/refactor as needed.
- Prefer end-to-end verify; if blocked, say what’s missing.

## Git
- Destructive ops forbidden unless explicit (reset --hard, clean, restore, rm, …).
- Commits: Conventional Commits (feat|fix|refactor|build|ci|chore|docs|style|perf|test).
- Commit format: `type(scope): subject` (scope required when applicable).

## Atomic Commits
- One logical change per commit; if message needs "and", split commits.
- No mixed concerns: behavior change, refactor, formatting, and renames in separate commits.
- Stage surgically: use `git add -p`; review staged diff before commit.
- Each commit must pass relevant checks for touched code (tests/lint/build subset).
- Keep commit size reviewable; large/mechanical changes isolated in dedicated commit.
- After any successful code change + verify, always make atomic commit(s) immediately.
- If worktree has unrelated changes, stage only touched hunks/files.
- If verify is blocked/failing, do not commit; report blocker.
