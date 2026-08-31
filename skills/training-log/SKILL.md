---
name: training-log
description: Maintain this repository's daily XCPC problem archive, including accepted source code, reflective notes, daily summaries, and topic indexes. Use when recording, updating, or reorganizing competitive-programming training completed by the user.
---

# XCPC Training Log

Archive training under `training/YYYY/MM/YYYY-MM-DD/<platform><problem-id>-<short-name>/`.

For each problem:

- When the user provides an implementation, put its cleaned, independently compilable accepted version in `solution.cpp`.
- Do not author `solution.cpp` on the user's behalf unless the user explicitly asks for an implementation. If no user implementation is available, archive the notes without `solution.cpp`, mark the code as pending, and add it later when the user supplies it.
- Put the learning record in `notes.md`. Reconstruct the user's actual thought process from the conversation when available; do not replace it with a generic editorial.
- Add `attempt.cpp`, generators, or brute-force programs only when they have lasting diagnostic value.
- Preserve the user's implementation style where reasonable, while removing local-only includes, forced file redirection, debug output, and other submission hazards from `solution.cpp`.

Maintain `training/PENDING.md` as the repository-wide code backlog. Add a problem only when the user has already worked out a complete solution verbally but has not yet provided code; do not use it for unsolved or merely attempted problems. When archiving notes without `solution.cpp`, add or update its backlog entry. Remove the entry after the user's implementation has been archived and verified.

The notes should record the problem link and metadata when known, the user's initial approach, the point where they became stuck, the key observation, the final algorithm, correctness reasoning, complexity, mistakes or subtle implementation points, and dated revisit entries. Clearly distinguish what the user discovered independently from hints or editorial knowledge.

When an agent helped teach the problem, include a concise `引导过程` section. Preserve the progression rather than only the final explanation: record the user's current obstacle, the small question or hint the agent gave next, the insight the user then reached, and any remaining obstacle that led to the next hint. Condense repeated exchanges and paraphrase instead of transcribing the conversation. The section should make the gradual path visible without becoming a long chat log. Do not invent user discoveries or claim that the user derived a step independently when it was directly supplied.

After changing a problem archive:

1. Update that day's `README.md` with one concise table row.
2. Update `training/PENDING.md` when the problem's code-pending status changes.
3. Update each relevant `topics/*.md` index with a link and the reusable insight. Use multiple indexes instead of duplicating the problem directory when a problem has several tags.
4. Update the root `README.md` only when repository-wide navigation or conventions change.
5. Check relative links and compile `solution.cpp` when one was added or changed.

Use Chinese for learning notes and indexes unless the surrounding content establishes another language. Prefer exact dates and stable platform IDs in paths. Do not manufacture solve time, submission status, or independence claims that the conversation does not establish; mark unknown fields as `未记录`.
