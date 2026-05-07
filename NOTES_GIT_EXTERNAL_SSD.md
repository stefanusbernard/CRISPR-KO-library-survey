# Git on External SSD — Important Notes

This repository lives on an external SSD and is used across multiple machines.
This has caused silent mass-deletion staging in the past. Read this before doing anything on a new machine.

## What went wrong before

All project files (scripts, figures, data) appeared **staged for deletion** in git, even though
the files were still on disk. This happened because macOS (or an IDE's git integration) silently
ran the equivalent of `git rm --cached` when the drive was mounted on a different machine.
No data was lost, but committing in that state would have wiped the entire git history.

## Rules when switching machines

1. **Run `git status` before opening any IDE or running any git command.**
   If you see staged deletions, immediately run:
   ```
   git restore --staged .
   ```

2. **Always eject the SSD properly** before unplugging. This ensures `.git/index` is fully
   flushed to disk and not left in a half-written state.

3. **Close the project in VS Code / RStudio before ejecting** on one machine.
   Open it only after the drive has fully mounted on the next machine.

## Why this happens

The SSD is formatted as **exFAT** (cross-platform format). When macOS mounts an exFAT drive,
it creates hidden `._*` resource fork files alongside your real files. IDE git panels can
misinterpret these or the index state and auto-stage deletions before you notice.

The `.gitignore` in this repo now includes `._*` and `.DS_Store` to suppress macOS metadata
files from ever being tracked.

## Quick reference

| Situation | Command |
|---|---|
| Check state before doing anything | `git status` |
| Undo staged deletions | `git restore --staged .` |
| Confirm repo matches remote | `git fetch origin && git log --oneline HEAD..origin/main` |
