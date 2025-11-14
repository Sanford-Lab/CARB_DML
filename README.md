# CARB

Working repo for the Sanford Lab CARB project. The repo is currently a clean slate with data stored outside of version control so analyses can evolve without exposing raw inputs.

## Repository Structure

- `code/` – core analysis scripts, functions, and any shared helpers (prefer tidyverse style).
- `data/input/` & `data/output/` – raw/confidential inputs and intermediate exports; the entire `data/` tree stays out of git.
- `writing/` – manuscripts, memos, slide decks, or other narrative artifacts tied to the project.
- `README.md` – project orientation, collaboration norms (this file).

Add additional subdirectories (e.g., `code/cleaning/`, `writing/figures/`) as the project grows; keep sensitive or bulky assets under `data/` unless they are safe to version.

## Git & Collaboration Workflow

- Clone: `git clone https://github.com/Sanford-Lab/CARB_DML.git`
- Keep `main` clean. Create feature branches (`git checkout -b feature-name`) for changes, then open pull requests.
- Before starting work, sync with `main`: `git pull origin main`.
- Stage and commit small logical chunks with clear messages.
- Push branches to origin and use GitHub pull requests to request reviews. Mention relevant collaborators.
- When reviewing, prefer tidyverse-style R code (pipe-friendly, snake_case objects, consistent indentation).
- If you need to share large data, place it in `data/` (ignored) and coordinate via shared drives instead of git.

## Getting Started

1. Ensure you have Git, R (with tidyverse), and the GitHub CLI installed.
2. Configure GitHub credentials (Keychain helper is already enabled on this machine).
3. Add any initial scripts or documentation, then `git add`, `git commit`, and `git push origin <branch>`.

Update this README as the repository structure and workflow evolve. yay.

## Cleaning Files 


## Analysis Files

## To Do List

- [ ] Create a to do list

