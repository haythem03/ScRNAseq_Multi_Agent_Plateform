---
description: Commit and push changes to GitHub
---

# Git Commit Workflow

Use this workflow when the user asks to commit changes to GitHub.

## Prerequisites
- Git must be installed (`git --version` to check)
- Repository must be initialized with `git init`
- Remote must be configured with `git remote add origin <url>`

## Steps

// turbo-all

1. Check git status to see what files have changed:
```bash
git status
```

2. Stage all changes:
```bash
git add .
```

3. Commit with a descriptive message (replace MESSAGE with actual description):
```bash
git commit -m "MESSAGE"
```

4. Push to GitHub:
```bash
git push origin main
```

## If push fails (first time or branch issues):
```bash
git push -u origin main --force
```

## Troubleshooting

If not authenticated, user needs to:
1. Run `git config --global user.email "you@example.com"`
2. Run `git config --global user.name "Your Name"`
3. Set up GitHub CLI or SSH keys for authentication
