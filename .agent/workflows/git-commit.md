---
description: Commit and push changes to GitHub
---

# Git Commit Workflow

Use this workflow when the user asks to commit changes to GitHub.

## Steps

// turbo-all

1. Check git status:
```powershell
& "C:\Program Files\Git\bin\git.exe" status
```

2. Stage all changes:
```powershell
& "C:\Program Files\Git\bin\git.exe" add .
```

3. Commit with message (replace MESSAGE with description from user or generate a descriptive one):
```powershell
& "C:\Program Files\Git\bin\git.exe" commit -m "MESSAGE"
```

4. Push to GitHub:
```powershell
& "C:\Program Files\Git\bin\git.exe" push origin main
```

## If push requires authentication:
User may see a browser popup for GitHub authentication. Wait for them to complete it.
