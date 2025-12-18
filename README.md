# 🧬 scRNA-seq Multi-Agent Analysis Platform

A comprehensive single-cell RNA-sequencing analysis platform powered by a multi-agent architecture. Upload raw sequencing data and receive publication-ready analysis through coordinated AI agents.

## ✨ Features

- **Multi-Agent Architecture**: 7 specialized agents coordinate analysis
- **Complete Pipeline**: QC → Filtering → Normalization → Clustering → Annotation
- **Publication-Ready**: Automated visualizations and reports
- **Real-Time Progress**: WebSocket-powered live updates
- **Biological Validation**: Control Agent ensures scientific rigor

## 🏗️ Architecture

| Agent | Role |
|-------|------|
| Program Manager | Pipeline orchestration & strategy |
| Execution Agent | Runs Scanpy analysis steps |
| Visualization Agent | Publication-quality plots |
| Control Agent | Biological validation |
| Audit Agent | Reproducibility tracking |
| Debug Agent | Error handling |

## 🚀 Quick Start

### Docker (Recommended)
```bash
docker-compose up --build
```

### Local Development

**Backend:**
```bash
cd backend
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

**Celery Worker:**
```bash
celery -A app.celery_worker worker --loglevel=info
```

**Frontend:**
```bash
cd frontend
npm install
npm run dev
```

## 📊 Pipeline Steps

1. **Quality Control** - Calculate %MT, gene counts
2. **Cell Filtering** - Remove low-quality cells
3. **Normalization** - LogNormalize/SCTransform
4. **HVG Selection** - Identify variable genes
5. **PCA** - Dimensionality reduction
6. **Clustering** - Leiden/Louvain
7. **UMAP** - 2D embedding
8. **Markers** - Differential expression
9. **Annotation** - CellTypist cell typing

## 🔧 Tech Stack

- **Backend**: FastAPI, Celery, Redis, Scanpy
- **Frontend**: React, Vite
- **Containerization**: Docker

## 📁 Project Structure

```
├── backend/
│   ├── app/
│   │   ├── agents/          # Multi-agent system
│   │   ├── api/             # FastAPI endpoints
│   │   ├── main.py          # App entry
│   │   └── tasks.py         # Celery tasks
│   └── requirements.txt
├── frontend/
│   └── src/
│       ├── components/      # React components
│       └── App.jsx
└── docker-compose.yml
```

## 📄 License

MIT License
