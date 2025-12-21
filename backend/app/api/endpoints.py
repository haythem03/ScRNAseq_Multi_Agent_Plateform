"""
API Endpoints for the scRNA-seq platform.
Supports both async (Celery) and background-sync execution modes.
"""
from fastapi import APIRouter, File, UploadFile, HTTPException, WebSocket, WebSocketDisconnect, BackgroundTasks, Body, Query
from fastapi.responses import JSONResponse, FileResponse
import shutil
import os
import uuid
import json
from datetime import datetime
from typing import Optional
from app.agents.program_manager import ProgramManager

router = APIRouter()

# Use local data directory for Windows compatibility
UPLOAD_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "uploads")
CHECKPOINT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "checkpoints")
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

# Store results for background-sync mode
sync_results: dict = {}

# Check if Celery/Redis is available
USE_CELERY = False
try:
    from app.tasks import run_pipeline, run_full_pipeline
    from celery.result import AsyncResult
    import redis
    r = redis.Redis(host='localhost', port=6379, socket_connect_timeout=1)
    r.ping()
    USE_CELERY = True
    print("Celery/Redis available - using async task queue")
except Exception as e:
    print(f"Celery/Redis not available ({str(e)[:50]}...) - using background-sync mode")

def run_sync_analysis(task_id: str, payload: dict):
    """Worker function for BackgroundTasks."""
    try:
        sync_results[task_id]["status"] = "PROGRESS"
        pm = ProgramManager()
        print(f"Starting background execution for task: {task_id}")
        result = pm.execute(payload)
        sync_results[task_id] = {"status": "SUCCESS", "result": result}
        print(f"Task {task_id} completed successfully")
    except Exception as e:
        print(f"Error in background task {task_id}: {str(e)}")
        sync_results[task_id] = {"status": "FAILURE", "error": str(e)}

@router.post("/upload")
async def upload_file(background_tasks: BackgroundTasks, file: UploadFile = File(...)):
    """Upload a data file and start QC analysis."""
    if not file:
        raise HTTPException(status_code=400, detail="No file sent")
    
    # Validate file extension
    valid_extensions = ['.h5ad', '.csv', '.h5', '.loom', '.mtx', '.txt']
    file_ext = os.path.splitext(file.filename)[1].lower()
    if file_ext not in valid_extensions:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid file type. Supported: {', '.join(valid_extensions)}"
        )
    
    # Secure filename and save
    file_id = str(uuid.uuid4())
    filename = f"{file_id}_{file.filename}"
    file_path = os.path.join(UPLOAD_DIR, filename)
    
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    if USE_CELERY:
        task = run_pipeline.delay("upload_and_qc", {"file_path": file_path})
        task_id = task.id
    else:
        task_id = f"sync-{file_id}"
        sync_results[task_id] = {"status": "PENDING"}
        background_tasks.add_task(run_sync_analysis, task_id, {
            "task_type": "upload_and_qc",
            "file_path": file_path
        })
    
    return {
        "filename": file.filename,
        "file_id": file_id,
        "task_id": task_id,
        "message": "Upload successful. QC analysis started."
    }

@router.post("/pipeline/start")
async def start_full_pipeline(
    background_tasks: BackgroundTasks, 
    file_id: str = Query(..., description="The ID of the uploaded file"), 
    config: Optional[dict] = Body(None)
):
    """Start the full analysis pipeline for an uploaded file."""
    print(f"API: Start pipeline for file_id={file_id}")
    
    # Find the file
    if not os.path.exists(UPLOAD_DIR):
        raise HTTPException(status_code=404, detail="Upload directory missing")
        
    files = os.listdir(UPLOAD_DIR)
    matching = [f for f in files if f.startswith(file_id)]
    
    if not matching:
        print(f"API: File not found for ID {file_id}. Available prefixes: {[f.split('_')[0] for f in files[:5]]}")
        raise HTTPException(status_code=404, detail=f"File not found for ID: {file_id}")
    
    file_path = os.path.join(UPLOAD_DIR, matching[0])
    
    if USE_CELERY:
        task = run_full_pipeline.delay(file_path, config or {})
        task_id = task.id
    else:
        task_id = f"sync-pipeline-{file_id}"
        sync_results[task_id] = {"status": "PENDING"}
        background_tasks.add_task(run_sync_analysis, task_id, {
            "task_type": "run_full_pipeline",
            "file_path": file_path,
            "config": config or {}
        })
    
    return {
        "file_id": file_id,
        "task_id": task_id,
        "message": "Full pipeline started."
    }

@router.get("/status/{task_id}")
def get_status(task_id: str):
    """Get the status of a running task."""
    if task_id.startswith("sync-"):
        if task_id in sync_results:
            return {
                "task_id": task_id,
                "status": sync_results[task_id]["status"],
                "result": sync_results[task_id].get("result"),
                "error": sync_results[task_id].get("error")
            }
        return {"task_id": task_id, "status": "PENDING", "result": None}
    
    if not USE_CELERY:
        return {"task_id": task_id, "status": "UNKNOWN", "result": None}
    
    task_result = AsyncResult(task_id)
    result = {"task_id": task_id, "status": task_result.status, "result": None}
    
    if task_result.ready():
        try:
            result["result"] = task_result.result
        except Exception as e:
            result["error"] = str(e)
    elif task_result.state == 'PROGRESS':
        result["progress"] = task_result.info
    
    return JSONResponse(content=result)

@router.get("/pipeline/{task_id}/step/{step_name}")
def get_step_result(task_id: str, step_name: str):
    """Get the result of a specific pipeline step."""
    if task_id.startswith("sync-") and task_id in sync_results:
        result = sync_results[task_id].get("result")
        if result and "steps" in result:
            step_result = result["steps"].get(step_name)
            if step_result:
                return {"status": "success", "step": step_name, "data": step_result}
        return {"status": "not_found", "message": f"Step '{step_name}' not found"}
    
    if not USE_CELERY:
        return {"status": "not_found", "message": "Task not found"}
    
    task_result = AsyncResult(task_id)
    if not task_result.ready():
        return {"status": "pending", "message": "Task still running"}
    
    result = task_result.result
    if result and "steps" in result:
        step_result = result["steps"].get(step_name)
        if step_result:
            return {"status": "success", "step": step_name, "data": step_result}
    
    return {"status": "not_found", "message": f"Step '{step_name}' not found"}

@router.get("/pipeline/{task_id}/plots")
def get_plots(task_id: str):
    """Get all visualization plots from a pipeline run."""
    if task_id.startswith("sync-") and task_id in sync_results:
        result = sync_results[task_id].get("result")
        plots = {}
        if result and "steps" in result:
            for step_name, step_data in result["steps"].items():
                if "plots" in step_data:
                    plots[step_name] = step_data["plots"]
                elif "plot" in step_data:
                    plots[step_name] = step_data["plot"]
        if result and "plots" in result:
            plots["qc"] = result["plots"]
        return {"status": "success", "plots": plots}
    
    if not USE_CELERY:
        return {"status": "not_found", "message": "Task not found"}
    
    task_result = AsyncResult(task_id)
    if not task_result.ready():
        return {"status": "pending", "message": "Task still running"}
    
    result = task_result.result
    plots = {}
    if result and "steps" in result:
        for step_name, step_data in result["steps"].items():
            if "plots" in step_data:
                plots[step_name] = step_data["plots"]
            elif "plot" in step_data:
                plots[step_name] = step_data["plot"]
    
    if result and "plots" in result:
        plots["qc"] = result["plots"]
    
    return {"status": "success", "plots": plots}

@router.get("/pipeline/config/defaults")
def get_default_config():
    """Get default pipeline configuration."""
    return {
        "qc": {"enabled": True},
        "filter": {
            "enabled": True,
            "min_genes": 200,
            "max_genes": 5000,
            "max_mito_pct": 20,
            "min_cells_per_gene": 3
        },
        "normalize": {"enabled": True, "method": "log_normalize"},
        "hvg": {"enabled": True, "n_top_genes": 2000, "flavor": "seurat"},
        "pca": {"enabled": True, "n_comps": 50},
        "neighbors": {"enabled": True, "n_neighbors": 15},
        "cluster": {"enabled": True, "resolution": 1.0, "method": "leiden"},
        "umap": {"enabled": True},
        "markers": {"enabled": True, "n_genes": 25, "method": "wilcoxon"},
        "annotate": {"enabled": True, "method": "celltypist"}
    }

@router.get("/mode")
def get_execution_mode():
    """Get current execution mode (sync, async, or background-sync)."""
    return {
        "mode": "celery" if USE_CELERY else "background-sync",
        "celery_available": USE_CELERY
    }


# Visualization directory for PNG files
VISUALIZATION_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "visualizations")
os.makedirs(VISUALIZATION_DIR, exist_ok=True)

@router.get("/checkpoints")
def list_checkpoints():
    """List all available checkpoint files."""
    if not os.path.exists(CHECKPOINT_DIR):
        return {"checkpoints": []}
    
    files = []
    for filename in os.listdir(CHECKPOINT_DIR):
        if filename.endswith('.h5ad'):
            filepath = os.path.join(CHECKPOINT_DIR, filename)
            stat = os.stat(filepath)
            files.append({
                "name": filename,
                "size_mb": round(stat.st_size / (1024 * 1024), 2),
                "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
            })
    
    # Sort by modification time (newest first)
    files.sort(key=lambda x: x["modified"], reverse=True)
    return {"checkpoints": files}


@router.get("/visualizations")
def list_visualizations():
    """List all available visualization PNG files."""
    if not os.path.exists(VISUALIZATION_DIR):
        return {"files": []}
    
    files = []
    for filename in os.listdir(VISUALIZATION_DIR):
        if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.svg')):
            filepath = os.path.join(VISUALIZATION_DIR, filename)
            stat = os.stat(filepath)
            files.append({
                "name": filename,
                "size_kb": round(stat.st_size / 1024, 2),
                "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
            })
    
    # Sort by modification time (newest first)
    files.sort(key=lambda x: x["modified"], reverse=True)
    return {"files": files}


@router.get("/visualizations/{filename}")
def get_visualization(filename: str):
    """Serve a visualization image file."""
    # Sanitize filename to prevent directory traversal
    safe_filename = os.path.basename(filename)
    filepath = os.path.join(VISUALIZATION_DIR, safe_filename)
    
    if not os.path.exists(filepath):
        raise HTTPException(status_code=404, detail=f"Visualization file not found: {safe_filename}")
    
    # Determine media type
    ext = os.path.splitext(safe_filename)[1].lower()
    media_types = {'.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg', '.svg': 'image/svg+xml'}
    media_type = media_types.get(ext, 'application/octet-stream')
    
    return FileResponse(path=filepath, media_type=media_type)


@router.get("/checkpoints/{filename}")
def download_checkpoint(filename: str):
    """Download a specific checkpoint file."""
    # Sanitize filename to prevent directory traversal
    safe_filename = os.path.basename(filename)
    filepath = os.path.join(CHECKPOINT_DIR, safe_filename)
    
    if not os.path.exists(filepath):
        raise HTTPException(status_code=404, detail=f"Checkpoint file not found: {safe_filename}")
    
    return FileResponse(
        path=filepath,
        filename=safe_filename,
        media_type="application/octet-stream"
    )


@router.get("/checkpoints/{filename}/visualize")
def visualize_checkpoint(filename: str):
    """Generate visualizations from a checkpoint h5ad file using scanpy."""
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import base64
    from io import BytesIO
    
    safe_filename = os.path.basename(filename)
    filepath = os.path.join(CHECKPOINT_DIR, safe_filename)
    
    if not os.path.exists(filepath):
        raise HTTPException(status_code=404, detail=f"Checkpoint file not found: {safe_filename}")
    
    try:
        adata = sc.read_h5ad(filepath)
        plots = {}
        
        # Generate UMAP plot if coordinates exist
        if 'X_umap' in adata.obsm:
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_facecolor('#1e293b')
            ax.set_facecolor('#0f172a')
            
            # Color by clusters if available
            if 'clusters' in adata.obs.columns:
                import numpy as np
                clusters = adata.obs['clusters'].astype('category')
                n_clusters = len(clusters.cat.categories)
                colors = plt.cm.tab20(np.linspace(0, 1, max(n_clusters, 20)))
                
                for i, cat in enumerate(clusters.cat.categories):
                    mask = clusters == cat
                    ax.scatter(
                        adata.obsm['X_umap'][mask, 0],
                        adata.obsm['X_umap'][mask, 1],
                        c=[colors[i % 20]],
                        label=f'Cluster {cat}',
                        s=15,
                        alpha=0.7
                    )
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                         facecolor='#1e293b', labelcolor='#f8fafc')
            else:
                ax.scatter(
                    adata.obsm['X_umap'][:, 0],
                    adata.obsm['X_umap'][:, 1],
                    c='#6366f1',
                    s=15,
                    alpha=0.7
                )
            
            ax.set_xlabel('UMAP 1', color='#f8fafc')
            ax.set_ylabel('UMAP 2', color='#f8fafc')
            ax.set_title(f'UMAP - {safe_filename}', color='#f8fafc', fontsize=14)
            ax.tick_params(colors='#f8fafc')
            ax.spines['bottom'].set_color('#f8fafc')
            ax.spines['left'].set_color('#f8fafc')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            buf = BytesIO()
            fig.savefig(buf, format='png', dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            plots['umap'] = base64.b64encode(buf.read()).decode('utf-8')
        
        # Generate PCA plot if coordinates exist
        if 'X_pca' in adata.obsm:
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_facecolor('#1e293b')
            ax.set_facecolor('#0f172a')
            
            if 'clusters' in adata.obs.columns:
                import numpy as np
                clusters = adata.obs['clusters'].astype('category')
                n_clusters = len(clusters.cat.categories)
                colors = plt.cm.tab20(np.linspace(0, 1, max(n_clusters, 20)))
                
                for i, cat in enumerate(clusters.cat.categories):
                    mask = clusters == cat
                    ax.scatter(
                        adata.obsm['X_pca'][mask, 0],
                        adata.obsm['X_pca'][mask, 1],
                        c=[colors[i % 20]],
                        label=f'Cluster {cat}',
                        s=15,
                        alpha=0.7
                    )
            else:
                ax.scatter(
                    adata.obsm['X_pca'][:, 0],
                    adata.obsm['X_pca'][:, 1],
                    c='#ec4899',
                    s=15,
                    alpha=0.7
                )
            
            ax.set_xlabel('PC1', color='#f8fafc')
            ax.set_ylabel('PC2', color='#f8fafc')
            ax.set_title(f'PCA - {safe_filename}', color='#f8fafc', fontsize=14)
            ax.tick_params(colors='#f8fafc')
            ax.spines['bottom'].set_color('#f8fafc')
            ax.spines['left'].set_color('#f8fafc')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            buf = BytesIO()
            fig.savefig(buf, format='png', dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            plots['pca'] = base64.b64encode(buf.read()).decode('utf-8')
        
        # Generate QC violin plots if QC metrics exist
        qc_cols = [c for c in adata.obs.columns if c in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']]
        if qc_cols:
            fig, axes = plt.subplots(1, len(qc_cols), figsize=(5 * len(qc_cols), 5))
            fig.patch.set_facecolor('#1e293b')
            if len(qc_cols) == 1:
                axes = [axes]
            
            colors = ['#6366f1', '#ec4899', '#10b981']
            for idx, col in enumerate(qc_cols):
                ax = axes[idx]
                ax.set_facecolor('#0f172a')
                data = adata.obs[col].dropna().values
                if len(data) > 0:
                    violin = ax.violinplot(data, showmedians=True)
                    for pc in violin['bodies']:
                        pc.set_facecolor(colors[idx % 3])
                        pc.set_alpha(0.7)
                    for partname in ['cmedians', 'cmins', 'cmaxes', 'cbars']:
                        if partname in violin:
                            violin[partname].set_color('#f8fafc')
                ax.set_title(col.replace('_', ' ').title(), color='#f8fafc')
                ax.tick_params(colors='#f8fafc')
                ax.set_xticks([])
                ax.spines['bottom'].set_color('#f8fafc')
                ax.spines['left'].set_color('#f8fafc')
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            buf = BytesIO()
            fig.savefig(buf, format='png', dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            plots['qc_violin'] = base64.b64encode(buf.read()).decode('utf-8')
        
        # Get summary info
        info = {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "obs_columns": list(adata.obs.columns),
            "var_columns": list(adata.var.columns),
            "obsm_keys": list(adata.obsm.keys()),
            "layers": list(adata.layers.keys()) if adata.layers else [],
            "uns_keys": list(adata.uns.keys()) if adata.uns else []
        }
        
        if 'clusters' in adata.obs.columns:
            info["n_clusters"] = int(adata.obs['clusters'].nunique())
            info["cluster_sizes"] = adata.obs['clusters'].value_counts().to_dict()
        
        return {"status": "success", "filename": safe_filename, "info": info, "plots": plots}
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Error visualizing file: {str(e)}")
