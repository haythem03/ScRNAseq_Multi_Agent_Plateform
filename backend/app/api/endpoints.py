"""
API Endpoints for the scRNA-seq platform.
Supports both async (Celery) and background-sync execution modes.
"""
from fastapi import APIRouter, File, UploadFile, HTTPException, WebSocket, WebSocketDisconnect, BackgroundTasks, Body, Query
from fastapi.responses import JSONResponse
import shutil
import os
import uuid
import json
from typing import Optional
from app.agents.program_manager import ProgramManager

router = APIRouter()

# Use local data directory for Windows compatibility
UPLOAD_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "uploads")
os.makedirs(UPLOAD_DIR, exist_ok=True)

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
