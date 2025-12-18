"""
API Endpoints for the scRNA-seq platform.
"""
from fastapi import APIRouter, File, UploadFile, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.responses import JSONResponse
import shutil
import os
import uuid
import json
from typing import Optional
from app.tasks import run_pipeline, run_full_pipeline
from celery.result import AsyncResult

router = APIRouter()

UPLOAD_DIR = "/data/uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)

# Store active WebSocket connections
active_connections: dict = {}

@router.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    """Upload a data file and start QC analysis."""
    if not file:
        raise HTTPException(status_code=400, detail="No file sent")
    
    # Validate file extension
    valid_extensions = ['.h5ad', '.csv', '.h5', '.loom', '.mtx']
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
    
    # Trigger QC Task
    task = run_pipeline.delay("upload_and_qc", {"file_path": file_path})
    
    return {
        "filename": file.filename,
        "file_id": file_id,
        "task_id": task.id,
        "message": "Upload successful. QC analysis started."
    }

@router.post("/pipeline/start")
async def start_full_pipeline(file_id: str, config: Optional[dict] = None):
    """Start the full analysis pipeline for an uploaded file."""
    # Find the file
    files = os.listdir(UPLOAD_DIR)
    matching = [f for f in files if f.startswith(file_id)]
    
    if not matching:
        raise HTTPException(status_code=404, detail="File not found")
    
    file_path = os.path.join(UPLOAD_DIR, matching[0])
    
    # Start full pipeline task
    task = run_full_pipeline.delay(file_path, config or {})
    
    return {
        "file_id": file_id,
        "task_id": task.id,
        "message": "Full pipeline started."
    }

@router.get("/status/{task_id}")
def get_status(task_id: str):
    """Get the status of a running task."""
    task_result = AsyncResult(task_id)
    
    result = {
        "task_id": task_id,
        "status": task_result.status,
        "result": None
    }
    
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
    
    # Also check top-level plots
    if result and "plots" in result:
        plots["qc"] = result["plots"]
    
    return {"status": "success", "plots": plots}

@router.post("/pipeline/{task_id}/configure")
async def configure_pipeline(task_id: str, config: dict):
    """Update pipeline configuration (for future steps)."""
    # This would update the configuration for a running pipeline
    # For now, return acknowledgment
    return {"status": "success", "message": "Configuration updated", "config": config}

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

@router.websocket("/ws/pipeline/{task_id}")
async def websocket_pipeline_updates(websocket: WebSocket, task_id: str):
    """WebSocket endpoint for real-time pipeline updates."""
    await websocket.accept()
    active_connections[task_id] = websocket
    
    try:
        while True:
            # Check task status periodically
            task_result = AsyncResult(task_id)
            
            status_update = {
                "task_id": task_id,
                "status": task_result.status,
            }
            
            if task_result.state == 'PROGRESS':
                status_update["progress"] = task_result.info
            elif task_result.ready():
                status_update["result"] = task_result.result
                await websocket.send_json(status_update)
                break
            
            await websocket.send_json(status_update)
            
            # Wait for client message or timeout
            try:
                data = await websocket.receive_text()
                if data == "close":
                    break
            except:
                pass
                
    except WebSocketDisconnect:
        pass
    finally:
        if task_id in active_connections:
            del active_connections[task_id]

@router.get("/files")
def list_uploaded_files():
    """List all uploaded files."""
    if not os.path.exists(UPLOAD_DIR):
        return {"files": []}
    
    files = []
    for filename in os.listdir(UPLOAD_DIR):
        file_path = os.path.join(UPLOAD_DIR, filename)
        if os.path.isfile(file_path):
            parts = filename.split('_', 1)
            files.append({
                "file_id": parts[0] if len(parts) > 1 else filename,
                "filename": parts[1] if len(parts) > 1 else filename,
                "size_bytes": os.path.getsize(file_path)
            })
    
    return {"files": files}

@router.delete("/files/{file_id}")
def delete_file(file_id: str):
    """Delete an uploaded file."""
    files = os.listdir(UPLOAD_DIR)
    matching = [f for f in files if f.startswith(file_id)]
    
    if not matching:
        raise HTTPException(status_code=404, detail="File not found")
    
    file_path = os.path.join(UPLOAD_DIR, matching[0])
    os.remove(file_path)
    
    return {"status": "success", "message": f"File {file_id} deleted"}
