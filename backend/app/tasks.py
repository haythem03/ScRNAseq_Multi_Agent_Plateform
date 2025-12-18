"""
Celery tasks for the scRNA-seq pipeline.
"""
from app.celery_worker import celery_app
from app.agents.program_manager import ProgramManager
from typing import Dict, Any

@celery_app.task(bind=True, name="run_pipeline")
def run_pipeline(self, task_type: str, payload: Dict[str, Any]):
    """Run a pipeline task (e.g., upload_and_qc)."""
    pm = ProgramManager()
    payload["task_type"] = task_type
    
    # Update task state for progress tracking
    self.update_state(state='PROGRESS', meta={'step': 'starting', 'message': f'Running {task_type}'})
    
    result = pm.execute(payload)
    return result

@celery_app.task(bind=True, name="run_full_pipeline")
def run_full_pipeline(self, file_path: str, config: Dict[str, Any]):
    """Run the complete analysis pipeline."""
    pm = ProgramManager()
    
    steps = ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster', 'umap', 'markers', 'annotate']
    total_steps = len(steps)
    
    # Update progress
    self.update_state(state='PROGRESS', meta={
        'current_step': 0,
        'total_steps': total_steps,
        'step_name': 'starting',
        'message': 'Initializing pipeline'
    })
    
    result = pm.execute({
        "task_type": "run_full_pipeline",
        "file_path": file_path,
        "config": config
    })
    
    return result

@celery_app.task(bind=True, name="run_single_step")
def run_single_step(self, step: str, params: Dict[str, Any]):
    """Run a single pipeline step."""
    pm = ProgramManager()
    
    self.update_state(state='PROGRESS', meta={'step': step, 'message': f'Running {step}'})
    
    result = pm.execute({
        "task_type": "run_step",
        "step": step,
        "params": params
    })
    
    return result
