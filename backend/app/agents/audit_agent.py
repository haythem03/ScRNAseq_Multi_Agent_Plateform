"""
Audit Agent - Ensures reproducibility, tracks parameters, and validates pipeline outputs.
"""
from typing import Dict, Any, List, Optional
import os
import json
import hashlib
from datetime import datetime
from .base_agent import BaseAgent

class AuditAgent(BaseAgent):
    """
    Tracks all pipeline parameters and outputs for reproducibility.
    Generates audit reports and validates output formats.
    """
    
    def __init__(self):
        super().__init__(agent_id="AUD-001", name="Audit Agent")
        self.audit_dir = os.path.join(os.getcwd(), 'data', 'audit')
        os.makedirs(self.audit_dir, exist_ok=True)
        self.current_run_id = None
        self.run_log = []
    
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        action = payload.get("action")
        
        if action == "start_run":
            return self._start_run(payload.get("metadata", {}))
        elif action == "log_step":
            return self._log_step(payload.get("step_name"), payload.get("parameters"), 
                                  payload.get("outputs"), payload.get("duration_ms"))
        elif action == "validate_output":
            return self._validate_output(payload.get("output_type"), payload.get("data"))
        elif action == "end_run":
            return self._end_run(payload.get("status"), payload.get("summary"))
        elif action == "get_run_report":
            return self._get_run_report(payload.get("run_id"))
        elif action == "compare_runs":
            return self._compare_runs(payload.get("run_id_1"), payload.get("run_id_2"))
        else:
            return {"status": "error", "message": f"Unknown action: {action}"}
    
    def _start_run(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """Initialize a new pipeline run for auditing."""
        self.current_run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_log = []
        
        run_info = {
            "run_id": self.current_run_id,
            "started_at": datetime.now().isoformat(),
            "metadata": metadata,
            "steps": []
        }
        
        # Save initial run file
        run_file = os.path.join(self.audit_dir, f"run_{self.current_run_id}.json")
        with open(run_file, 'w') as f:
            json.dump(run_info, f, indent=2)
        
        self.log(f"Started audit run: {self.current_run_id}")
        return {"status": "success", "run_id": self.current_run_id}
    
    def _log_step(self, step_name: str, parameters: Dict[str, Any], 
                  outputs: Dict[str, Any], duration_ms: Optional[int] = None) -> Dict[str, Any]:
        """Log a pipeline step execution."""
        if not self.current_run_id:
            return {"status": "error", "message": "No active run. Call start_run first."}
        
        step_entry = {
            "step_name": step_name,
            "timestamp": datetime.now().isoformat(),
            "parameters": parameters,
            "outputs_summary": self._summarize_outputs(outputs),
            "duration_ms": duration_ms,
            "parameter_hash": self._hash_params(parameters)
        }
        
        self.run_log.append(step_entry)
        
        # Update run file
        run_file = os.path.join(self.audit_dir, f"run_{self.current_run_id}.json")
        if os.path.exists(run_file):
            with open(run_file, 'r') as f:
                run_info = json.load(f)
            run_info['steps'].append(step_entry)
            with open(run_file, 'w') as f:
                json.dump(run_info, f, indent=2)
        
        self.log(f"Logged step: {step_name}")
        return {"status": "success", "step_logged": step_name}
    
    def _validate_output(self, output_type: str, data: Any) -> Dict[str, Any]:
        """Validate that output matches expected format."""
        validators = {
            "qc_metrics": self._validate_qc_output,
            "clustering": self._validate_clustering_output,
            "markers": self._validate_markers_output,
            "adata": self._validate_adata_output,
        }
        
        validator = validators.get(output_type)
        if not validator:
            return {"status": "warning", "message": f"No validator for output type: {output_type}"}
        
        try:
            is_valid, issues = validator(data)
            return {
                "status": "success" if is_valid else "invalid",
                "valid": is_valid,
                "issues": issues
            }
        except Exception as e:
            return {"status": "error", "message": str(e)}
    
    def _validate_qc_output(self, data: Dict[str, Any]) -> tuple:
        """Validate QC output format."""
        issues = []
        required_fields = ['n_cells', 'n_genes']
        
        for field in required_fields:
            if field not in data:
                issues.append(f"Missing required field: {field}")
        
        if 'n_cells' in data and not isinstance(data['n_cells'], int):
            issues.append("n_cells should be an integer")
        
        return len(issues) == 0, issues
    
    def _validate_clustering_output(self, data: Dict[str, Any]) -> tuple:
        """Validate clustering output format."""
        issues = []
        
        if 'clusters' not in data:
            issues.append("Missing 'clusters' field")
        elif not isinstance(data['clusters'], list):
            issues.append("'clusters' should be a list")
        
        if 'n_clusters' not in data:
            issues.append("Missing 'n_clusters' field")
        
        return len(issues) == 0, issues
    
    def _validate_markers_output(self, data: Dict[str, Any]) -> tuple:
        """Validate marker genes output format."""
        issues = []
        
        if not isinstance(data, dict):
            issues.append("Markers should be a dictionary keyed by cluster")
        
        return len(issues) == 0, issues
    
    def _validate_adata_output(self, data: Dict[str, Any]) -> tuple:
        """Validate AnnData object summary."""
        issues = []
        
        if 'path' in data and not os.path.exists(data['path']):
            issues.append(f"AnnData file not found: {data['path']}")
        
        return len(issues) == 0, issues
    
    def _end_run(self, status: str, summary: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Finalize the current run."""
        if not self.current_run_id:
            return {"status": "error", "message": "No active run"}
        
        run_file = os.path.join(self.audit_dir, f"run_{self.current_run_id}.json")
        if os.path.exists(run_file):
            with open(run_file, 'r') as f:
                run_info = json.load(f)
            
            run_info['ended_at'] = datetime.now().isoformat()
            run_info['final_status'] = status
            run_info['summary'] = summary or {}
            run_info['total_steps'] = len(self.run_log)
            
            # Calculate total duration
            if run_info.get('started_at'):
                start = datetime.fromisoformat(run_info['started_at'])
                end = datetime.fromisoformat(run_info['ended_at'])
                run_info['total_duration_seconds'] = (end - start).total_seconds()
            
            with open(run_file, 'w') as f:
                json.dump(run_info, f, indent=2)
        
        self.log(f"Ended run: {self.current_run_id} with status: {status}")
        completed_run_id = self.current_run_id
        self.current_run_id = None
        self.run_log = []
        
        return {"status": "success", "run_id": completed_run_id, "report_path": run_file}
    
    def _get_run_report(self, run_id: Optional[str] = None) -> Dict[str, Any]:
        """Retrieve a run report."""
        target_id = run_id or self.current_run_id
        if not target_id:
            return {"status": "error", "message": "No run ID specified"}
        
        run_file = os.path.join(self.audit_dir, f"run_{target_id}.json")
        if not os.path.exists(run_file):
            return {"status": "error", "message": f"Run not found: {target_id}"}
        
        with open(run_file, 'r') as f:
            return {"status": "success", "report": json.load(f)}
    
    def _compare_runs(self, run_id_1: str, run_id_2: str) -> Dict[str, Any]:
        """Compare two runs for reproducibility check."""
        report1 = self._get_run_report(run_id_1)
        report2 = self._get_run_report(run_id_2)
        
        if report1.get('status') != 'success' or report2.get('status') != 'success':
            return {"status": "error", "message": "Could not load one or both runs"}
        
        differences = []
        r1, r2 = report1['report'], report2['report']
        
        # Compare step counts
        if r1.get('total_steps') != r2.get('total_steps'):
            differences.append(f"Different step counts: {r1.get('total_steps')} vs {r2.get('total_steps')}")
        
        # Compare parameter hashes
        steps1 = {s['step_name']: s.get('parameter_hash') for s in r1.get('steps', [])}
        steps2 = {s['step_name']: s.get('parameter_hash') for s in r2.get('steps', [])}
        
        for step_name in set(steps1.keys()) | set(steps2.keys()):
            h1, h2 = steps1.get(step_name), steps2.get(step_name)
            if h1 != h2:
                differences.append(f"Parameter difference in step '{step_name}'")
        
        return {
            "status": "success",
            "identical": len(differences) == 0,
            "differences": differences
        }
    
    def _summarize_outputs(self, outputs: Dict[str, Any]) -> Dict[str, Any]:
        """Create a summary of outputs without full data."""
        summary = {}
        for key, value in outputs.items():
            if isinstance(value, (list, tuple)):
                summary[key] = f"<list: {len(value)} items>"
            elif isinstance(value, dict):
                summary[key] = f"<dict: {len(value)} keys>"
            elif isinstance(value, str) and len(value) > 100:
                summary[key] = f"<string: {len(value)} chars>"
            else:
                summary[key] = value
        return summary
    
    def _hash_params(self, parameters: Dict[str, Any]) -> str:
        """Create a hash of parameters for reproducibility tracking."""
        param_str = json.dumps(parameters, sort_keys=True, default=str)
        return hashlib.md5(param_str.encode()).hexdigest()[:12]
