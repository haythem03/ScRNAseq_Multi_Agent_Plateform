"""
Program Manager Agent - Strategic lead for pipeline orchestration.
"""
from typing import Dict, Any, List, Optional
from .base_agent import BaseAgent
from .execution_agent import ExecutionAgent
from .visualization_agent import VisualizationAgent
from .control_agent import ControlAgent
from .audit_agent import AuditAgent

class ProgramManager(BaseAgent):
    """
    Orchestrates the complete scRNA-seq analysis pipeline.
    Coordinates all agents and makes strategic decisions.
    """
    
    # Default pipeline configuration
    DEFAULT_PIPELINE_CONFIG = {
        "qc": {"enabled": True},
        "filter": {
            "enabled": True,
            "min_genes": 200,
            "max_genes": 5000,
            "max_mito_pct": 20
        },
        "normalize": {"enabled": True, "method": "log_normalize"},
        "hvg": {"enabled": True, "n_top_genes": 2000},
        "pca": {"enabled": True, "n_comps": 50},
        "neighbors": {"enabled": True, "n_neighbors": 15},
        "cluster": {"enabled": True, "resolution": 1.0, "method": "leiden"},
        "umap": {"enabled": True},
        "markers": {"enabled": True, "n_genes": 25},
        "annotate": {"enabled": True, "method": "celltypist"}
    }
    
    PIPELINE_STEPS = [
        "qc", "filter", "normalize", "hvg", "pca", 
        "neighbors", "cluster", "umap", "markers", "annotate"
    ]
    
    def __init__(self):
        super().__init__(agent_id="PM-001", name="Program Manager")
        self.execution_agent = ExecutionAgent()
        self.visualization_agent = VisualizationAgent()
        self.control_agent = ControlAgent()
        self.audit_agent = AuditAgent()
        self.pipeline_state = {}
        self.current_step_index = 0
    
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Orchestrate the analysis pipeline."""
        task_type = payload.get("task_type")
        self.log(f"Received task: {task_type}")
        
        if task_type == "upload_and_qc":
            return self._handle_upload_and_qc(payload)
        elif task_type == "run_full_pipeline":
            return self._run_full_pipeline(payload)
        elif task_type == "run_step":
            return self._run_single_step(payload)
        elif task_type == "get_pipeline_status":
            return self._get_pipeline_status()
        elif task_type == "configure_pipeline":
            return self._configure_pipeline(payload.get("config", {}))
        else:
            return {"status": "error", "message": f"Unknown task type: {task_type}"}

    def _handle_upload_and_qc(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Handle initial upload and QC."""
        file_path = payload.get("file_path")
        self.log(f"Planning QC for file: {file_path}")
        
        # Start audit
        self.audit_agent.execute({
            "action": "start_run",
            "metadata": {"file_path": file_path}
        })
        
        # Execute QC
        qc_result = self.execution_agent.execute({
            "action": "load_and_qc",
            "file_path": file_path
        })
        
        if qc_result.get("status") != "success":
            return {
                "status": "error",
                "step": "qc",
                "file": file_path,
                "qc_results": qc_result
            }
        
        # Validate QC with Control Agent
        qc_data = qc_result.get("data", {})
        validation = self.control_agent.execute({
            "action": "validate_qc",
            "qc_metrics": {
                "n_cells": qc_data.get("n_cells"),
                "n_genes": qc_data.get("n_genes"),
                "pct_counts_mt": qc_data.get("qc_plots", {}).get("pct_counts_mt", []),
                "n_genes_by_counts": qc_data.get("qc_plots", {}).get("n_genes_by_counts", [])
            }
        })
        
        # Generate QC visualizations
        viz_result = self.visualization_agent.execute({
            "action": "generate_qc_plots",
            "qc_data": qc_data.get("qc_plots", {})
        })
        
        # Log to audit
        self.audit_agent.execute({
            "action": "log_step",
            "step_name": "qc",
            "parameters": {"file_path": file_path},
            "outputs": {"n_cells": qc_data.get("n_cells"), "n_genes": qc_data.get("n_genes")}
        })
        
        # Update pipeline state
        self.pipeline_state["qc"] = {
            "completed": True,
            "result": qc_result,
            "validation": validation,
            "plots": viz_result.get("plots", {})
        }
        self.current_step_index = 1
        
        return {
            "status": "completed",
            "step": "qc",
            "file": file_path,
            "steps": {
                "qc": {
                    "result": qc_result,
                    "validation": validation,
                    "plots": viz_result.get("plots", {})
                }
            },
            "recommendations": validation.get("recommendation", {})
        }

    def _run_full_pipeline(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Run the complete analysis pipeline."""
        file_path = payload.get("file_path")
        config = payload.get("config", self.DEFAULT_PIPELINE_CONFIG.copy())
        
        self.log("Starting full pipeline execution")
        results = {"steps": {}, "status": "running"}
        
        # Start audit
        self.audit_agent.execute({
            "action": "start_run",
            "metadata": {"file_path": file_path, "config": config}
        })
        
        # Step 1: Load and QC
        if config.get("qc", {}).get("enabled", True):
            qc_result = self.execution_agent.execute({
                "action": "load_and_qc",
                "file_path": file_path
            })
            if qc_result.get("status") != "success":
                return {"status": "error", "step": "qc", "error": qc_result}
            
            # Validate
            qc_data = qc_result.get("data", {})
            validation = self.control_agent.execute({
                "action": "validate_qc",
                "qc_metrics": {
                    "n_cells": qc_data.get("n_cells"),
                    "n_genes": qc_data.get("n_genes"),
                    "pct_counts_mt": qc_data.get("qc_plots", {}).get("pct_counts_mt", []),
                    "n_genes_by_counts": qc_data.get("qc_plots", {}).get("n_genes_by_counts", [])
                }
            })
            
            if validation.get("decision") == "STOP":
                return {"status": "stopped", "step": "qc", "reason": validation.get("issues")}
            
            # Generate plots
            viz_result = self.visualization_agent.execute({
                "action": "generate_qc_plots",
                "qc_data": qc_data.get("qc_plots", {})
            })
            
            results["steps"]["qc"] = {
                "result": qc_result,
                "validation": validation,
                "plots": viz_result.get("plots", {})
            }
        
        # Step 2: Filter
        if config.get("filter", {}).get("enabled", True):
            filter_params = config.get("filter", {})
            filter_result = self.execution_agent.execute({
                "action": "filter_cells",
                "params": filter_params
            })
            if filter_result.get("status") != "success":
                return {"status": "error", "step": "filter", "error": filter_result}
            
            # Validate filtering
            filter_data = filter_result.get("data", {})
            validation = self.control_agent.execute({
                "action": "validate_filtering",
                "before": filter_data.get("before", {}),
                "after": filter_data.get("after", {})
            })
            
            if validation.get("decision") == "STOP":
                return {"status": "stopped", "step": "filter", "reason": validation.get("issues")}
            
            results["steps"]["filter"] = {"result": filter_result, "validation": validation}
        
        # Step 3: Normalize
        if config.get("normalize", {}).get("enabled", True):
            norm_method = config.get("normalize", {}).get("method", "log_normalize")
            norm_result = self.execution_agent.execute({
                "action": "normalize",
                "method": norm_method
            })
            if norm_result.get("status") != "success":
                return {"status": "error", "step": "normalize", "error": norm_result}
            results["steps"]["normalize"] = {"result": norm_result}
        
        # Step 4: Find HVGs
        if config.get("hvg", {}).get("enabled", True):
            hvg_params = config.get("hvg", {})
            hvg_result = self.execution_agent.execute({
                "action": "find_hvg",
                "params": hvg_params
            })
            if hvg_result.get("status") != "success":
                return {"status": "error", "step": "hvg", "error": hvg_result}
            results["steps"]["hvg"] = {"result": hvg_result}
        
        # Step 5: PCA
        if config.get("pca", {}).get("enabled", True):
            n_comps = config.get("pca", {}).get("n_comps", 50)
            pca_result = self.execution_agent.execute({
                "action": "run_pca",
                "n_comps": n_comps
            })
            if pca_result.get("status") != "success":
                return {"status": "error", "step": "pca", "error": pca_result}
            results["steps"]["pca"] = {"result": pca_result}
        
        # Step 6: Neighbors
        if config.get("neighbors", {}).get("enabled", True):
            n_neighbors = config.get("neighbors", {}).get("n_neighbors", 15)
            neighbors_result = self.execution_agent.execute({
                "action": "run_neighbors",
                "n_neighbors": n_neighbors
            })
            if neighbors_result.get("status") != "success":
                return {"status": "error", "step": "neighbors", "error": neighbors_result}
            results["steps"]["neighbors"] = {"result": neighbors_result}
        
        # Step 7: Cluster
        if config.get("cluster", {}).get("enabled", True):
            cluster_config = config.get("cluster", {})
            cluster_result = self.execution_agent.execute({
                "action": "cluster",
                "resolution": cluster_config.get("resolution", 1.0),
                "method": cluster_config.get("method", "leiden")
            })
            if cluster_result.get("status") != "success":
                return {"status": "error", "step": "cluster", "error": cluster_result}
            
            # Validate clustering
            cluster_data = cluster_result.get("data", {})
            validation = self.control_agent.execute({
                "action": "validate_clustering",
                "cluster_data": {
                    "n_clusters": cluster_data.get("n_clusters"),
                    "n_cells": len(cluster_data.get("clusters", [])),
                    "cluster_sizes": list(cluster_data.get("cluster_sizes", {}).values())
                }
            })
            
            results["steps"]["cluster"] = {"result": cluster_result, "validation": validation}
        
        # Step 8: UMAP
        if config.get("umap", {}).get("enabled", True):
            umap_result = self.execution_agent.execute({"action": "run_umap"})
            if umap_result.get("status") != "success":
                return {"status": "error", "step": "umap", "error": umap_result}
            
            # Generate UMAP plot
            cluster_data = results.get("steps", {}).get("cluster", {}).get("result", {}).get("data", {})
            umap_data = umap_result.get("data", {})
            
            viz_result = self.visualization_agent.execute({
                "action": "generate_umap_plot",
                "embedding": umap_data.get("umap_coords", []),
                "clusters": cluster_data.get("clusters", [])
            })
            
            results["steps"]["umap"] = {"result": umap_result, "plot": viz_result.get("plot")}
        
        # Step 9: Markers
        if config.get("markers", {}).get("enabled", True):
            markers_params = config.get("markers", {})
            markers_result = self.execution_agent.execute({
                "action": "find_markers",
                "params": markers_params
            })
            if markers_result.get("status") != "success":
                return {"status": "error", "step": "markers", "error": markers_result}
            
            # Validate markers
            markers_data = markers_result.get("data", {}).get("markers", {})
            validation = self.control_agent.execute({
                "action": "validate_markers",
                "markers": {k: v.get("genes", []) for k, v in markers_data.items()}
            })
            
            results["steps"]["markers"] = {"result": markers_result, "validation": validation}
        
        # Step 10: Annotate
        if config.get("annotate", {}).get("enabled", True):
            annotate_method = config.get("annotate", {}).get("method", "celltypist")
            annotate_result = self.execution_agent.execute({
                "action": "annotate_cells",
                "method": annotate_method
            })
            if annotate_result.get("status") != "success":
                # Annotation failure is non-fatal
                results["steps"]["annotate"] = {"result": annotate_result, "warning": "Annotation failed"}
            else:
                results["steps"]["annotate"] = {"result": annotate_result}
        
        # End audit
        self.audit_agent.execute({
            "action": "end_run",
            "status": "completed",
            "summary": {"steps_completed": list(results["steps"].keys())}
        })
        
        results["status"] = "completed"
        self.log("Full pipeline completed successfully")
        
        return results

    def _run_single_step(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Run a single pipeline step."""
        step = payload.get("step")
        params = payload.get("params", {})
        
        action_map = {
            "qc": "load_and_qc",
            "filter": "filter_cells",
            "normalize": "normalize",
            "hvg": "find_hvg",
            "pca": "run_pca",
            "neighbors": "run_neighbors",
            "cluster": "cluster",
            "umap": "run_umap",
            "markers": "find_markers",
            "annotate": "annotate_cells"
        }
        
        if step not in action_map:
            return {"status": "error", "message": f"Unknown step: {step}"}
        
        result = self.execution_agent.execute({
            "action": action_map[step],
            **params
        })
        
        return {"status": "success", "step": step, "result": result}

    def _get_pipeline_status(self) -> Dict[str, Any]:
        """Get current pipeline status."""
        return {
            "status": "success",
            "pipeline_state": self.pipeline_state,
            "current_step_index": self.current_step_index,
            "steps": self.PIPELINE_STEPS,
            "data_summary": self.execution_agent.execute({"action": "get_data_summary"})
        }

    def _configure_pipeline(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Update pipeline configuration."""
        for step, step_config in config.items():
            if step in self.DEFAULT_PIPELINE_CONFIG:
                self.DEFAULT_PIPELINE_CONFIG[step].update(step_config)
        
        return {"status": "success", "config": self.DEFAULT_PIPELINE_CONFIG}
