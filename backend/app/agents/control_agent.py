"""
Control Agent - Biology Expert for validation and quality control at each pipeline step.
"""
from typing import Dict, Any, List, Optional
from .base_agent import BaseAgent

class ControlAgent(BaseAgent):
    """
    Validates scientific rigor at each analysis step.
    Can STOP the pipeline if quality issues are detected.
    """
    
    # Configurable thresholds
    DEFAULT_THRESHOLDS = {
        'min_cells': 100,
        'max_mito_pct': 20.0,
        'min_mito_pct': 0.0,
        'min_genes_per_cell': 200,
        'max_genes_per_cell': 8000,
        'min_hvgs': 500,
        'max_clusters': 30,  # For <10K cells
        'min_clusters': 2,
        'max_doublet_rate': 0.15,
    }
    
    def __init__(self):
        super().__init__(agent_id="CTRL-001", name="Control Agent")
        self.thresholds = self.DEFAULT_THRESHOLDS.copy()
        self.validation_log = []
    
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        action = payload.get("action")
        
        if action == "validate_qc":
            return self._validate_qc(payload.get("qc_metrics"))
        elif action == "validate_filtering":
            return self._validate_filtering(payload.get("before"), payload.get("after"))
        elif action == "validate_clustering":
            return self._validate_clustering(payload.get("cluster_data"))
        elif action == "validate_markers":
            return self._validate_markers(payload.get("markers"))
        elif action == "validate_annotation":
            return self._validate_annotation(payload.get("annotations"))
        elif action == "set_thresholds":
            return self._set_thresholds(payload.get("thresholds"))
        elif action == "get_validation_log":
            return {"status": "success", "log": self.validation_log}
        else:
            return {"status": "error", "message": f"Unknown action: {action}"}
    
    def _log_validation(self, step: str, result: str, details: str):
        """Log validation decision."""
        entry = {"step": step, "result": result, "details": details}
        self.validation_log.append(entry)
        self.log(f"[{step}] {result}: {details}")
    
    def _validate_qc(self, qc_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate QC metrics against biological expectations.
        """
        issues = []
        warnings = []
        
        n_cells = qc_metrics.get('n_cells', 0)
        n_genes = qc_metrics.get('n_genes', 0)
        pct_mito = qc_metrics.get('pct_counts_mt', [])
        n_genes_per_cell = qc_metrics.get('n_genes_by_counts', [])
        
        # Check cell count
        if n_cells < self.thresholds['min_cells']:
            issues.append(f"Too few cells: {n_cells} < {self.thresholds['min_cells']}. "
                         "Analysis may not be statistically meaningful.")
        
        # Check mitochondrial percentage distribution
        if pct_mito:
            median_mito = self._median(pct_mito)
            max_mito = max(pct_mito)
            
            if median_mito > self.thresholds['max_mito_pct']:
                issues.append(f"High median %MT ({median_mito:.1f}%). "
                             "May indicate cell damage or apoptosis.")
            elif median_mito > 10:
                warnings.append(f"Elevated median %MT ({median_mito:.1f}%). "
                               "Consider tissue type - some tissues have higher baseline.")
        
        # Check genes per cell distribution
        if n_genes_per_cell:
            median_genes = self._median(n_genes_per_cell)
            if median_genes < self.thresholds['min_genes_per_cell']:
                warnings.append(f"Low median genes/cell ({median_genes:.0f}). "
                               "Data may be low quality or highly sparse.")
        
        # Determine decision
        if issues:
            decision = "STOP"
            self._log_validation("QC", "STOP", "; ".join(issues))
        elif warnings:
            decision = "PROCEED_WITH_CAUTION"
            self._log_validation("QC", "WARNING", "; ".join(warnings))
        else:
            decision = "PROCEED"
            self._log_validation("QC", "PASS", "QC metrics within acceptable ranges")
        
        return {
            "status": "success",
            "decision": decision,
            "issues": issues,
            "warnings": warnings,
            "recommendation": self._get_qc_recommendation(qc_metrics)
        }
    
    def _validate_filtering(self, before: Dict[str, int], after: Dict[str, int]) -> Dict[str, Any]:
        """Validate that filtering retained sufficient data."""
        issues = []
        warnings = []
        
        cells_before = before.get('n_cells', 1)
        cells_after = after.get('n_cells', 0)
        retention_rate = cells_after / cells_before if cells_before > 0 else 0
        
        if cells_after < self.thresholds['min_cells']:
            issues.append(f"Too few cells remaining after filtering: {cells_after}")
        
        if retention_rate < 0.5:
            warnings.append(f"Filtered out {(1-retention_rate)*100:.1f}% of cells. "
                           "Consider relaxing thresholds.")
        
        if issues:
            decision = "STOP"
            self._log_validation("Filtering", "STOP", "; ".join(issues))
        elif warnings:
            decision = "PROCEED_WITH_CAUTION"
            self._log_validation("Filtering", "WARNING", "; ".join(warnings))
        else:
            decision = "PROCEED"
            self._log_validation("Filtering", "PASS", 
                               f"Retained {retention_rate*100:.1f}% of cells ({cells_after})")
        
        return {
            "status": "success",
            "decision": decision,
            "issues": issues,
            "warnings": warnings,
            "cells_retained": cells_after,
            "retention_rate": retention_rate
        }
    
    def _validate_clustering(self, cluster_data: Dict[str, Any]) -> Dict[str, Any]:
        """Validate clustering results."""
        issues = []
        warnings = []
        
        n_clusters = cluster_data.get('n_clusters', 0)
        n_cells = cluster_data.get('n_cells', 0)
        cluster_sizes = cluster_data.get('cluster_sizes', [])
        
        # Check for over/under-clustering
        if n_clusters < self.thresholds['min_clusters']:
            issues.append(f"Under-clustering detected: only {n_clusters} clusters. "
                         "Consider lowering resolution.")
        
        cells_per_cluster_threshold = 10000  # For large datasets, more clusters expected
        max_expected = self.thresholds['max_clusters']
        if n_cells < cells_per_cluster_threshold and n_clusters > max_expected:
            warnings.append(f"Potential over-clustering: {n_clusters} clusters for {n_cells} cells. "
                           "Consider increasing resolution.")
        
        # Check for imbalanced clusters
        if cluster_sizes:
            min_size = min(cluster_sizes)
            max_size = max(cluster_sizes)
            if min_size < 10:
                warnings.append(f"Very small cluster detected ({min_size} cells). "
                               "May represent doublets or technical artifacts.")
        
        if issues:
            decision = "STOP"
            self._log_validation("Clustering", "STOP", "; ".join(issues))
        elif warnings:
            decision = "PROCEED_WITH_CAUTION"
            self._log_validation("Clustering", "WARNING", "; ".join(warnings))
        else:
            decision = "PROCEED"
            self._log_validation("Clustering", "PASS", 
                               f"{n_clusters} clusters identified")
        
        return {
            "status": "success",
            "decision": decision,
            "issues": issues,
            "warnings": warnings
        }
    
    def _validate_markers(self, markers: Dict[str, Any]) -> Dict[str, Any]:
        """Validate marker gene identification."""
        issues = []
        warnings = []
        
        # Check if markers were found
        for cluster_id, cluster_markers in markers.items():
            if not cluster_markers or len(cluster_markers) == 0:
                warnings.append(f"No significant markers for cluster {cluster_id}")
            else:
                # Check for housekeeping/ribosomal dominance
                top_genes = cluster_markers[:5] if isinstance(cluster_markers, list) else []
                ribo_count = sum(1 for g in top_genes if g.startswith(('RPS', 'RPL', 'MT-', 'mt-')))
                if ribo_count >= 3:
                    warnings.append(f"Cluster {cluster_id} markers dominated by ribosomal/MT genes")
        
        if issues:
            decision = "STOP"
        elif warnings:
            decision = "PROCEED_WITH_CAUTION"
            self._log_validation("Markers", "WARNING", "; ".join(warnings))
        else:
            decision = "PROCEED"
            self._log_validation("Markers", "PASS", "Marker genes identified for all clusters")
        
        return {
            "status": "success",
            "decision": decision,
            "issues": issues,
            "warnings": warnings
        }
    
    def _validate_annotation(self, annotations: Dict[str, Any]) -> Dict[str, Any]:
        """Validate cell type annotations."""
        issues = []
        warnings = []
        
        confidence_scores = annotations.get('confidence_scores', {})
        cell_types = annotations.get('cell_types', {})
        
        # Check for low confidence annotations
        low_confidence = [ct for ct, score in confidence_scores.items() if score < 0.5]
        if low_confidence:
            warnings.append(f"Low confidence annotations for: {', '.join(low_confidence)}")
        
        # Check for unexpected cell types (would need tissue context)
        if not cell_types:
            warnings.append("No cell types annotated")
        
        if issues:
            decision = "STOP"
        elif warnings:
            decision = "PROCEED_WITH_CAUTION"
            self._log_validation("Annotation", "WARNING", "; ".join(warnings))
        else:
            decision = "PROCEED"
            self._log_validation("Annotation", "PASS", f"{len(cell_types)} cell types identified")
        
        return {
            "status": "success",
            "decision": decision,
            "issues": issues,
            "warnings": warnings
        }
    
    def _set_thresholds(self, thresholds: Dict[str, Any]) -> Dict[str, Any]:
        """Update validation thresholds."""
        for key, value in thresholds.items():
            if key in self.thresholds:
                self.thresholds[key] = value
                self.log(f"Updated threshold {key} = {value}")
        return {"status": "success", "thresholds": self.thresholds}
    
    def _get_qc_recommendation(self, qc_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """Suggest filtering thresholds based on data distribution."""
        pct_mito = qc_metrics.get('pct_counts_mt', [])
        n_genes = qc_metrics.get('n_genes_by_counts', [])
        
        recommendations = {}
        
        if pct_mito:
            # Suggest threshold at 95th percentile or 20%, whichever is lower
            p95_mito = self._percentile(pct_mito, 95)
            recommendations['max_mito_pct'] = min(p95_mito, 20.0)
        
        if n_genes:
            # Filter bottom 5% and top 1%
            recommendations['min_genes'] = max(self._percentile(n_genes, 5), 200)
            recommendations['max_genes'] = self._percentile(n_genes, 99)
        
        return recommendations
    
    def _median(self, values: List[float]) -> float:
        """Calculate median of a list."""
        sorted_vals = sorted(values)
        n = len(sorted_vals)
        if n == 0:
            return 0
        mid = n // 2
        if n % 2 == 0:
            return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2
        return sorted_vals[mid]
    
    def _percentile(self, values: List[float], p: float) -> float:
        """Calculate percentile of a list."""
        sorted_vals = sorted(values)
        n = len(sorted_vals)
        if n == 0:
            return 0
        idx = int(n * p / 100)
        return sorted_vals[min(idx, n - 1)]
