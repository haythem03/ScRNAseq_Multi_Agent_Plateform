"""
Execution Agent - Implements all scRNA-seq analysis steps.
"""
from typing import Dict, Any, Optional, List
import scanpy as sc
import pandas as pd
import numpy as np
import os
from .base_agent import BaseAgent

class ExecutionAgent(BaseAgent):
    """
    Executes all analysis steps as directed by the Program Manager.
    Maintains state of the current AnnData object throughout the pipeline.
    """
    
    def __init__(self):
        super().__init__(agent_id="EXE-001", name="Execution Agent")
        self.adata = None
        self.checkpoint_dir = os.path.join(os.getcwd(), 'data', 'checkpoints')
        os.makedirs(self.checkpoint_dir, exist_ok=True)
    
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        action = payload.get("action")
        
        action_map = {
            "load_and_qc": lambda: self._load_and_qc(payload.get("file_path")),
            "filter_cells": lambda: self._filter_cells(payload.get("params", {})),
            "normalize": lambda: self._normalize(payload.get("method", "log_normalize")),
            "find_hvg": lambda: self._find_hvg(payload.get("params", {})),
            "run_pca": lambda: self._run_pca(payload.get("n_comps", 50)),
            "run_neighbors": lambda: self._run_neighbors(payload.get("n_neighbors", 15)),
            "cluster": lambda: self._cluster(payload.get("resolution", 1.0), payload.get("method", "leiden")),
            "run_umap": lambda: self._run_umap(),
            "find_markers": lambda: self._find_markers(payload.get("params", {})),
            "annotate_cells": lambda: self._annotate_cells(payload.get("method", "celltypist")),
            "save_checkpoint": lambda: self._save_checkpoint(payload.get("name", "checkpoint")),
            "load_checkpoint": lambda: self._load_checkpoint(payload.get("name")),
            "get_data_summary": lambda: self._get_data_summary(),
        }
        
        if action in action_map:
            return action_map[action]()
        else:
            return {"status": "error", "message": f"Unknown action: {action}"}

    def _load_and_qc(self, file_path: str) -> Dict[str, Any]:
        """Load data and calculate QC metrics."""
        self.log(f"Loading file from {file_path}")
        
        try:
            if not os.path.exists(file_path):
                return {"status": "error", "message": "File not found"}

            # Load based on file type
            if file_path.endswith('.h5ad'):
                self.adata = sc.read_h5ad(file_path)
            elif file_path.endswith('.csv'):
                # Try to determine if genes are rows or columns
                df_preview = pd.read_csv(file_path, nrows=5, index_col=0)
                # Heuristic: if many numeric columns, genes are likely columns
                self.adata = sc.read_csv(file_path)
                # Transpose if needed (genes should be in var)
                if self.adata.n_vars < self.adata.n_obs:
                    self.log("Transposing data: detected genes as observations")
                    self.adata = self.adata.T
            elif file_path.endswith('.loom'):
                self.adata = sc.read_loom(file_path)
            elif file_path.endswith('.h5'):
                self.adata = sc.read_10x_h5(file_path)
            elif os.path.isdir(file_path):
                # 10X directory format
                self.adata = sc.read_10x_mtx(file_path)
            else:
                self.adata = sc.read(file_path)

            self.adata.var_names_make_unique()
            self.adata.obs_names_make_unique()
            
            # Detect mitochondrial genes (human: MT-, mouse: mt-)
            self.adata.var['mt'] = (
                self.adata.var_names.str.startswith('MT-') | 
                self.adata.var_names.str.startswith('mt-')
            )
            
            # Detect ribosomal genes
            self.adata.var['ribo'] = (
                self.adata.var_names.str.startswith('RPS') | 
                self.adata.var_names.str.startswith('RPL') |
                self.adata.var_names.str.startswith('Rps') |
                self.adata.var_names.str.startswith('Rpl')
            )
            
            # Calculate QC metrics
            sc.pp.calculate_qc_metrics(
                self.adata, 
                qc_vars=['mt', 'ribo'], 
                percent_top=None, 
                log1p=False, 
                inplace=True
            )
            
            # Save initial checkpoint
            self._save_checkpoint("initial_load")
            
            results = {
                "n_cells": int(self.adata.n_obs),
                "n_genes": int(self.adata.n_vars),
                "obs_head": self.adata.obs.head().to_dict(),
                "qc_plots": {
                    "n_genes_by_counts": self.adata.obs['n_genes_by_counts'].tolist(),
                    "total_counts": self.adata.obs['total_counts'].tolist(),
                    "pct_counts_mt": self.adata.obs['pct_counts_mt'].tolist()
                },
                "qc_summary": {
                    "median_genes": float(self.adata.obs['n_genes_by_counts'].median()),
                    "median_counts": float(self.adata.obs['total_counts'].median()),
                    "median_mito_pct": float(self.adata.obs['pct_counts_mt'].median()),
                    "mt_genes_detected": int(self.adata.var['mt'].sum()),
                }
            }
            self.log("QC execution successful")
            return {"status": "success", "data": results}

        except Exception as e:
            self.log(f"Error during QC: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _filter_cells(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Filter cells based on QC thresholds."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded. Run load_and_qc first."}
        
        try:
            before_cells = self.adata.n_obs
            before_genes = self.adata.n_vars
            
            # Extract filter parameters with defaults
            min_genes = params.get('min_genes', 200)
            max_genes = params.get('max_genes', None)
            min_counts = params.get('min_counts', None)
            max_counts = params.get('max_counts', None)
            max_mito_pct = params.get('max_mito_pct', 20)
            min_cells = params.get('min_cells_per_gene', 3)
            
            # Filter cells by gene count
            sc.pp.filter_cells(self.adata, min_genes=min_genes)
            if max_genes:
                self.adata = self.adata[self.adata.obs['n_genes_by_counts'] < max_genes, :]
            
            # Filter by counts
            if min_counts:
                self.adata = self.adata[self.adata.obs['total_counts'] > min_counts, :]
            if max_counts:
                self.adata = self.adata[self.adata.obs['total_counts'] < max_counts, :]
            
            # Filter by mitochondrial percentage
            if max_mito_pct is not None:
                self.adata = self.adata[self.adata.obs['pct_counts_mt'] < max_mito_pct, :]
            
            # Filter genes by cell count
            sc.pp.filter_genes(self.adata, min_cells=min_cells)
            
            self._save_checkpoint("after_filtering")
            
            results = {
                "before": {"n_cells": before_cells, "n_genes": before_genes},
                "after": {"n_cells": int(self.adata.n_obs), "n_genes": int(self.adata.n_vars)},
                "params_used": {
                    "min_genes": min_genes,
                    "max_genes": max_genes,
                    "max_mito_pct": max_mito_pct,
                    "min_cells_per_gene": min_cells
                }
            }
            
            self.log(f"Filtered: {before_cells} → {self.adata.n_obs} cells, {before_genes} → {self.adata.n_vars} genes")
            return {"status": "success", "data": results}
            
        except Exception as e:
            self.log(f"Error during filtering: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _normalize(self, method: str = "log_normalize") -> Dict[str, Any]:
        """Normalize the data."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            # Store raw counts
            self.adata.layers['counts'] = self.adata.X.copy()
            
            if method == "log_normalize":
                # Standard log normalization
                sc.pp.normalize_total(self.adata, target_sum=1e4)
                sc.pp.log1p(self.adata)
                self.log("Applied log normalization (normalize_total + log1p)")
                
            elif method == "scran":
                # scran pooling-based normalization (requires R via rpy2 or pre-computed factors)
                # Fallback to standard normalization if scran not available
                try:
                    sc.pp.normalize_total(self.adata, target_sum=1e4)
                    sc.pp.log1p(self.adata)
                    self.log("Applied standard normalization (scran fallback)")
                except:
                    sc.pp.normalize_total(self.adata, target_sum=1e4)
                    sc.pp.log1p(self.adata)
                    
            else:
                return {"status": "error", "message": f"Unknown normalization method: {method}"}
            
            self._save_checkpoint("after_normalization")
            
            return {
                "status": "success", 
                "data": {
                    "method": method,
                    "n_cells": int(self.adata.n_obs),
                    "n_genes": int(self.adata.n_vars)
                }
            }
            
        except Exception as e:
            self.log(f"Error during normalization: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _find_hvg(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Find highly variable genes."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            n_top_genes = params.get('n_top_genes', 2000)
            flavor = params.get('flavor', 'seurat')
            
            sc.pp.highly_variable_genes(
                self.adata, 
                n_top_genes=n_top_genes,
                flavor=flavor,
                subset=False  # Keep all genes but mark HVGs
            )
            
            hvg_count = int(self.adata.var['highly_variable'].sum())
            top_hvgs = self.adata.var[self.adata.var['highly_variable']].index.tolist()[:50]
            
            self._save_checkpoint("after_hvg")
            
            self.log(f"Identified {hvg_count} highly variable genes")
            return {
                "status": "success",
                "data": {
                    "n_hvgs": hvg_count,
                    "top_hvgs": top_hvgs,
                    "flavor": flavor
                }
            }
            
        except Exception as e:
            self.log(f"Error finding HVGs: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _run_pca(self, n_comps: int = 50) -> Dict[str, Any]:
        """Run PCA dimensionality reduction."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            # Subset to HVGs if available
            if 'highly_variable' in self.adata.var.columns:
                sc.pp.scale(self.adata, max_value=10)
                sc.tl.pca(self.adata, n_comps=min(n_comps, self.adata.n_vars - 1), use_highly_variable=True)
            else:
                sc.pp.scale(self.adata, max_value=10)
                sc.tl.pca(self.adata, n_comps=min(n_comps, self.adata.n_vars - 1))
            
            # Get variance ratios for elbow plot
            variance_ratio = self.adata.uns['pca']['variance_ratio'].tolist()
            
            self._save_checkpoint("after_pca")
            
            self.log(f"PCA completed with {n_comps} components")
            return {
                "status": "success",
                "data": {
                    "n_components": n_comps,
                    "variance_ratio": variance_ratio,
                    "cumulative_variance": np.cumsum(variance_ratio).tolist()
                }
            }
            
        except Exception as e:
            self.log(f"Error running PCA: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _run_neighbors(self, n_neighbors: int = 15) -> Dict[str, Any]:
        """Compute neighborhood graph."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            n_pcs = min(40, self.adata.obsm['X_pca'].shape[1] if 'X_pca' in self.adata.obsm else 40)
            sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            
            self.log(f"Computed neighbors with n_neighbors={n_neighbors}, n_pcs={n_pcs}")
            return {
                "status": "success",
                "data": {"n_neighbors": n_neighbors, "n_pcs": n_pcs}
            }
            
        except Exception as e:
            self.log(f"Error computing neighbors: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _cluster(self, resolution: float = 1.0, method: str = "leiden") -> Dict[str, Any]:
        """Perform clustering."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            if method == "leiden":
                sc.tl.leiden(self.adata, resolution=resolution, key_added='clusters')
            elif method == "louvain":
                sc.tl.louvain(self.adata, resolution=resolution, key_added='clusters')
            else:
                return {"status": "error", "message": f"Unknown clustering method: {method}"}
            
            cluster_counts = self.adata.obs['clusters'].value_counts().to_dict()
            n_clusters = len(cluster_counts)
            
            self._save_checkpoint("after_clustering")
            
            self.log(f"Clustering complete: {n_clusters} clusters found")
            return {
                "status": "success",
                "data": {
                    "method": method,
                    "resolution": resolution,
                    "n_clusters": n_clusters,
                    "cluster_sizes": cluster_counts,
                    "clusters": self.adata.obs['clusters'].tolist()
                }
            }
            
        except Exception as e:
            self.log(f"Error during clustering: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _run_umap(self) -> Dict[str, Any]:
        """Run UMAP embedding."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            sc.tl.umap(self.adata)
            
            umap_coords = self.adata.obsm['X_umap'].tolist()
            
            self._save_checkpoint("after_umap")
            
            self.log("UMAP embedding complete")
            return {
                "status": "success",
                "data": {
                    "umap_coords": umap_coords,
                    "n_cells": len(umap_coords)
                }
            }
            
        except Exception as e:
            self.log(f"Error running UMAP: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _find_markers(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Find marker genes for each cluster."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            groupby = params.get('groupby', 'clusters')
            method = params.get('method', 'wilcoxon')
            n_genes = params.get('n_genes', 25)
            
            if groupby not in self.adata.obs.columns:
                return {"status": "error", "message": f"Grouping column '{groupby}' not found"}
            
            sc.tl.rank_genes_groups(self.adata, groupby=groupby, method=method, n_genes=n_genes)
            
            # Extract results
            markers = {}
            for cluster in self.adata.obs[groupby].unique():
                try:
                    genes = sc.get.rank_genes_groups_df(self.adata, group=str(cluster))
                    markers[str(cluster)] = {
                        "genes": genes['names'].tolist(),
                        "scores": genes['scores'].tolist(),
                        "pvals": genes['pvals'].tolist(),
                        "logfoldchanges": genes['logfoldchanges'].tolist()
                    }
                except:
                    markers[str(cluster)] = {"genes": [], "scores": [], "pvals": [], "logfoldchanges": []}
            
            self._save_checkpoint("after_markers")
            
            self.log(f"Marker gene analysis complete for {len(markers)} groups")
            return {
                "status": "success",
                "data": {
                    "groupby": groupby,
                    "method": method,
                    "markers": markers
                }
            }
            
        except Exception as e:
            self.log(f"Error finding markers: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _annotate_cells(self, method: str = "celltypist") -> Dict[str, Any]:
        """Annotate cell types."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        try:
            if method == "celltypist":
                try:
                    import celltypist
                    from celltypist import models
                    
                    # Download model if needed
                    models.download_models(force_update=False, model='Immune_All_Low.pkl')
                    model = models.Model.load(model='Immune_All_Low.pkl')
                    
                    predictions = celltypist.annotate(self.adata, model=model, majority_voting=True)
                    self.adata = predictions.to_adata()
                    
                    cell_types = self.adata.obs['majority_voting'].value_counts().to_dict()
                    
                    self.log(f"CellTypist annotation complete: {len(cell_types)} cell types")
                    return {
                        "status": "success",
                        "data": {
                            "method": "celltypist",
                            "model": "Immune_All_Low",
                            "cell_types": cell_types,
                            "annotations": self.adata.obs['majority_voting'].tolist()
                        }
                    }
                except ImportError:
                    self.log("CellTypist not available, using cluster labels")
                    return {
                        "status": "success",
                        "data": {
                            "method": "clusters",
                            "message": "CellTypist not installed, using cluster IDs as cell type labels",
                            "cell_types": self.adata.obs['clusters'].value_counts().to_dict() if 'clusters' in self.adata.obs else {}
                        }
                    }
            else:
                return {"status": "error", "message": f"Unknown annotation method: {method}"}
                
        except Exception as e:
            self.log(f"Error annotating cells: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _save_checkpoint(self, name: str) -> Dict[str, Any]:
        """Save current state as checkpoint."""
        if self.adata is None:
            return {"status": "error", "message": "No data to checkpoint."}
        
        try:
            checkpoint_path = os.path.join(self.checkpoint_dir, f"{name}.h5ad")
            self.adata.write_h5ad(checkpoint_path)
            self.log(f"Checkpoint saved: {name}")
            return {"status": "success", "path": checkpoint_path}
        except Exception as e:
            self.log(f"Error saving checkpoint: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _load_checkpoint(self, name: str) -> Dict[str, Any]:
        """Load state from checkpoint."""
        try:
            checkpoint_path = os.path.join(self.checkpoint_dir, f"{name}.h5ad")
            if not os.path.exists(checkpoint_path):
                return {"status": "error", "message": f"Checkpoint not found: {name}"}
            
            self.adata = sc.read_h5ad(checkpoint_path)
            self.log(f"Checkpoint loaded: {name}")
            return {"status": "success", "n_cells": int(self.adata.n_obs), "n_genes": int(self.adata.n_vars)}
        except Exception as e:
            self.log(f"Error loading checkpoint: {str(e)}")
            return {"status": "error", "message": str(e)}

    def _get_data_summary(self) -> Dict[str, Any]:
        """Get current data state summary."""
        if self.adata is None:
            return {"status": "error", "message": "No data loaded."}
        
        return {
            "status": "success",
            "data": {
                "n_cells": int(self.adata.n_obs),
                "n_genes": int(self.adata.n_vars),
                "obs_columns": list(self.adata.obs.columns),
                "var_columns": list(self.adata.var.columns),
                "obsm_keys": list(self.adata.obsm.keys()) if self.adata.obsm else [],
                "uns_keys": list(self.adata.uns.keys()) if self.adata.uns else [],
                "layers": list(self.adata.layers.keys()) if self.adata.layers else []
            }
        }
