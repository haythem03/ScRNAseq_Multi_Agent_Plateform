"""
Visualization Agent - Generates publication-quality plots and reports.
"""
from typing import Dict, Any, List, Optional
import os
import base64
from io import BytesIO
from .base_agent import BaseAgent

class VisualizationAgent(BaseAgent):
    def __init__(self):
        super().__init__(agent_id="VIZ-001", name="Visualization Agent")
        self.output_dir = os.path.join(os.getcwd(), 'data', 'visualizations')
        os.makedirs(self.output_dir, exist_ok=True)
    
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        action = payload.get("action")
        
        if action == "generate_qc_plots":
            return self._generate_qc_plots(payload.get("qc_data"))
        elif action == "generate_umap_plot":
            return self._generate_umap_plot(payload.get("embedding"), payload.get("clusters"))
        elif action == "generate_marker_plots":
            return self._generate_marker_plots(payload.get("markers"), payload.get("adata_path"))
        elif action == "generate_violin_plot":
            return self._generate_violin_plot(payload.get("data"), payload.get("gene"))
        else:
            return {"status": "error", "message": f"Unknown action: {action}"}
    
    def _generate_qc_plots(self, qc_data: Dict[str, Any]) -> Dict[str, Any]:
        """Generate QC visualization plots (violin, scatter) using pure matplotlib."""
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
            import numpy as np
            
            if not qc_data:
                self.log("No QC data provided for plotting")
                return {"status": "error", "message": "No QC data provided", "plots": {}}
            
            plots = {}
            
            # Create figure with 3 subplots for QC metrics
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            fig.patch.set_facecolor('#1e293b')
            
            for ax in axes:
                ax.set_facecolor('#0f172a')
                ax.tick_params(colors='#f8fafc')
                ax.spines['bottom'].set_color('#f8fafc')
                ax.spines['left'].set_color('#f8fafc')
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
            
            # Violin plot for n_genes_by_counts
            if 'n_genes_by_counts' in qc_data and qc_data['n_genes_by_counts']:
                data = [float(x) for x in qc_data['n_genes_by_counts'] if x is not None]
                if data:
                    violin1 = axes[0].violinplot(data, showmedians=True)
                    for pc in violin1['bodies']:
                        pc.set_facecolor('#6366f1')
                        pc.set_alpha(0.7)
                    violin1['cmedians'].set_color('#f8fafc')
                    violin1['cmins'].set_color('#f8fafc')
                    violin1['cmaxes'].set_color('#f8fafc')
                    violin1['cbars'].set_color('#f8fafc')
                axes[0].set_title('Genes per Cell', color='#f8fafc', fontsize=12)
                axes[0].set_ylabel('Count', color='#f8fafc')
                axes[0].set_xticks([])
            
            # Violin plot for total_counts
            if 'total_counts' in qc_data and qc_data['total_counts']:
                data = [float(x) for x in qc_data['total_counts'] if x is not None]
                if data:
                    violin2 = axes[1].violinplot(data, showmedians=True)
                    for pc in violin2['bodies']:
                        pc.set_facecolor('#ec4899')
                        pc.set_alpha(0.7)
                    violin2['cmedians'].set_color('#f8fafc')
                    violin2['cmins'].set_color('#f8fafc')
                    violin2['cmaxes'].set_color('#f8fafc')
                    violin2['cbars'].set_color('#f8fafc')
                axes[1].set_title('Total Counts per Cell', color='#f8fafc', fontsize=12)
                axes[1].set_ylabel('Count', color='#f8fafc')
                axes[1].set_xticks([])
            
            # Violin plot for pct_counts_mt
            if 'pct_counts_mt' in qc_data and qc_data['pct_counts_mt']:
                data = [float(x) for x in qc_data['pct_counts_mt'] if x is not None]
                if data:
                    violin3 = axes[2].violinplot(data, showmedians=True)
                    for pc in violin3['bodies']:
                        pc.set_facecolor('#10b981')
                        pc.set_alpha(0.7)
                    violin3['cmedians'].set_color('#f8fafc')
                    violin3['cmins'].set_color('#f8fafc')
                    violin3['cmaxes'].set_color('#f8fafc')
                    violin3['cbars'].set_color('#f8fafc')
                axes[2].set_title('% Mitochondrial', color='#f8fafc', fontsize=12)
                axes[2].set_ylabel('Percentage', color='#f8fafc')
                axes[2].set_xticks([])
            
            plt.tight_layout()
            
            # Save and encode
            qc_path = os.path.join(self.output_dir, 'qc_violin.png')
            fig.savefig(qc_path, dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            
            # Convert to base64 for frontend
            with open(qc_path, 'rb') as f:
                plots['qc_violin'] = base64.b64encode(f.read()).decode('utf-8')
            
            # Generate scatter plot (n_genes vs total_counts)
            if ('n_genes_by_counts' in qc_data and qc_data['n_genes_by_counts'] and
                'total_counts' in qc_data and qc_data['total_counts']):
                
                fig2, ax2 = plt.subplots(figsize=(8, 6))
                fig2.patch.set_facecolor('#1e293b')
                ax2.set_facecolor('#0f172a')
                ax2.spines['bottom'].set_color('#f8fafc')
                ax2.spines['left'].set_color('#f8fafc')
                ax2.spines['top'].set_visible(False)
                ax2.spines['right'].set_visible(False)
                
                x_data = [float(x) for x in qc_data['total_counts'] if x is not None]
                y_data = [float(x) for x in qc_data['n_genes_by_counts'] if x is not None]
                
                # Color by MT% if available
                if 'pct_counts_mt' in qc_data and qc_data['pct_counts_mt']:
                    c_data = [float(x) for x in qc_data['pct_counts_mt'] if x is not None]
                    scatter = ax2.scatter(x_data, y_data, c=c_data, cmap='viridis', alpha=0.6, s=5)
                    cbar = plt.colorbar(scatter, ax=ax2)
                    cbar.set_label('% Mitochondrial', color='#f8fafc')
                    cbar.ax.tick_params(colors='#f8fafc')
                else:
                    ax2.scatter(x_data, y_data, c='#6366f1', alpha=0.6, s=5)
                
                ax2.set_xlabel('Total Counts', color='#f8fafc')
                ax2.set_ylabel('Number of Genes', color='#f8fafc')
                ax2.set_title('Genes vs Counts (colored by %MT)', color='#f8fafc')
                ax2.tick_params(colors='#f8fafc')
                
                scatter_path = os.path.join(self.output_dir, 'qc_scatter.png')
                fig2.savefig(scatter_path, dpi=150, facecolor='#1e293b', bbox_inches='tight')
                plt.close(fig2)
                
                with open(scatter_path, 'rb') as f:
                    plots['qc_scatter'] = base64.b64encode(f.read()).decode('utf-8')
            
            self.log("QC plots generated successfully")
            return {"status": "success", "plots": plots, "paths": [qc_path]}
            
        except Exception as e:
            self.log(f"Error generating QC plots: {str(e)}")
            import traceback
            traceback.print_exc()
            return {"status": "error", "message": str(e), "plots": {}}
    
    def _generate_umap_plot(self, embedding: List[List[float]], clusters: List[int]) -> Dict[str, Any]:
        """Generate UMAP visualization with cluster coloring."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import numpy as np
            
            fig, ax = plt.subplots(figsize=(10, 8))
            fig.patch.set_facecolor('#1e293b')
            ax.set_facecolor('#0f172a')
            
            embedding = np.array(embedding)
            clusters = np.array(clusters)
            
            # Use a colorful palette
            n_clusters = len(np.unique(clusters))
            colors = plt.cm.tab20(np.linspace(0, 1, max(n_clusters, 20)))
            
            for i in np.unique(clusters):
                mask = clusters == i
                ax.scatter(
                    embedding[mask, 0], 
                    embedding[mask, 1],
                    c=[colors[i % 20]],
                    label=f'Cluster {i}',
                    s=5,
                    alpha=0.7
                )
            
            ax.set_xlabel('UMAP 1', color='#f8fafc')
            ax.set_ylabel('UMAP 2', color='#f8fafc')
            ax.set_title('UMAP Clustering', color='#f8fafc', fontsize=14)
            ax.tick_params(colors='#f8fafc')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                     facecolor='#1e293b', labelcolor='#f8fafc')
            
            umap_path = os.path.join(self.output_dir, 'umap_clusters.png')
            fig.savefig(umap_path, dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            
            with open(umap_path, 'rb') as f:
                plot_b64 = base64.b64encode(f.read()).decode('utf-8')
            
            self.log("UMAP plot generated successfully")
            return {"status": "success", "plot": plot_b64, "path": umap_path}
            
        except Exception as e:
            self.log(f"Error generating UMAP plot: {str(e)}")
            return {"status": "error", "message": str(e)}
    
    def _generate_marker_plots(self, markers: Dict[str, Any], adata_path: str) -> Dict[str, Any]:
        """Generate marker gene dot plots and heatmaps."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import scanpy as sc
            
            # Load adata for expression data
            adata = sc.read_h5ad(adata_path) if adata_path else None
            
            # For now, create a simple bar plot of marker significance
            if markers and 'top_genes' in markers:
                fig, ax = plt.subplots(figsize=(12, 6))
                fig.patch.set_facecolor('#1e293b')
                ax.set_facecolor('#0f172a')
                
                genes = markers['top_genes'][:20]  # Top 20 markers
                scores = markers.get('scores', list(range(len(genes))))[:20]
                
                bars = ax.barh(genes, scores, color='#6366f1')
                ax.set_xlabel('Score', color='#f8fafc')
                ax.set_title('Top Marker Genes', color='#f8fafc')
                ax.tick_params(colors='#f8fafc')
                
                marker_path = os.path.join(self.output_dir, 'markers_bar.png')
                fig.savefig(marker_path, dpi=150, facecolor='#1e293b', bbox_inches='tight')
                plt.close(fig)
                
                with open(marker_path, 'rb') as f:
                    plot_b64 = base64.b64encode(f.read()).decode('utf-8')
                
                return {"status": "success", "plot": plot_b64, "path": marker_path}
            
            return {"status": "success", "plots": {}, "message": "No markers to plot"}
            
        except Exception as e:
            self.log(f"Error generating marker plots: {str(e)}")
            return {"status": "error", "message": str(e)}
    
    def _generate_violin_plot(self, data: List[float], gene: str) -> Dict[str, Any]:
        """Generate a single violin plot for a gene."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots(figsize=(6, 5))
            fig.patch.set_facecolor('#1e293b')
            ax.set_facecolor('#0f172a')
            
            violin = ax.violinplot(data, showmedians=True)
            for pc in violin['bodies']:
                pc.set_facecolor('#6366f1')
                pc.set_alpha(0.7)
            
            ax.set_title(f'{gene} Expression', color='#f8fafc')
            ax.set_ylabel('Expression', color='#f8fafc')
            ax.tick_params(colors='#f8fafc')
            
            path = os.path.join(self.output_dir, f'violin_{gene}.png')
            fig.savefig(path, dpi=150, facecolor='#1e293b', bbox_inches='tight')
            plt.close(fig)
            
            with open(path, 'rb') as f:
                plot_b64 = base64.b64encode(f.read()).decode('utf-8')
            
            return {"status": "success", "plot": plot_b64, "path": path}
            
        except Exception as e:
            return {"status": "error", "message": str(e)}
