"""Lightweight stub of the `scanpy` package for tests.

This module provides simple placeholders that are intended to be patched
by the test suite (using `unittest.mock.patch`). It avoids requiring the
full `scanpy` installation during unit tests.
"""
def read_h5ad(*args, **kwargs):
    raise ImportError("scanpy.read_h5ad is not implemented in stub; tests should patch this")

def read(*args, **kwargs):
    raise ImportError("scanpy.read is not implemented in stub; tests should patch this")


class _PP:
    @staticmethod
    def calculate_qc_metrics(adata, qc_vars=None, percent_top=None, log1p=False, inplace=True):
        """Minimal implementation to populate expected QC fields on `adata.obs`.

        This implementation is intentionally lightweight: it sets a few
        columns that `ExecutionAgent._load_and_qc` expects. Tests typically
        patch `scanpy.read_h5ad` to return a MagicMock AnnData, so this
        function attempts to work with either a real AnnData or a mock.
        """
        try:
            n = int(getattr(adata, 'n_obs', 0) or 0)
        except Exception:
            n = 0

        # Prepare simple default columns
        try:
            # If obs acts like a dict/DataFrame, attempt to assign lists
            adata.obs['n_genes_by_counts'] = [1] * n
            adata.obs['total_counts'] = [100] * n
            adata.obs['pct_counts_mt'] = [0.0] * n
        except Exception:
            # Best-effort: set attributes that tests may read
            try:
                setattr(adata, 'obs', {'n_genes_by_counts': [1] * n,
                                       'total_counts': [100] * n,
                                       'pct_counts_mt': [0.0] * n})
            except Exception:
                pass


pp = _PP()
