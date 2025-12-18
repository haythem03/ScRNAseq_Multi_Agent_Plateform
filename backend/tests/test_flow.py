import pytest
from app.agents.program_manager import ProgramManager
from app.agents.execution_agent import ExecutionAgent
from unittest.mock import MagicMock, patch

def test_program_manager_orchestration():
    pm = ProgramManager()
    
    # Mock ExecutionAgent to avoid real file I/O and Scanpy overhead
    with patch('app.agents.program_manager.ExecutionAgent') as MockExecAgent:
        mock_instance = MockExecAgent.return_value
        mock_instance.execute.return_value = {"status": "success", "data": {"n_cells": 100, "n_genes": 200}}
        
        payload = {
            "task_type": "upload_and_qc",
            "file_path": "/dummy/path/data.h5ad"
        }
        
        result = pm.execute(payload)
        
        assert result["status"] == "completed"
        assert result["step"] == "qc"
        assert result["qc_results"]["data"]["n_cells"] == 100
        
        # Verify call arguments
        mock_instance.execute.assert_called_with({
            "action": "load_and_qc",
            "file_path": "/dummy/path/data.h5ad"
        })

def test_execution_agent_logic():
    # Test valid action
    agent = ExecutionAgent()
    # We mock os.path.exists and sc.read
    with patch('os.path.exists') as mock_exists, \
         patch('scanpy.read_h5ad') as mock_read, \
         patch('scanpy.read') as mock_read_generic:
        
        mock_exists.return_value = True
        
        # Mock AnnData
        mock_adata = MagicMock()
        mock_adata.n_obs = 50
        mock_adata.n_vars = 10
        mock_adata.obs = MagicMock()
        
        mock_read.return_value = mock_adata
        
        result = agent.execute({"action": "load_and_qc", "file_path": "test.h5ad"})
        
        assert result["status"] == "success"
        assert result["data"]["n_cells"] == 50
