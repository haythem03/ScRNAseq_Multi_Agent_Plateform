from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import traceback

class BaseAgent(ABC):
    def __init__(self, agent_id: str, name: str):
        self.agent_id = agent_id
        self.name = name
    
    def execute_safe(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Wrapper to ensure errors are caught and reported.
        """
        try:
            return self.execute(payload)
        except Exception as e:
            error_msg = f"{str(e)}\n{traceback.format_exc()}"
            self.log(f"Error executing agent: {error_msg}")
            
            # Dynamically import to avoid circular dependency issues at module level
            from .debug_agent import DebugAgent
            # Prevent infinite recursion if DebugAgent itself fails
            if not isinstance(self, DebugAgent):
                try:
                    debugger = DebugAgent()
                    debugger.report_error(self.name, error_msg)
                except Exception as log_err:
                    print(f"CRITICAL: Failed to report error to DebugAgent: {log_err}")
            
            return {"status": "error", "message": str(e), "details": error_msg}

    @abstractmethod
    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute the agent's primary function.
        :param payload: input data/parameters
        :return: result dictionary
        """
        pass
    
    def log(self, message: str):
        print(f"[{self.name}] {message}")
