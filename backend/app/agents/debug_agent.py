from typing import Dict, Any, Optional
import os
import datetime
import traceback
from .base_agent import BaseAgent

class DebugAgent(BaseAgent):
    def __init__(self):
        # Avoid recursion by not calling super().__init__ if it uses DebugAgent, 
        # but BaseAgent is simple so it is fine.
        super().__init__(agent_id="DBG-001", name="Debug Agent")
        self.log_file = os.path.join(os.getcwd(), 'data', 'debug.log')
        os.makedirs(os.path.dirname(self.log_file), exist_ok=True)
        self.log("Debug Agent Initialized")

    def execute(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Debug agent can respond to status checks or specific debug commands.
        """
        action = payload.get("action")
        if action == "status":
            return {"status": "active", "log_file": self.log_file}
        elif action == "report_error":
            source = payload.get("source", "unknown")
            error = payload.get("error", "No error provided")
            self.report_error(source, error)
            return {"status": "recorded"}
        
        return {"status": "ignored", "reason": "unknown action"}

    def report_error(self, source: str, error: str):
        """
        Log an error from another agent.
        """
        timestamp = datetime.datetime.now().isoformat()
        log_entry = f"[{timestamp}] [ERROR] [{source}] {error}\n"
        
        try:
            with open(self.log_file, "a") as f:
                f.write(log_entry)
            # Also print to stdout for immediate visibility in terminal
            print(f"!!! DEBUG AGENT CAUGHT ERROR !!! {log_entry.strip()}")
        except Exception as e:
            print(f"CRITICAL: DebugAgent failed to log error: {e}")

    def monitor(self):
        """
        Active monitoring method.
        In a real deployment, this could be run as a periodic task (e.g. via Celery beat)
        to check system health.
        """
        timestamp = datetime.datetime.now().isoformat()
        heartbeat_msg = f"[{timestamp}] [INFO] [System] DebugAgent is active and monitoring."
        try:
           with open(self.log_file, "a") as f:
                f.write(heartbeat_msg + "\n")
        except Exception:
            pass # resilient to logging errors
