import sys
import os

# Add app to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app.agents.base_agent import BaseAgent
from app.agents.debug_agent import DebugAgent

class FaultyAgent(BaseAgent):
    def __init__(self):
        super().__init__("FAULTY-001", "Faulty Agent")
    
    def execute(self, payload):
        raise ValueError("Intentional Failure for Testing")

def test_debug_agent():
    print("Testing Debug Agent...")
    
    # 1. Create faulty agent
    agent = FaultyAgent()
    
    # 2. Execute safely
    result = agent.execute_safe({})
    
    print(f"Result: {result}")
    
    # 3. Check if error was reported
    log_file = os.path.join(os.getcwd(), 'data', 'debug.log')
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = f.read()
            if "Intentional Failure for Testing" in content:
                print("SUCCESS: Error found in debug log.")
            else:
                print("FAILURE: Error NOT found in debug log.")
                print(f"Log content: {content}")
    else:
        print(f"FAILURE: Log file {log_file} does not exist.")

if __name__ == "__main__":
    test_debug_agent()
