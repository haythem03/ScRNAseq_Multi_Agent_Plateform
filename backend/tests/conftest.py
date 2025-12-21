"""
Pytest configuration for backend tests.
Ensures the test scanpy stub is used instead of the real scanpy package.
"""
import sys
import os

# Add the scanpy_stub to sys.path BEFORE the real scanpy
stub_path = os.path.join(os.path.dirname(__file__), 'scanpy_stub')
if stub_path not in sys.path:
    sys.path.insert(0, os.path.dirname(stub_path))

# Also ensure the parent directory (backend) is in path
backend_path = os.path.dirname(os.path.dirname(__file__))
if backend_path not in sys.path:
    sys.path.insert(0, backend_path)
