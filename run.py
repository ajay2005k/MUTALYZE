import subprocess
import sys
import os
import time
from threading import Thread

def run_backend():
    os.chdir('api')
    subprocess.run([sys.executable, '-m', 'uvicorn', 'main:app', '--reload', '--host', '0.0.0.0', '--port', '8000'])

def run_frontend():
    os.chdir('..')
    subprocess.run([sys.executable, '-m', 'streamlit', 'run', 'app.py'])

if __name__ == "__main__":
    # Create necessary directories
    os.makedirs('models', exist_ok=True)
    
    # Start backend in a separate thread
    backend_thread = Thread(target=run_backend)
    backend_thread.daemon = True
    backend_thread.start()
    
    # Wait for backend to start
    print("Starting backend server...")
    time.sleep(5)
    
    # Start frontend
    print("Starting frontend server...")
    run_frontend() 