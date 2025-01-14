from concurrent.futures import ThreadPoolExecutor
import time

# Example function to be executed by threads
def worker_function(task_id):
    print(f"Task {task_id} is starting.")
    time.sleep(2)  # Simulate a long-running task
    print(f"Task {task_id} is complete.")

# Create a ThreadPoolExecutor with 4 threads
with ThreadPoolExecutor(max_workers=4) as executor:
    # Submit tasks to the thread pool
    tasks = [executor.submit(worker_function, i) for i in range(10)]

print("All tasks have been submitted.")
