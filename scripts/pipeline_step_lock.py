"""Serialize concurrent invocations of a pipeline step within one project."""
from contextlib import contextmanager
import fcntl
from config import INTERIM
@contextmanager
def step_lock(number):
    path=INTERIM/f'.step_{number}.lock'
    with path.open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print(f'Step {number}: waiting for the concurrent invocation in this project',flush=True)
            fcntl.flock(handle,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(handle,fcntl.LOCK_UN)
