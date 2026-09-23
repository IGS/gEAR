# resource_limits.py - Process memory-limit management.

"""
resource_limits.py - Process memory-limit management.
"""

import functools
import os
import sys
import typing


def set_memory_limit_from_cgroup(fraction: float = 0.9) -> None:
    """
    Self-impose a process memory ceiling based on the container's cgroup memory limit,
    so that approaching the real limit raises a catchable MemoryError instead of the
    kernel OOM-killer sending an uncatchable SIGKILL (exit code 137).

    Reads the cgroup v2 limit first (/sys/fs/cgroup/memory.max), falling back to
    cgroup v1 (/sys/fs/cgroup/memory/memory.limit_in_bytes). If no bounded limit is
    found (unbounded, or the files aren't present), this is a no-op and processes
    remain subject to the OS OOM-killer as before.

    This is best-effort defensive setup: any failure here is caught and logged
    rather than raised, so it can never prevent a caller from starting up.

    Parameters
    ----------
    fraction : float, optional (default: 0.9)
        Fraction of the detected cgroup memory limit to set as the RLIMIT_AS ceiling.
    """
    import resource

    # Cgroup v1 reports an implausibly large number (close to 2**63, platform max)
    # to mean "no limit" rather than a sentinel string like v2's "max".
    UNBOUNDED_V1_THRESHOLD = 2**62

    limit_bytes = None

    try:
        cgroup_v2_path = "/sys/fs/cgroup/memory.max"
        cgroup_v1_path = "/sys/fs/cgroup/memory/memory.limit_in_bytes"

        if os.path.exists(cgroup_v2_path):
            with open(cgroup_v2_path) as f:
                value = f.read().strip()
            if value != "max":
                limit_bytes = int(value)
        elif os.path.exists(cgroup_v1_path):
            with open(cgroup_v1_path) as f:
                value = int(f.read().strip())
            if value < UNBOUNDED_V1_THRESHOLD:
                limit_bytes = value

        if limit_bytes is None:
            print(
                "set_memory_limit_from_cgroup: no bounded cgroup memory limit found; "
                "not setting a self-imposed RLIMIT_AS ceiling.",
                file=sys.stderr,
            )
            return

        target = int(limit_bytes * fraction)
        resource.setrlimit(resource.RLIMIT_AS, (target, target))
        print(
            f"set_memory_limit_from_cgroup: cgroup limit is {limit_bytes} bytes; "
            f"set RLIMIT_AS to {target} bytes ({fraction:.0%}).",
            file=sys.stderr,
        )
    except Exception as e:
        print(f"set_memory_limit_from_cgroup: failed to set memory limit: {e}", file=sys.stderr)


def catch_memory_error() -> typing.Callable:
    """
    A decorator factory that catches MemoryError exceptions in the decorated function.

    Returns:
        Callable: A decorator that wraps the target function. If a MemoryError is raised during
        execution, it prints an error message to stderr and returns a tuple containing a result
        dictionary and a 500 status code.

    Example:
        @catch_memory_error()
        def my_function():
            # function implementation
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except MemoryError as e:
                print(f"Exceeded memory in {func.__name__}: {e}", file=sys.stderr)

                result = {
                    "message": "Exceeded memory limit",
                    "success": -1,
                    "error": str(e),
                }

                return result, 500

        return wrapper

    return decorator
