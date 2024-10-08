# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import time
from time import sleep
from functools import wraps
from typing import Callable, Any
import asyncio

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

def timeit(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to measure and log the execution time of a function."""
    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = time.time()
        if asyncio.iscoroutinefunction(func):
            async def async_wrapper():
                result = await func(*args, **kwargs)
                end_time = time.time()
                LOG.debug(f"{func.__name__} executed in {end_time - start_time:.2f} seconds")
                return result
            return async_wrapper()
        else:
            result = func(*args, **kwargs)
            end_time = time.time()
            LOG.debug(f"{func.__name__} executed in {end_time - start_time:.2f} seconds")
            return result
    return wrapper
