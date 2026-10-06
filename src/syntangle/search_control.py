"""Solve-local cooperative deadline and feasible-component progress channel."""
from contextvars import ContextVar
from functools import wraps
from time import monotonic

_deadline = ContextVar('syntangle_deadline', default=None)
_progress = ContextVar('syntangle_component_progress', default=None)


def expired():
    deadline = _deadline.get()
    return deadline is not None and monotonic() >= deadline


def emit(state):
    callback = _progress.get()
    if callback is not None:
        callback(state)


def controlled(solve):
    @wraps(solve)
    def wrapped(*args, **kwargs):
        seconds = kwargs.get('time_limit_seconds')
        if seconds is not None and seconds <= 0:
            raise ValueError('time_limit_seconds must be positive')
        token = _deadline.set(None if seconds is None else monotonic() + seconds)
        try:
            return solve(*args, **kwargs)
        finally:
            _deadline.reset(token)
    return wrapped
