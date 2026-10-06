"""Solve-local cooperative deadline and feasible-component progress channel."""
from contextvars import ContextVar
from functools import wraps
from time import monotonic

_deadline = ContextVar('syntangle_deadline', default=None)
_coupled_size = ContextVar('syntangle_coupled_bound_size', default=0)
_cluster_size = ContextVar('syntangle_bound_cluster_size', default=0)
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
        coupled_token = _coupled_size.set(kwargs.get('coupled_bound_size', 0))
        cluster_token = _cluster_size.set(kwargs.get('bound_cluster_size', 0))
        token = _deadline.set(None if seconds is None else monotonic() + seconds)
        try:
            return solve(*args, **kwargs)
        finally:
            _deadline.reset(token)
            _cluster_size.reset(cluster_token)
            _coupled_size.reset(coupled_token)
    return wrapped
