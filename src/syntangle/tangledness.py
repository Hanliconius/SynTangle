"""Presentation vs intrinsic crossing metrics under one fixed legal model."""
from operator import index


def _count(name, value):
    if isinstance(value, bool): raise ValueError(name+' must be an integer crossing count')
    try: value=index(value)
    except TypeError as exc: raise ValueError(name+' must be an integer crossing count') from exc
    if value<0: raise ValueError(name+' must be nonnegative')
    return value


def tangledness_metrics(display_crossings, best_feasible_crossings, lower_bound=0):
    """Report exact values only when a valid global bound closes the gap.

    The caller supplies scores of legal layouts on identical evidence and a
    global (not neighborhood-conditional) lower bound. The display is itself a
    feasible upper bound; it tightens an inferior supplied candidate's bound.
    This arithmetic does not validate or manufacture a solver certificate.
    """
    display=_count('display_crossings',display_crossings)
    candidate=_count('best_feasible_crossings',best_feasible_crossings)
    lower=_count('lower_bound',lower_bound)
    upper=min(display,candidate)
    if lower>upper:raise ValueError('Global lower bound exceeds a known feasible layout')
    exact=lower==upper
    return dict(
        visual_tangledness=display,
        intrinsic_tangledness=upper if exact else None,
        intrinsic_tangledness_bounds=dict(lower=lower,upper=upper),
        excess_layout_tangledness=display-upper if exact else None,
        excess_layout_tangledness_bounds=dict(lower=display-upper,upper=display-lower),
        demonstrated_avoidable_crossings=display-upper,
        best_retained_visual_tangledness=upper,
        best_retained_excess_bounds=dict(lower=0,upper=upper-lower),
        optimality_gap=upper-lower,
        optimality_established=exact,
        objective='unweighted adjacent-species homology-link crossings',
        scope='identical evidence, species layer order and hard constraints; legal whole-chromosome moves/flips',
        interpretation='Intrinsic under the layout model, not an evolutionary rearrangement count or a reader-perception measurement')
