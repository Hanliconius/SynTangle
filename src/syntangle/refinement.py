"""Validated, never-worse refinement of a supplied legal reference layout.

GENESPACE is one source of a reference layout; this module does not invoke it
or claim that synthetic reference layouts are native GENESPACE results.
"""
from .audit import fixture_fingerprint
from .layout import score_crossings
from .orientation_space import orientation_basis
from .saved_layout import decode_saved_layout,validate_saved_layout
from .tangledness import tangledness_metrics


def select_refinement(fixture,reference,candidate=None,*,lower_bound=0,reference_name='GENESPACE'):
    """Reject illegal states; retain the reference on ties and regressions.

    Both layouts are scored canonically against the same fixture. A supplied
    bound must be global for that fixture; conditional bounds are inadmissible.
    Returning the reference on failure is not evidence of a solver improvement.
    """
    basis=orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,reference,basis)
    before=score_crossings(fixture,reference).crossings
    raw=None
    if candidate is not None:
        validate_saved_layout(fixture,candidate,basis)
        raw=score_crossings(fixture,candidate).crossings
    improved=raw is not None and raw<before
    selected=candidate if improved else reference
    after=raw if improved else before
    outcome=('improved' if improved else 'no candidate; reference retained' if raw is None
             else 'tie; reference retained' if raw==before else 'regression; reference retained')
    return dict(reference_name=reference_name,input_fingerprint=fixture_fingerprint(fixture),
        reference_crossings=before,raw_candidate_crossings=raw,returned_crossings=after,
        improved=improved,outcome=outcome,reference_retained=not improved,
        selected_state=selected.to_dict(),
        reference_metrics=tangledness_metrics(before,after,lower_bound),
        returned_metrics=tangledness_metrics(after,after,lower_bound),
        comparison_scope='Same crossing objective and legal evidence; reference-species row is not fixed by this adapter')


def refine_reference_layout(fixture,reference,seconds=30,*,method='hybrid_3',seed=1,
                            reference_name='GENESPACE',**options):
    """One fresh hybrid run from a supplied reference, followed by checked selection.

    A layout is an incumbent, not a previous branch frontier. This adapter never
    claims to import exclusions from an earlier optimization run. Defaults do
    not enable the currently unbenchmarked mirror/policy experiment.
    """
    from .hybrid_experiment import run_hybrid
    basis=orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,reference,basis)
    fingerprint=fixture_fingerprint(fixture)
    saved_reference=reference.to_dict()
    # A separate state prevents an optimizer from mutating the caller's baseline.
    start=decode_saved_layout(fixture,saved_reference)
    raw=run_hybrid(fixture,start,method,seconds,seed,**options)
    if fixture_fingerprint(fixture)!=fingerprint:
        raise AssertionError('Refinement changed input evidence')
    retained=decode_saved_layout(fixture,saved_reference)
    candidate=decode_saved_layout(fixture,raw['optimized_state'])
    validate_saved_layout(fixture,candidate,basis)
    score=score_crossings(fixture,candidate).crossings
    if score!=raw['upper_bound']:raise AssertionError('Solver upper bound disagrees with canonical candidate score')
    decision=select_refinement(fixture,retained,candidate,lower_bound=raw['lower_bound'],reference_name=reference_name)
    decision['solver_result']=raw
    return decision
