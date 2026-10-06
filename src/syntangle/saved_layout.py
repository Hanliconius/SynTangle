"""Decode and validate complete saved display states without changing evidence."""
from .layout import LayoutState


def decode_saved_layout(fixture, data):
    refs = {(r.species_id, r.chromosome_id): r for r in fixture.chromosome_refs}
    labels = {r.label: r for r in fixture.chromosome_refs}
    return LayoutState({sp: tuple(refs[sp, chrom] for chrom in order)
                        for sp, order in data['chromosome_order'].items()},
                       {labels[label]: sign for label, sign in data['chromosome_orientation'].items()})


def validate_saved_layout(fixture, state, basis):
    if set(state.chromosome_order) != set(fixture.species_ids):
        raise ValueError('Saved layout species do not match fixture')
    for species, order in state.chromosome_order.items():
        expected = {r for r in fixture.chromosome_refs if r.species_id == species}
        if len(order) != len(expected) or set(order) != expected:
            raise ValueError('Saved layout must contain every chromosome exactly once')
    if set(state.chromosome_orientation) != set(fixture.chromosome_refs):
        raise ValueError('Saved layout orientations do not match chromosomes')
    if any(sign not in (-1, 1) for sign in state.chromosome_orientation.values()):
        raise ValueError('Saved layout orientation must be +1 or -1')
    for group in basis.free_flip_groups:
        if len({state.chromosome_orientation[r] * basis.base_assignment[r] for r in group}) != 1:
            raise ValueError('Saved layout violates hard orientation equations')
