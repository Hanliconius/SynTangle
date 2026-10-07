"""Decision components include hard couplings, without changing homology evidence."""
from .incidence import build_incidence_graph, chromosome_node_id


def optimization_component_map(fixture):
    """Merge homology components connected by hard orientation equations.

    Independence requires both objective and constraint separation. Raw incidence
    components/cycle metrics remain biological quantities and are not modified.
    The merged node sets retain every original node for scoring/reconstruction.
    """
    raw=build_incidence_graph(fixture).connected_components()
    by_node={n:i for i,nodes in enumerate(raw) for n in nodes}
    parent=list(range(len(raw)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for c in fixture.orientation_constraints:
        a,b=root(by_node[chromosome_node_id(c.a)]),root(by_node[chromosome_node_id(c.b)])
        if a!=b:parent[max(a,b)]=min(a,b)
    merged={}
    for i,nodes in enumerate(raw):merged.setdefault(root(i),set()).update(nodes)
    components=tuple(sorted((frozenset(nodes) for nodes in merged.values()),key=min))
    index={n:i for i,nodes in enumerate(components) for n in nodes}
    return components,{ref:index[chromosome_node_id(ref)] for ref in fixture.chromosome_refs}
