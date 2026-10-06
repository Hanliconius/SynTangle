"""Admissible graph-factor relaxation coupling relative orders across layers.

A legal chromosome layout assigns every pair-order bit consistently. Relaxing
transitivity enlarges that space; joint minimization within disjoint factor
buckets remains a lower bound. Neighbor orders are never fixed to an incumbent.
"""
from itertools import combinations, product
from .crossing_cache import _decision
from .layout import _occurrences_by_species_homology, AmbiguousHomologyError
from .incidence import chromosome_node_id


def bucket_tables(factors, size):
    """Partition variables; assign every nonnegative factor exactly once."""
    remaining = {v for scope, _ in factors for v in scope}
    adjacency = {v:{} for v in remaining}
    for scope,values in factors:
        weight=max(values.values())-min(values.values())
        for a,b in combinations(scope,2):
            adjacency[a][b]=adjacency[a].get(b,0)+weight
            adjacency[b][a]=adjacency[b].get(a,0)+weight
    clusters=[]
    owner={}
    while remaining:
        first=min(remaining)
        cluster={first}; affinity={}
        def add_neighbors(v):
            for other,weight in adjacency[v].items():
                if other in remaining and other not in cluster:
                    affinity[other]=affinity.get(other,0)+weight
        add_neighbors(first)
        while len(cluster)<size and affinity:
            candidate=max(affinity,key=lambda v:(affinity[v],-v))
            if not affinity[candidate]:
                break
            del affinity[candidate]
            cluster.add(candidate);add_neighbors(candidate)
        remaining-=cluster
        scope=tuple(sorted(cluster))
        for v in scope:
            owner[v]=len(clusters)
        clusters.append(scope)
    inside_by_cluster={};pending=[]
    for scope,table in factors:
        owners={owner[v] for v in scope}
        if len(owners)==1:
            inside_by_cluster.setdefault(next(iter(owners)),[]).append((scope,table))
        else:
            pending.append((scope,table))
    output=[]
    for i,scope in enumerate(clusters):
        inside=inside_by_cluster.get(i,[])
        if not inside:
            continue
        position = {v:i for i,v in enumerate(scope)}
        table = {bits:sum(t[tuple(bits[position[v]] for v in s)] for s,t in inside)
                 for bits in product((0,1), repeat=len(scope))}
        def minimum(bits):
            if bits not in table:
                i=bits.index(None)
                table[bits]=min(minimum(bits[:i]+(b,)+bits[i+1:]) for b in (0,1))
            return table[bits]
        for bits in product((None,0,1),repeat=len(scope)):
            minimum(bits)
        output.append((scope,table))
    # Factors crossing clusters retain their own joint two-bit minimum, not
    # duplicated into either cluster. Partial entries include correlated bits.
    for scope, full in pending:
        table={bits:min(value for complete,value in full.items()
                       if all(bit is None or bit==v for bit,v in zip(bits,complete)))
               for bits in product((None,0,1),repeat=len(scope))}
        output.append((scope,table))
    return tuple(output)


class CoupledCrossingBound:
    def __init__(self, fixture, nodes, independent, size=6):
        if not 1 <= size <= 8:
            raise ValueError('Coupled bound size must be 1..8')
        self.independent = independent
        self.basis = independent.basis
        self.order_variables = {}
        next_variable = len(self.basis.free_flip_groups)
        index = _occurrences_by_species_homology(fixture)
        aggregate = {}
        for left,right in zip(fixture.species_ids,fixture.species_ids[1:]):
            links=[]
            for homology in sorted(set(index.get(left,{})) & set(index.get(right,{}))):
                a,b=index[left][homology],index[right][homology]
                if len(a)!=1 or len(b)!=1:
                    raise AmbiguousHomologyError('Coupled bound requires one-to-one homology')
                if all(chromosome_node_id(occ[0].ref) in nodes for occ in (a[0],b[0])):
                    links.append((a[0],b[0]))
            for a,b in combinations(links,2):
                decisions=[]
                for endpoint in (0,1):
                    decision,relations=_decision(a[endpoint],b[endpoint])
                    if len(decision)==1:
                        ref=decision[0]
                        variable=independent.group_index_by_ref[ref]
                        if self.basis.base_assignment[ref]==-1:
                            relations=tuple(reversed(relations))
                    else:
                        if decision not in self.order_variables:
                            self.order_variables[decision]=next_variable
                            next_variable+=1
                        variable=self.order_variables[decision]
                    decisions.append((variable,relations))
                (va,ra),(vb,rb)=decisions
                scope=tuple(sorted({va,vb}))
                table=aggregate.setdefault(scope,{bits:0 for bits in product((0,1),repeat=len(scope))})
                for bits in table:
                    assignment=dict(zip(scope,bits))
                    table[bits]+=int(ra[assignment[va]]*rb[assignment[vb]]<0)
        self.variable_count=next_variable
        self.factor_count=len(aggregate)
        self.constant=sum(next(iter(values.values())) for values in aggregate.values()
                          if len(set(values.values()))==1)
        self.tables=bucket_tables([(scope,values) for scope,values in aggregate.items()
                                   if len(set(values.values()))>1],size)

    def lower_bound(self,bits,orders=None):
        if len(bits)!=len(self.basis.free_flip_groups):
            raise ValueError('Orientation bit vector does not match basis')
        assignment=list(bits)+[None]*(self.variable_count-len(bits))
        ranks={ref:i for order in (orders or {}).values() for i,ref in enumerate(order)}
        for (a,b),variable in (self.order_variables.items() if ranks else ()):
            if a in ranks and b in ranks:
                assignment[variable]=int(ranks[a]>ranks[b])
        coupled=self.constant+sum(table[tuple(assignment[v] for v in scope)] for scope,table in self.tables)
        # Both are bounds on the SAME complete objective. Never add them.
        return max(self.independent.lower_bound(bits),coupled)

    def for_orders(self,orientation,orders):
        bits=tuple(int(orientation[group[0]]!=self.basis.base_assignment[group[0]])
                   for group in self.basis.free_flip_groups)
        return self.lower_bound(bits,orders)
