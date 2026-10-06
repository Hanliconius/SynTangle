"""Experimental joint order/flip models, dispatched before combinatorial search.

Shares incidence/GF(2) reductions. Not a restart/resume API for reduced search
frontiers: callers may supply a feasible incumbent, never a prior search ledger.
SDP objectives are numerical diagnostics, not certified pruning bounds.
"""
from itertools import combinations
import warnings
import math
import time
import numpy as np
from .heuristic import _component_map
from .orientation_space import orientation_basis
from .bounds import build_relaxed_crossing_bound
from .coupled_bound import CoupledCrossingBound
from .layout import LayoutState, canonicalize_component_order, score_crossings
from .saved_layout import validate_saved_layout


class JointModel:
    def __init__(self, fixture, nodes, refs):
        self.fixture, self.nodes = fixture, nodes
        self.basis = orientation_basis(fixture, tuple(sorted(refs)))
        independent = build_relaxed_crossing_bound(fixture, nodes, self.basis)
        self.graph = CoupledCrossingBound(fixture, nodes, independent, prepare_tables=False)
        self.orders = dict(self.graph.order_variables)
        self.rows = {sp: tuple(sorted(r for r in refs if r.species_id == sp))
                     for sp in fixture.species_ids}
        n = self.graph.variable_count
        # Zero-cost pair decisions are required for complete total-order consistency.
        for row in self.rows.values():
            for pair in combinations(row, 2):
                if pair not in self.orders:
                    self.orders[pair] = n
                    n += 1
        self.n = n
        self.constant = self.graph.constant
        self.linear = np.zeros(n)
        self.quadratic = {}
        for scope, table in self.graph.factors:
            if len(scope) == 1:
                a, = scope
                self.constant += table[0,]
                self.linear[a] += table[1,] - table[0,]
            else:
                a, b = scope
                t00,t01,t10,t11 = (table[x] for x in ((0,0),(0,1),(1,0),(1,1)))
                self.constant += t00
                self.linear[a] += t10-t00
                self.linear[b] += t01-t00
                q = t11-t10-t01+t00
                if q:
                    self.quadratic[a,b] = q
        self.triangles = tuple((self.orders[a,b],self.orders[b,c],self.orders[a,c])
                              for row in self.rows.values() for a,b,c in combinations(row,3))

    def encode(self, state):
        bits = np.zeros(self.n)
        for i, group in enumerate(self.basis.free_flip_groups):
            bits[i] = int(state.chromosome_orientation[group[0]] != self.basis.base_assignment[group[0]])
        ranks = {ref:i for row in state.chromosome_order.values() for i,ref in enumerate(row)}
        for (a,b),i in self.orders.items():
            bits[i] = int(ranks[a] > ranks[b])
        return bits

    def objective(self, bits):
        return self.constant + self.linear@bits + sum(q*bits[a]*bits[b] for (a,b),q in self.quadratic.items())

    def decode(self, bits, state):
        orders = dict(state.chromosome_order)
        for sp,row in self.rows.items():
            ranks = {r:0 for r in row}
            for (a,b),i in self.orders.items():
                if a.species_id == sp:
                    ranks[a if bits[i] > .5 else b] += 1
            if sorted(ranks.values()) != list(range(len(row))):
                raise AssertionError('Solver returned a cyclic/incomplete chromosome ordering')
            replacement = iter(sorted(row,key=ranks.get))
            orders[sp] = tuple(next(replacement) if r in ranks else r for r in orders[sp])
        signs = dict(state.chromosome_orientation)
        for i,group in enumerate(self.basis.free_flip_groups):
            for r in group:
                signs[r] = self.basis.base_assignment[r] * (-1 if bits[i] > .5 else 1)
        result = LayoutState(orders,signs)
        validate_saved_layout(self.fixture,result,orientation_basis(self.fixture,self.fixture.chromosome_refs))
        if abs(self.objective(np.rint(bits))-score_crossings(self.fixture,result,restrict_component_nodes=self.nodes).crossings) > 1e-6:
            raise AssertionError('Joint model differs from canonical crossing scorer')
        return result

    def linear_model(self):
        if hasattr(self,'_linear_cache'):
            return self._linear_cache
        matrix_started=time.perf_counter()
        from scipy.sparse import coo_matrix
        pairs = list(self.quadratic)
        c = np.r_[self.linear, [self.quadratic[p] for p in pairs]]
        rows,cols,values,lower,upper = [],[],[],[],[]
        def add(entries,lo,hi):
            row = len(lower)
            for col,value in entries:
                rows.append(row);cols.append(col);values.append(value)
            lower.append(lo);upper.append(hi)
        for a,b,d in self.triangles:
            add(((a,1),(b,1),(d,-1)),0,1)
        for z,(a,b) in enumerate(pairs,self.n):
            add(((z,1),(a,-1)),-np.inf,0)
            add(((z,1),(b,-1)),-np.inf,0)
            add(((z,1),(a,-1),(b,-1)),-1,np.inf)
        matrix = coo_matrix((values,(rows,cols)),shape=(len(lower),len(c))).tocsr()
        self.linear_model_build_seconds=time.perf_counter()-matrix_started
        self._linear_cache=(c,matrix,np.asarray(lower),np.asarray(upper),pairs)
        return self._linear_cache

    def milp(self, incumbent, seconds, fixed=None):
        build_started=time.perf_counter()
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import vstack, csr_matrix
        c,matrix,lower,upper,pairs=self.linear_model()
        warm=self.encode(incumbent)
        matrix=vstack([matrix,csr_matrix(c.reshape(1,-1))],format='csc')
        lower=np.r_[lower,-np.inf]
        upper=np.r_[upper,self.objective(warm)-self.constant]
        lb,ub=np.zeros(len(c)),np.ones(len(c))
        for i,value in (fixed or {}).items():lb[i]=ub[i]=value
        api_preparation_seconds=time.perf_counter()-build_started
        solver_started=time.perf_counter()
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',message='Unrecognized options detected.*')
            result=milp(c,integrality=np.r_[np.ones(self.n),np.zeros(len(pairs))],
                bounds=Bounds(lb,ub),constraints=LinearConstraint(matrix,lower,upper),
                options={'time_limit':max(.001,seconds-(time.perf_counter()-build_started)),
                         'mip_rel_gap':0.0,'threads':1})
        return self._finish_milp(incumbent,result.x,getattr(result,'mip_dual_bound',None),dict(
            status=int(result.status),message=result.message,
            api_preparation_seconds=api_preparation_seconds,solver_seconds=time.perf_counter()-solver_started,
            linear_model_build_seconds=self.linear_model_build_seconds,
            nodes=int(getattr(result,'mip_node_count',0) or 0),variables=len(c),constraints=len(lower),
            bound_scope='restricted neighborhood' if fixed else 'whole component',
            solver='scipy.optimize.milp/HiGHS'))

    def _finish_milp(self, incumbent, solution, dual, info):
        state=incumbent;warm=self.encode(incumbent)
        if solution is not None:
            bits=np.asarray(solution[:self.n])
            if np.max(abs(bits-np.rint(bits))) > 1e-5:
                raise AssertionError('Nonintegral MILP incumbent')
            candidate=self.decode(bits,incumbent)
            if self.objective(self.encode(candidate)) < self.objective(warm):state=candidate
        bound=max(0,math.ceil(float(dual)+self.constant-1e-5)) if dual is not None and np.isfinite(dual) else 0
        score=int(round(self.objective(self.encode(state))))
        if bound>score:raise AssertionError('MILP lower bound exceeds retained incumbent')
        info.update(raw_dual_without_constant=float(dual) if dual is not None and np.isfinite(dual) else None,
                    objective_constant=float(self.constant),component_lower=bound,component_upper=score)
        return state,bound,info

    def highs(self, incumbent, seconds, *, strict_proof=False, progress=None, use_mip_start=True):
        """One global solve with a real feasible MIP start; no restart/handoff."""
        build_started=time.perf_counter()
        import highspy
        from scipy.sparse import vstack, csr_matrix
        c,matrix,lower,upper,pairs=self.linear_model()
        bits=self.encode(incumbent);cutoff=self.objective(bits)-self.constant-int(strict_proof)
        matrix=vstack([matrix,csr_matrix(c.reshape(1,-1))],format='csr')
        lower=np.r_[lower,-np.inf];upper=np.r_[upper,cutoff]
        h=highspy.Highs()
        def checked(status):
            if status==highspy.HighsStatus.kError:raise RuntimeError('HiGHS API error')
        for key,value in [('output_flag',False),('threads',1),('mip_rel_gap',0.0)]:checked(h.setOptionValue(key,value))
        checked(h.addCols(len(c),c,np.zeros(len(c)),np.ones(len(c)),0,np.zeros(len(c)+1,dtype=np.int32),np.array([],dtype=np.int32),np.array([],dtype=float)))
        checked(h.addRows(len(lower),lower,upper,len(matrix.data),matrix.indptr.astype(np.int32),matrix.indices.astype(np.int32),matrix.data))
        checked(h.changeColsIntegrality(self.n,np.arange(self.n,dtype=np.int32),np.ones(self.n,dtype=np.uint8)))
        if not strict_proof and use_mip_start:
            solution=np.r_[bits,[bits[a]*bits[b] for a,b in pairs]]
            checked(h.setSolution(len(c),np.arange(len(c),dtype=np.int32),solution))
        # Export canonically scored legal incumbents while the global tree runs.
        if progress and not strict_proof:
            best=[self.objective(bits)]
            def checkpoint(event):
                raw=event.data_out.mip_solution
                if len(raw)<self.n:return
                v=np.rint(np.asarray(raw[:self.n]))
                value=self.objective(v)
                if value<best[0]:
                    candidate=self.decode(v,incumbent)
                    best[0]=value;progress(candidate)
            h.cbMipImprovingSolution += checkpoint
        checked(h.setOptionValue('time_limit',max(.001,seconds-(time.perf_counter()-build_started))))
        api_preparation_seconds=time.perf_counter()-build_started
        solver_started=time.perf_counter()
        checked(h.run())
        solver_seconds=time.perf_counter()-solver_started
        info=h.getInfo();status=h.getModelStatus()
        if strict_proof:
            solution=h.getSolution()
            counterexample=None
            if solution.value_valid:
                candidate=self.decode(np.asarray(solution.col_value[:self.n]),incumbent)
                counterexample=int(round(self.objective(self.encode(candidate))))
            return dict(counterexample_crossings=counterexample,no_better_proven=status==highspy.HighsModelStatus.kInfeasible,
                        status=h.modelStatusToString(status),cutoff_crossings=int(round(cutoff+self.constant)),
                        nodes=int(info.mip_node_count),seconds=time.perf_counter()-build_started)
        solution=h.getSolution()
        vector=solution.col_value if solution.value_valid else None
        dual=info.mip_dual_bound if info.valid else None
        return self._finish_milp(incumbent,vector,dual,dict(status=h.modelStatusToString(status),
            nodes=int(info.mip_node_count),variables=len(c),constraints=len(lower),
            api_preparation_seconds=api_preparation_seconds,solver_seconds=solver_seconds,
            linear_model_build_seconds=self.linear_model_build_seconds,
            bound_scope='whole component',solver='highspy/HiGHS',solver_version=h.version(),feasible_mip_start=use_mip_start))

    def sdp(self, incumbent, seconds, seed=1, block_size=32, full_threshold=128):
        """Shared-moment block SDP; all pair costs retained, no bucket minima.

        Full Shor SDP for small models; block PSD plus 3x3 pair PSD for larger
        models. This is a relaxation, not the exact published SDP formulation.
        """
        build_started = time.perf_counter()
        import cvxpy as cp
        n = self.n
        z = cp.Variable(n)
        pairs = list(self.quadratic)
        y = cp.Variable(len(pairs)) if pairs else None
        positions = {p:i for i,p in enumerate(pairs)}
        constant = self.constant + self.linear.sum()/2 + sum(self.quadratic.values())/4
        linear = self.linear.copy()/2
        for (a,b),q in self.quadratic.items():
            linear[a]+=q/4;linear[b]+=q/4
        objective = constant+linear@z
        if pairs:
            objective += np.array([self.quadratic[p]/4 for p in pairs])@y
        # Connected greedy blocks; factor moments outside blocks stay shared
        # through the same global z and 3x3 PSD constraints.
        remaining=set(range(n));blocks=[]
        adjacency={i:set() for i in range(n)}
        for a,b in pairs:
            adjacency[a].add(b);adjacency[b].add(a)
        size=n if n <= full_threshold else block_size
        if n <= full_threshold:
            blocks=[list(range(n))]
            remaining.clear()
        while remaining:
            block={min(remaining)}
            while len(block)<size:
                candidates=set().union(*(adjacency[i] for i in block)) & remaining-block
                if not candidates:break
                block.add(max(candidates,key=lambda i:(len(adjacency[i]&block),-i)))
            remaining-=block;blocks.append(sorted(block))
        constraints=[z>=-1,z<=1]
        owner={v:k for k,block in enumerate(blocks) for v in block}
        for k,block in enumerate(blocks):
            matrix=cp.Variable((len(block)+1,len(block)+1),symmetric=True)
            constraints += [matrix >> 0,cp.diag(matrix)==1,matrix[0,1:]==z[block]]
            pos={v:i+1 for i,v in enumerate(block)}
            for a,b in pairs:
                if owner[a]==k and owner[b]==k:
                    constraints.append(matrix[pos[a],pos[b]]==y[positions[a,b]])
        for (a,b),j in positions.items():
            if owner[a]!=owner[b]:
                constraints.append(cp.bmat([[1,z[a],z[b]],[z[a],1,y[j]],[z[b],y[j],1]]) >> 0)
            # McCormick constraints strengthen pair moments for binary values.
            constraints += [y[j]>=z[a]+z[b]-1,y[j]>=-z[a]-z[b]-1,
                            y[j]<=1+z[a]-z[b],y[j]<=1-z[a]+z[b]]
        if self.triangles:
            a,b,c=np.asarray(self.triangles).T
            constraints += [z[a]+z[b]-z[c]>=-1,z[a]+z[b]-z[c]<=1]
        problem=cp.Problem(cp.Minimize(objective),constraints)
        problem.get_problem_data(cp.SCS)
        compilation_seconds=time.perf_counter()-build_started
        problem.solve(solver='SCS',time_limit_secs=max(.01,seconds-compilation_seconds),
                      max_iters=5000,eps=1e-4,verbose=False)
        # Fractional precedence values need projection onto legal whole-row orders.
        state=incumbent;best=self.objective(self.encode(state))
        if z.value is not None and np.all(np.isfinite(z.value)):
            rng=np.random.default_rng(seed)
            for trial in range(32):
                v=np.asarray(z.value).copy()
                if trial:v+=rng.normal(0,.5,n)
                orders=dict(state.chromosome_order)
                for sp,row in self.rows.items():
                    ranks={r:0.0 for r in row}
                    for (a,b),i in self.orders.items():
                        if a.species_id==sp:
                            ranks[a]+=v[i];ranks[b]-=v[i]
                    replacement=iter(sorted(row,key=lambda r:(ranks[r],r)))
                    orders[sp]=tuple(next(replacement) if r in ranks else r for r in orders[sp])
                signs=dict(state.chromosome_orientation)
                for i,group in enumerate(self.basis.free_flip_groups):
                    for r in group:signs[r]=self.basis.base_assignment[r]*(-1 if v[i]>0 else 1)
                candidate=LayoutState(orders,signs)
                # Canonical scoring, not fractional/rounded SDP objective.
                score=score_crossings(self.fixture,candidate,restrict_component_nodes=self.nodes).crossings
                if score<best:state,best=candidate,score
        return state,dict(status=problem.status,numerical_relaxation_objective=problem.value,
            certified_lower_bound=None,psd_block_sizes=[len(b)+1 for b in blocks],
            crossing_pair_moments=len(pairs),formulation='full Shor' if n<=full_threshold else 'shared-moment block Shor',
            compilation_seconds=compilation_seconds,
            note='Numerical diagnostic only; never used for pruning or optimality claims')


def run_experiment(fixture, starting_state, method, seconds, seed=1, progress=None):
    started=time.perf_counter();deadline=started+seconds
    whole=orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,starting_state,whole)
    state=canonicalize_component_order(fixture,starting_state)
    components,component_of=_component_map(fixture)
    diagnostics=[];lower=0;preparation=0
    for i,nodes in enumerate(components):
        if time.perf_counter()>=deadline:
            diagnostics.append(dict(component=i,status='deadline unstarted'));continue
        tick=time.perf_counter()
        refs=tuple(r for r in fixture.chromosome_refs if component_of[r]==i)
        model=JointModel(fixture,nodes,refs)
        preparation+=time.perf_counter()-tick
        record=dict(component=i,decision_variables=model.n,crossing_factors=model.graph.factor_count,
                    hard_orientation_groups=[[r.label for r in g] for g in model.basis.free_flip_groups])
        budget=max(.001,(deadline-time.perf_counter())/(len(components)-i))
        if method=='milp':
            state,bound,info=model.milp(state,budget);lower+=bound;record.update(info)
        elif method=='lns':
            # Fixing outside a temporary neighborhood is never exported as a
            # global reduction. Every move remains within the shared hard basis.
            stop=time.perf_counter()+budget;rounds=0
            while time.perf_counter()<stop:
                center=rounds % len(fixture.species_ids)
                active=set(fixture.species_ids[max(0,center-1):center+2])
                bits=model.encode(state)
                free={j for (a,b),j in model.orders.items() if a.species_id in active}
                free.update(j for j,g in enumerate(model.basis.free_flip_groups)
                            if any(r.species_id in active for r in g))
                fixed={j:bits[j] for j in range(model.n) if j not in free}
                state,_,info=model.milp(state,min(10,stop-time.perf_counter()),fixed)
                rounds+=1
                if progress:progress(state)
                if info['status']==0 and not fixed:break
            record.update(rounds=rounds,bound_scope='neighborhood only; global bound not exported')
        elif method=='sdp':
            state,info=model.sdp(state,budget,seed);record.update(info)
        else:raise ValueError(method)
        diagnostics.append(record)
        if progress:progress(state)
    validate_saved_layout(fixture,state,whole)
    upper=score_crossings(fixture,state).crossings
    if upper>score_crossings(fixture,starting_state).crossings:raise AssertionError('Incumbent lost')
    return dict(optimized_state=state.to_dict(),upper_bound=upper,lower_bound=lower,
        optimality_gap=upper-lower,optimality_status='proven optimum' if upper==lower else 'bounded best known',
        seconds=time.perf_counter()-started,preparation_seconds=preparation,method=method,
        component_diagnostics=diagnostics,starting_state=starting_state.to_dict(),
        whole_chromosome_order_changes={sp:[r.chromosome_id for r in state.chromosome_order[sp]]
            for sp in fixture.species_ids if state.chromosome_order[sp]!=starting_state.chromosome_order[sp]},
        whole_chromosome_flips=[r.label for r in fixture.chromosome_refs
            if state.chromosome_orientation[r]!=starting_state.chromosome_orientation[r]],
        reduction_policy='shared graph/GF2 preprocessing; independent experiment, no search-ledger restart')
