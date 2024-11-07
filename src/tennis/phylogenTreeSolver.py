"""
BSD 3-Clause License

Copyright (c) 2024, Xiaofei Carl Zang, Ke Chen, Mingfu Shao, and The Pennsylvania State University

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

from .util import *
from pysat.card import EncType, CardEnc
from pysat.solvers import Solver
from pysat.formula import Atom, IDPool, CNF, PYSAT_TRUE, PYSAT_FALSE #, Neg
import sys
import time
import random
import itertools
from typing import List
from threading import Timer

class PhylogenTreeSolver():
    def __init__(self, knownTx: List[List[int]],
                 maxAddlNodes:int = 4,
                 maxEdged: int = 1,
                 solver: str = '',
                 formulation: str = 'SATSimple',
                 time_limit = 900,
                 printAllVar: bool = False) -> None:

        self.knownTx = knownTx              # a 2D binary matrix for known transcripts
        self.maxEdgeDist = maxEdged 
        self.numKnownTx = len(knownTx)
        self.numExons = len(knownTx[0])
        self.minAddlNodes = -1
        self.maxAddlNodes = maxAddlNodes
        self.solver = solver 
        self.novelTx = []
        self.novelTxInfo = []
        self.formulation = formulation
        self.printAllVar = printAllVar
        self.timed_out = False
        self.time_limit = time_limit
        if solver == '':
            self.solver = 'Glucose4' if (formulation == 'SAT' or
                                         formulation == 'SATSimple') else 'ERROR'
        self.__sanity_check()
        self.__get_maxAddlNodes()
        self.__satPhyloTree()
    
    def is_timed_out(self) -> bool:
        return self.timed_out

    def get_novelTx_and_info(self) -> List[List[int]]:
        if not self.is_feasible():
            assert self.novelTx == []
        if self.minAddlNodes == -1:
            assert self.novelTx == []
        assert len(self.novelTx) == len(self.novelTxInfo)
        return self.novelTx, self.novelTxInfo

    def get_minAddlNodes(self) -> int:
        return self.minAddlNodes 
    
    def is_feasible(self) -> bool:
        return self.minAddlNodes != -1

    # Input: multi sols as 3D array
    # Populate self.novelTx and self.novelTxInfo by flattening novel tx info
    # return: number of non-dup groups
    def __remove_dup_novel_solution_and_get_mult_sol_info(self, multi_sols) -> int:
        if len(multi_sols) == 0:
            return 0
        if len(multi_sols[0]) == 0:
            return 0
        addlNodes = len(multi_sols[0])
        UnDupNovelTx = dict()   # 2D chains of transcripts, not multi groups
        # UnDupNovelTxHashs = set()
        groupHashs = set() 
        transHashs2groups = dict()
        for novelGroup in multi_sols:
            assert addlNodes == len(novelGroup)

            TxHashs = [hash_as_str(tx) for tx in novelGroup]
            TxHashs.sort()
            groupHash = hash_as_str(TxHashs)
            # whether whole group is dup
            if groupHash in groupHashs:
                continue
            groupHashs.add(groupHash)

            # record tx to non-dup group
            group_index = len(groupHashs)   # because some groups are removed due to whole-group-duplication
            # for h in TxHashs:
            for tx in novelGroup:
                h = hash_as_str(tx)
                assert h in TxHashs
                if h not in transHashs2groups:
                    assert h not in UnDupNovelTx
                    UnDupNovelTx[h] = tx
                    transHashs2groups[h] = [group_index]
                else:
                    assert h in UnDupNovelTx
                    transHashs2groups[h].append(group_index)        
        # store info
        assert len(self.novelTx)     == 0
        assert len(self.novelTxInfo) == 0
        for h, tx in UnDupNovelTx.items():
            self.novelTx.append(tx)
            groups = transHashs2groups[h]
            info = dict()
            info['solutions_in'] = groups
            info['solutions_total_num'] = len(groupHashs)
            info['solution_unique'] = (len(groupHashs) == 1)
            info['mandatory'] = (len(groups) == len(groupHashs))
            info['PctIn'] = (len(groups)/len(groupHashs))
            info['timed_out'] = self.timed_out
            self.novelTxInfo.append(info)
        assert(len(self.novelTx)) == len(transHashs2groups)
        return len(groupHashs)

    def __sanity_check(self):
        assert self.numKnownTx >= 1
        for tx in self.knownTx:
            assert len(tx) == self.numExons
            for x in tx:
                assert isinstance(x, int)
                assert is_binary(x)
        # check formulation
        assert self.formulation in ['SAT', 'SATSimple']
        return 0

    def __get_maxAddlNodes(self) -> int:
        b = [sum(i) for i in self.knownTx]
        x = self.numExons * self.numKnownTx - sum(b)
        self.maxAddlNodes = min(x, self.maxAddlNodes)
        return self.maxAddlNodes

    # return `minNum`` needed missing internal nodes to construct a tree
    def __satPhyloTree(self) -> int:      
        for i in range(self.maxAddlNodes + 1):
            # print("Current add'l nodes: ", i)
            if self.formulation == 'SATSimple':
                is_solved = self.__ilpPhyloTreeSatSimple(i)
                if self.timed_out and not is_solved: 
                    break    # time out and not solved
                elif is_solved:
                    self.minAddlNodes = i
                    # print(f"Found solution with {i} add'l nodes")
                    return i
                else:
                    continue
            else:
                raise RuntimeError("model not implemented")
        print("Cannot contruct phylo tree! ")
        if self.timed_out:
            print(f"timed out when trying {i} nodes")
        return -1

    # where a phylo tree can be constructed with distance-maxEdgeDist edges
    # using addlN missing internal nodes
    # this version does not consider alternative TSS/TES, namely, all bits participate in the distance calculation
    def __ilpPhyloTreeSatSimple(self, addlN) -> bool:
        startTime = time.time()
        # use when the atoms in the list v maybe constant
        def nameOfAtom(x:Atom) -> int:
            if x == PYSAT_TRUE:
                return 1
            elif x == PYSAT_FALSE:
                return 0
            else:
                return x.name
            
        def Neg(x:Atom) -> int:
            if x == PYSAT_TRUE:
                return 0
            elif x == PYSAT_FALSE:
                return 1
            else:
                return -x.name

        # this version assumes the Atoms/Negs are already processed into ints (by the above two functions),
        # 0/1 are reserved for PYSAT_FALSE/TRUE respectively   
        def listOr(v:list) -> list:
            result = [0]*len(v)
            resultLen = 0
            for x in v:
                assert(x != -1)
                if x == 1:
                    return [] # [] is treated as false in the solver but we will remove them from the final cnf
                elif x == 0:
                    continue
                else:
                    result[resultLen] = x
                    resultLen += 1
            return result[:resultLen]

        
        N = self.numKnownTx + addlN    # total num nodes
        assert(N >= 1 and self.numExons >= 1)
        minExonNum = min([sum(tx) for tx in self.knownTx])
        assert(minExonNum >= 1)
        assert self.minAddlNodes <= 0   # no sol found yet
        idPool = IDPool(start_from = 2) # 0 and 1 reserved for PYSAT_FALSE and PYSAT_TRUE

        treeSatisfy = []
        
        # vN: all known and unknown tx, represented as binary vector
        vN = [[PYSAT_TRUE if x==1 else PYSAT_FALSE for x in row] for row in self.knownTx]
        for i in range(addlN):
            vN.append([Atom(idPool.id(f'v_{i+self.numKnownTx}_{k}')) for k in range(self.numExons)])
        assert len(vN) == N
        for x in vN:
            assert len(x) == self.numExons

        addlNRange = range(self.numKnownTx, N)
        
        CONSIDER_MULTI_PARTIAL_INTRON = True
        Eq = [[None]*self.numExons for _ in range(N)]
        if CONSIDER_MULTI_PARTIAL_INTRON:
            #helper variables Eq_ik == True iff vN[i][k] == vN[i][k-1]
            for i in range(N):
                for k in range(self.numExons):
                    if k <= 0:
                        Eq[i][k] = PYSAT_FALSE
                    elif i < self.numKnownTx: # and k >= 1
                        Eq[i][k] = PYSAT_TRUE if vN[i][k] == vN[i][k-1] else PYSAT_FALSE
                    else: # i >= self.numKnownTx and k >= 1
                        assert k >= 1
                        Eq[i][k] = Atom(idPool.id(f'Eq_{i}_{k}'))
                        # if vik==vjk, then Eqik must be True
                        treeSatisfy.extend([listOr([nameOfAtom(vN[i][k]), nameOfAtom(vN[i][k - 1]), nameOfAtom(Eq[i][k])]),
                                            listOr([Neg(vN[i][k]), Neg(vN[i][k - 1]), nameOfAtom(Eq[i][k])])])
                        # if vik!=vjk, then Eqik must be False
                        treeSatisfy.extend([listOr([nameOfAtom(vN[i][k]), Neg(vN[i][k - 1]), Neg(Eq[i][k])]),
                                            listOr([Neg(vN[i][k]), nameOfAtom(vN[i][k - 1]), Neg(Eq[i][k])])])
        else:
            # should not impact following computations
            Eq = [[PYSAT_FALSE]*self.numExons for _ in range(N)]
        for Eqi in Eq:
            for Eqik in Eqi:
                assert Eqik is not None

        # d[i,j,k]: binary, 0 if pos k of transcripts i and j agree, 1 otherwise
        # namely, d[i,j,k] = vN[i,k]!=vN[j,k] and (not Eq[i][k] or not Eq[j][k])
        # D[i,j]: integer, distance between transcripts i and j,
        # namely, D[i,j] = lpSum(d[i,j]), not used in SAT formulation (only computed between known transcripts)

        d = [[[None]*self.numExons for _ in range(N)] for _ in range(N)]
        # only calculate for the known transcripts
        D = [[0]*self.numKnownTx for _ in range(self.numKnownTx)]                          

        for i in range(N):
            for j in range(N):
                for k in range(self.numExons):
                    if i < self.numKnownTx and j < self.numKnownTx:
                        if vN[i][k] != vN[j][k] and (Eq[i][k] == PYSAT_FALSE or Eq[j][k] == PYSAT_FALSE):
                            d[i][j][k] = PYSAT_TRUE
                            D[i][j] += 1
                        else:
                            assert vN[i][k] == vN[j][k] or (Eq[i][k] == PYSAT_TRUE and Eq[j][k] == PYSAT_TRUE)
                            d[i][j][k] = PYSAT_FALSE
                    elif i < j:
                        d[i][j][k] = Atom(idPool.id(f'd_{i}_{j}_{k}'))
                        # if vik==vjk, then dijk must be false
                        treeSatisfy.extend([listOr([nameOfAtom(vN[i][k]), nameOfAtom(vN[j][k]), Neg(d[i][j][k])]),
                                            listOr([Neg(vN[i][k]), Neg(vN[j][k]), Neg(d[i][j][k])])])
                        # if Eq_ik and Eq_jk, then dijk must be false
                        treeSatisfy.extend([listOr([Neg(Eq[i][k]), Neg(Eq[j][k]), Neg(d[i][j][k])])])
                        # if vik!=vjk and (not Eq_ik or not Eq_jk), then dijk must be true
                        treeSatisfy.extend([
                                            listOr([nameOfAtom(vN[i][k]), Neg(vN[j][k]), nameOfAtom(d[i][j][k]), nameOfAtom(Eq[i][k])]),
                                            listOr([nameOfAtom(vN[i][k]), Neg(vN[j][k]), nameOfAtom(d[i][j][k]), nameOfAtom(Eq[j][k])]),
                                            listOr([Neg(vN[i][k]), nameOfAtom(vN[j][k]), nameOfAtom(d[i][j][k]), nameOfAtom(Eq[i][k])]),
                                            listOr([Neg(vN[i][k]), nameOfAtom(vN[j][k]), nameOfAtom(d[i][j][k]), nameOfAtom(Eq[j][k])])
                                            ])
                    elif i == j:
                        d[i][j][k] = PYSAT_FALSE
                    else:
                        d[i][j][k] = d[j][i][k]


        # Helper variables:
        # D_eq_1[i,j]: bianry, 1 iff D[i,j] == 1
        D_eq_1 = [[None]*N for _ in range(N)]

        for i in range(N):
            for j in range(N):
                if i < self.numKnownTx and j < self.numKnownTx:
                    D_eq_1[i][j] = PYSAT_TRUE if D[i][j] == 1 else PYSAT_FALSE
                elif i == j:
                    D_eq_1[i][j] = PYSAT_FALSE
                elif i < j:
                    dijVarIds = [x.name for x in d[i][j]]

                    D_eq_1[i][j] = Atom(idPool.id(f'Deq1_{i}_{j}'))
                    # if sum(d[i][j][k]) is not 1, then Deq1[i][j] must be false
                    treeSatisfy.append([-D_eq_1[i][j].name] + dijVarIds) # if sum is 0
                    treeSatisfy.extend([[-D_eq_1[i][j].name, -x, -y] for x,y in itertools.combinations(dijVarIds,2)]) # if sum >= 2
                    # if sum(d[i][j][k]) is 1, then Deq1[i][j] must be true
                    treeSatisfy.extend([[D_eq_1[i][j].name, -dijVarIds[k]] + dijVarIds[:k] + dijVarIds[k+1:] for k in range(self.numExons)])
                else: # i > j, can borrow existing variables
                    D_eq_1[i][j] = D_eq_1[j][i]
                

        # T[i,k]: whether node i on level k
        T = [[Atom(idPool.id(f'T_{i}_{k}')) for k in range(N)] for i in range(N)]

        # Constraint 1: A node appears exactly once
        for i in range(N):
            treeSatisfy.extend(CardEnc.equals(lits=[x.name for x in T[i]], bound=1, encoding=EncType.pairwise).clauses)

        # Constraint 2: if node i is on level g (g>=1), then on level g-1 exists node j with D_eq_1[j][i] = 1
        # C = [[[Atom(idPool.id(f'C_{i}_{j}_{g1}')) for g1 in range(N - 1)] for j in range(N)] for i in range(N)]
        C = [[[None]*(N-1) for _ in range(N)] for _ in range(N)]
        for i in range(N):
            for j in range(N):
                if i == j or D_eq_1[j][i] == PYSAT_FALSE:
                    for g1 in range(N - 1):
                        C[i][j][g1] = PYSAT_FALSE
                else:
                    for g1 in range(N - 1):
                        C[i][j][g1] = Atom(idPool.id(f'C_{i}_{j}_{g1}'))
                        # if i is not on g1+1 or j is not on g1, C[i][j][g1] is trivially false
                        treeSatisfy.extend([[-C[i][j][g1].name, T[i][g1+1].name],
                                            [-C[i][j][g1].name, T[j][g1].name]])
                    # if D_eq_1[j][i] is false, then C[i][j][k] must be false
                    if D_eq_1[j][i] != PYSAT_TRUE:
                        treeSatisfy.extend([[-C[i][j][g1].name, D_eq_1[j][i].name] for g1 in range(N - 1)])                    

        for i in range(N):
            level1_or_g = [cijk.name for cij in C[i] for cijk in cij if cijk != PYSAT_FALSE]
            level1_or_g.append(T[i][0].name)        
            treeSatisfy.extend(CardEnc.equals(lits=level1_or_g, bound=1, encoding=EncType.pairwise).clauses)

        # Constraint 3: Must have only one node on root level (enforcing a tree, not a forest)
        treeSatisfy.extend(CardEnc.equals(lits=[T[i][0].name for i in range(N)], bound=1, encoding=EncType.pairwise).clauses)
        
        curTime = time.time()
        # print(f'Building the model takes {curTime - startTime:.3f} seconds')

        # solve
        startTime = time.time()
        
        satSolver = None
        
        if self.solver == 'Glucose3':
            satSolver = Solver(name = 'g3')
        elif self.solver == 'Glucose4':
            satSolver = Solver(name = 'g4')
        else:
            print('Cannot recognize solver or solver is not available, proceed with Glucose4', file=sys.stderr)
            satSolver = Solver(name = 'g4')

        # may have added [] as a clause which is treated as false, need to remove them before passing to solver
        satSolver.append_formula([c for c in treeSatisfy if len(c) > 0])

        timed_out = False
        
        def sat_interrupt(s: Solver):
            s.interrupt()
            timed_out = True
            
        timer = Timer(self.time_limit, sat_interrupt, [satSolver])

        satSolver.start_mode(warm=True)
        multNovelTx = []
        timer.start()
        sol_counter = 0
        while satSolver.solve_limited(expect_interrupt=True) == True:
            # gather solution
            sol_counter += 1
            truthAssignment = satSolver.get_model()
            currentSolution = []
            for i in range(addlN):
                binary_tst = [1 if truthAssignment[x.name - 1] > 0 else 0 for x in vN[self.numKnownTx + i]]
                currentSolution.append(binary_tst)
            multNovelTx.append(currentSolution)

            # ban the current solution and ask solver to solve again
            satSolver.add_clause([-truthAssignment[x.name - 1] for i in range(addlN) for x in vN[self.numKnownTx + i]])

            if self.printAllVar:
                print(f'solution {len(multNovelTx)}')
                for i in range(1, idPool.top):
                    print(f'{idPool.obj(i+1)}: {1 if truthAssignment[i] > 0 else 0}')
            if addlN == 0:
                break   # only care when adding new isoforms
        timer.cancel()
                
        curTime = time.time()
        self.timed_out = timed_out
        if sol_counter > 0:
            print(f"\tUsing {addlN} additional node(s)")
            print(f"\tSAT found {sol_counter} solution(s)")
            print(f"\tSolver took {curTime - startTime:.3f} seconds. Timed out? {timed_out}")

        satSolver.delete()

        if sol_counter > 0:
            assert addlN == 0 or sol_counter == len(multNovelTx)
        else:
            assert len(multNovelTx) == 0

        self.__remove_dup_novel_solution_and_get_mult_sol_info(multNovelTx)
        return sol_counter > 0



def test(cols = 20, rows = 10):
    while True:
        randomBinrayMatrix = [[random.randint(0, 1) for _ in range(cols)] for _ in range(rows)]
        y = randomBinrayMatrix
        y = [x for x in randomBinrayMatrix if sum(x) >= cols - 4]
        if len(y) <= 0:
            continue
        if len(y[0]) <= 0:
            continue
        break
    print('simulated matrix done')
    for __num_tests in range(1):    
        x = PhylogenTreeSolver(y)
        print('Success ilpTree!')
        print(f"test #{__num_tests}, addlNodes {x.get_minAddlNodes()}")
    return 0

def test_array(formulation: str = 'ILPOriginal'):
    # variables names must be different
    x = [[1, 1, 0], [1, 1, 1], [0, 0, 1], [0, 1, 0]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array1', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[1 ,0 ,0], [0, 1, 0]] 
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array2', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[1 ,0 ,0], [0, 1, 1]] 
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)    
    print('test_array3', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[0 ,1 ,0], [1, 0, 1], [1, 1, 1]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array4', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[0 ,1 ,0, 0], [1, 1, 1, 1]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array5', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[0 ,1 ,0], [1, 0, 1]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array6', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[1 ,0 ,0], [1, 1, 1]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array7', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)

    x = [[0,0,1,0,1,1,0,1,0,0,0],
         [0,0,0,1,0,1,0,1,0,0,0],
         [0,0,0,1,1,1,1,0,0,0,0]]
    y = PhylogenTreeSolver(x, formulation = formulation, printAllVar = True)
    print('test_array8', y.is_feasible(), y.get_minAddlNodes())
    print(y.novelTx)
    
    return 0


if __name__ == "__main__":
    raise RuntimeError("This module is not meant to be run as a script")
