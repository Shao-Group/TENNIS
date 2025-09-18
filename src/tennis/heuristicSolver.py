"""
BSD 3-Clause License

Copyright (c) 2025, Xiaofei Carl Zang, Ke Chen, Mingfu Shao, and The Pennsylvania State University

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

import numpy as np
from typing import List, Tuple, Dict
import itertools

try:
    from .bound import boundComputer
    from .util import *
except ImportError:
    # For standalone execution
    from bound import boundComputer
    from util import *


class HeuristicSolver:
    def __init__(self, matrix: List[List[int]], bc: boundComputer = None) -> None:
        self.matrix = np.array(matrix)
        self.bound_computer = bc if bc is not None else boundComputer(self.matrix)
        self.connected_components = self.bound_computer.get_connected_components()
        self.num_components = len(self.connected_components)
        
        self.is_feasible = False
        self.min_additional_nodes = -1
        self.all_solutions = []
        self.solution_info = []
        
        self._solve()
    
    def _solve(self) -> None:
        """Main solving logic based on number of connected components."""
        if self.num_components <= 1:
            # Already connected, no additional nodes needed
            self.is_feasible = True
            self.min_additional_nodes = 0
            self.all_solutions = []
            self.solution_info = []
        elif self.num_components == 2:
            # Use 2-CC heuristic
            self._solve_two_components()
        else:
            # For now, only handle 2 CCs
            self.is_feasible = False
            self.min_additional_nodes = -1
    
    def _solve_two_components(self) -> None:
        """Solve the case with exactly 2 connected components."""
        cc1, cc2 = self.connected_components[0], self.connected_components[1]
        
        # Find minimum hamming distance between the two components
        min_distance = float('inf')
        best_pairs = []
        
        # Vectorized computation of all pairwise distances
        matrix_cc1 = self.matrix[cc1]  # Shape: (len(cc1), num_features)
        matrix_cc2 = self.matrix[cc2]  # Shape: (len(cc2), num_features)

        # Broadcast and compute all pairwise distances at once
        distances = np.sum(matrix_cc1[:, np.newaxis] != matrix_cc2[np.newaxis, :], axis=2)

        # Find minimum distance and all pairs with that distance
        min_distance = np.min(distances)
        min_indices = np.where(distances == min_distance)
        best_pairs = [(cc1[i], cc2[j]) for i, j in zip(min_indices[0], min_indices[1])]
        
        # Number of additional nodes needed = min_distance - 1
        self.min_additional_nodes = max(0, min_distance - 1)
        self.is_feasible = True
        
        if self.min_additional_nodes == 0:
            # Components are already distance 1 apart, no intermediate nodes needed
            self.all_solutions = []
            self.solution_info = []
        else:
            # Generate all possible minimum solutions
            self._generate_all_paths(best_pairs, min_distance)
    
    def _generate_all_paths(self, best_pairs: List[Tuple[int, int]], distance: int) -> None:
        """Generate all possible paths between the two components."""
        all_solutions = []
        all_info = []
        
        for node1_idx, node2_idx in best_pairs:
            vec1 = self.matrix[node1_idx]
            vec2 = self.matrix[node2_idx]
            
            # Find all positions where vectors differ
            diff_positions = np.where(vec1 != vec2)[0]
            assert len(diff_positions) == distance
            
            # Generate all possible intermediate paths
            paths = self._enumerate_paths(vec1, vec2, diff_positions)
            
            for path in paths:
                all_solutions.append(path)
                info = {
                    'start_node': node1_idx,
                    'end_node': node2_idx,
                    'start_vector': vec1.tolist(),
                    'end_vector': vec2.tolist(),
                    'path_length': distance,
                    'intermediate_nodes': len(path)
                }
                all_info.append(info)
        
        self.all_solutions = all_solutions
        self.solution_info = all_info
    
    def _enumerate_paths(self, start_vec: np.ndarray, end_vec: np.ndarray, 
                        diff_positions: np.ndarray) -> List[List[List[int]]]:
        """Enumerate all possible intermediate paths between two vectors."""
        if len(diff_positions) <= 1:
            return []  # No intermediate nodes needed
        
        paths = []
        
        # Generate all possible orderings of bit flips
        for flip_order in itertools.permutations(diff_positions):
            path = []
            current_vec = start_vec.copy()
            
            # Apply flips one by one to create intermediate nodes
            for i in range(len(flip_order) - 1):  # Exclude the last flip (goes to end_vec)
                pos = flip_order[i]
                current_vec[pos] = end_vec[pos]  # Flip to target value
                path.append(current_vec.copy().tolist())
            
            paths.append(path)
        
        return paths
    
    def get_is_feasible(self) -> bool:
        """Return whether the problem is feasible."""
        return self.is_feasible
    
    def get_min_additional_nodes(self) -> int:
        """Return minimum number of additional nodes needed."""
        return self.min_additional_nodes
    
    def get_all_solutions(self) -> List[List[List[int]]]:
        """Return all possible minimum solutions."""
        return self.all_solutions
    
    def get_solution_info(self) -> List[Dict]:
        """Return detailed information about each solution."""
        return self.solution_info
    
    def get_num_solutions(self) -> int:
        """Return total number of solutions found."""
        return len(self.all_solutions)


def test_heuristic_solver():
    """Test function for the heuristic solver."""
    # Test case 1: Two components with distance 3
    print("Test Matrix:")
    test_matrix = [
        [1, 0, 0, 0],  # CC1
        [1, 1, 0, 0],  # CC1 (connected to first by 1 bit)
        [0, 0, 1, 1],  # CC2 (distance 3 from CC1)
        [0, 1, 1, 1],  # CC2 (connected to third by 1 bit)
    ]
    
    solver = HeuristicSolver(test_matrix)
    
    print("Test Matrix:")
    for i, row in enumerate(test_matrix):
        print(f"  Vector {i}: {row}")
    
    print(f"\nNumber of connected components: {solver.num_components}")
    print(f"Connected components: {solver.connected_components}")
    print(f"Is feasible: {solver.get_is_feasible()}")
    print(f"Minimum additional nodes: {solver.get_min_additional_nodes()}")
    print(f"Number of solutions: {solver.get_num_solutions()}")
    
    if solver.get_num_solutions() > 0:
        print("\nAll solutions:")
        for i, (solution, info) in enumerate(zip(solver.get_all_solutions(), solver.get_solution_info())):
            print(f"  Solution {i+1}:")
            print(f"    Connects vector {info['start_node']} to vector {info['end_node']}")
            print(f"    Intermediate nodes: {solution}")


if __name__ == "__main__":
    test_heuristic_solver()