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

from pysat.formula import Atom, PYSAT_TRUE

# input: exon list of tuple
# must be inclusive, non-overlapping, and sorted
def exon_list_valid(exs):
    assert len(exs) >= 1 
    for i in range(len(exs)):
        assert len(exs[i]) == 2  # each exon is expressed as a tuple
        assert exs[i][0] <= exs[i][1] 
    
    for i in range(len(exs) - 1):
        ex1 = exs[i][1]
        ex2 = exs[i + 1][0]
        assert ex1 < ex2
    return None

def is_binary(x) -> bool:
    return is_0(x) or is_1(x)

def is_0(x) -> bool:
    return x < 0.00001 and x > -0.00001

def is_1(x) -> bool:
    return x > 0.99999 and x < 1.00001

# index of the first nonzero element in a list in the forward direction
# or the backward direction 
def first_non_zero(x:list, forwardDir:bool=True) -> int:
    itr = range(len(x)) if forwardDir else range(len(x)-1, -1, -1)
    for i in itr:
        if isinstance(x[i], Atom):
            if x[i] == PYSAT_TRUE:
                return i
        elif not is_0(x[i]):
            return i
    return -1

def hash_as_str(x):
    return hash(str(x))
