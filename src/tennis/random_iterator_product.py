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

import itertools
import random
from copy import copy
import sys

# similar to itertools.product(l, repeat=n), but randomized
# works better when n is very big
class random_self_cartesian_product:
    def __init__(self, l, repeat):
        self.l = l
        self.r = repeat
        self.recur_iterable = []
        self.non_empty_iterable = list(range(len(self.l)))
        self.y = self.random_yield(self.l, self.r)

    def random_yield(self, list_to_randomize, remaining_time: int):
        if self.recur_iterable == []:
            self.recur_iterable = [random_self_cartesian_product(self.l, self.r - 1) for i in range(len(self.l))] if self.r >= 2 else []
        assert remaining_time >= 1
        if remaining_time <= 1:
            list_len = len(list_to_randomize)
            for i in random.sample(list_to_randomize, list_len):
                # print('l, r, yield', self.l, self.r, [i] )
                yield [i]
        else:
            while len(self.non_empty_iterable) >= 1:
                z = random.choice(self.non_empty_iterable)
                is_yielded = False
                x = self.recur_iterable[z].__peek()
                if x is None:
                    self.non_empty_iterable.remove(z)
                    continue
                x.append(self.l[z])                
                yield x
        return None

    def __iter__(self):
        return self.y

    def __peek(self):
        if self.recur_iterable == []:
            self.recur_iterable = [random_self_cartesian_product(self.l, self.r - 1) for i in range(len(self.l))] if self.r >= 2 else []
        try:
            return next(self.y)
        except StopIteration:
            return None
    
    def peek(self):
        x = self.__peek()
        if x is None:
            raise StopIteration
        else:
            return x
        

def produce_random(l, r):
    random.seed(2024)
    recursion_limit = sys.getrecursionlimit()
    if 10 * r >= recursion_limit:
        sys.setrecursionlimit(10 * r) 
    y = random_self_cartesian_product(l, r)
    while True:
        try:
            x = y.peek()
            yield x
        except StopIteration:
            break
    sys.setrecursionlimit(recursion_limit)
    return None


if __name__ == "__main__":
    raise RuntimeError("This module is not supposed to be executed as a standalone script.")