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

from sys import argv
from GTF import parse as parse_gtf_line
from os.path import basename, exists
from collections import defaultdict


def quick_scan(gtf):
    # gid_span = {gid:(first appearance, last appearance)}
    # we split file at lines that is not inside any gid_span
    gid_span = dict() 

    f = open(gtf, 'r')
    line_counter = 0
    splits = []
    for line in f.readlines():
        if line.startswith('#'): continue
        line_counter += 1
        fields = parse_gtf_line(line)
        gid = fields['gene_id']
        if gid not in gid_span:
            gid_span[gid] = (line_counter, line_counter)
        else:
            s, t = gid_span[gid]
            gid_span[gid] = (min(s, line_counter), max(t, line_counter))
    
    sorted_span = sorted(list(gid_span.values()))
    sorted_start = sorted([x[0] for x in list(gid_span.values())])
    sorted_end   = sorted([x[1] for x in list(gid_span.values())], reverse=True)
    
    GENE_NUM_PER_FILE = 500
    for x in range(1, len(sorted_span), GENE_NUM_PER_FILE):
        s,t = sorted_span[x]
        assert s == sorted_start[x]
        if t == sorted_end[x]:
            splits.append(t)
        elif t < sorted_end[x]:
            print('previous gene spans into this one')
        else:
            print('this gene spans into later ones')
    f.close()
    return splits


# assuming gtf is sorted by gene ids
def split(gtf):
    global file_index
    line_batch = []
    gene_ids_permanent = set()
    gene_ids_batch = set()
    GENE_NUM_PER_FILE = 500

    splits = quick_scan(gtf)

    f = open(gtf, 'r')
    line_counter = 0
    split_counter = 0
    for line in f.readlines():
        if line.startswith('#'): continue
        line_counter += 1
        if line_counter > splits[split_counter]:
            # flush out previous genes
            assert (line_counter - 1 ==splits[split_counter])
            split_counter += 1
            write_file(line_batch)
            line_batch = []
            gene_ids_batch = set()
            file_index += 1
        fields = parse_gtf_line(line)
        gid = fields['gene_id']
        # new gid
        if gid not in gene_ids_batch:  
            # record new genes
            if gid in gene_ids_permanent:
                print(gid)
            assert  gid not in gene_ids_permanent# otherwise gtf is not sorted, a gene might be in different files
            gene_ids_permanent.add(gid)
            gene_ids_batch.add(gid)
        line_batch.append(line)
    write_file(line_batch)
    f.close()

def write_file(lines):
    global file_index
    global output_prefix
    if (file_index % 10 == 0 ): 
        print(file_index)
    with open(output_prefix + "-" + str(file_index) + ".gtf",'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    gtf = argv[1]
    print(gtf)
    # some dangerous GLOBAL variables
    file_index = 0
    output_prefix = basename(gtf) 
    output_prefix = (argv[2]) if len(argv) >= 3 else basename(gtf)
    print(output_prefix)
    split(gtf)
