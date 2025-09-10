#!/usr/bin/env python

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

from .random_iterator_product import produce_random
from .phylogenTreeSolver import PhylogenTreeSolver
from .GTF import parse as parse_gtf_line
from .GTF import get_xi_counts
from .util import *
from typing import List
from collections import defaultdict
from collections import Counter
from sys import argv
from pprint import pprint, pformat
import pickle
from os.path import basename, exists
import sys
import os.path
from time import sleep 
from datetime import datetime
import random
import argparse

class GeneChainToTree():
    # given chains of a gene (2D matrix)
    # make binary representation
    # then construct tree
    def __init__(self, chains, maxAddlNodes, formulation:str = 'SATSimple', args = None):
        self.chains = chains
        self.__maxAddlNodes = maxAddlNodes
        self.formulation = formulation
        self.args = args
        self.__exons = []
        self.__binaries = []
        self.__trivial_cols = dict() # col: val
        self.__get_exons(with_introns = True)
        self.__chains_to_binaries()

        self.__timed_out = False
        self.__is_feasible  = False
        self.__minAddlNodes = -1
        self.__novel_binaries = [] # multiple solutions (3D) of 2D array
        self.__novel_info = []
        if formulation == "RandomAll" or formulation == "Random1":
            self.__random_sol(1)
        elif formulation == "RandomX":
            self.__random_sol(self.__maxAddlNodes)
        else:
            self.__contruct_tree()
        assert self.__is_feasible == (self.__minAddlNodes >= 0)

    def is_feasible(self):
        return self.__is_feasible
    
    def get_minAddlNodes(self):
        return self.__minAddlNodes

    # use novel binary matrix to return novel feature chains
    def get_novel_transcripts(self):
        tsts = [[]] * len(self.__novel_binaries)
        for i, tx_binary in enumerate(self.__novel_binaries):
            tsts[i] = [self.__exons[j] for j in range(len(self.__exons)) if is_1(tx_binary[j])]
        if not self.__is_feasible:
            assert self.__minAddlNodes == -1
            assert len(tsts) == 0
        return tsts
    
    def get_tx_info(self):
        return self.__novel_info
    
    def delivery(self):
        print('print self binaries')
        pprint(self.__binaries)

    def __get_exons(self, with_introns: bool = True):
        exons = set()
        for chain in self.chains:
            for exon in chain:
                exons.add(exon)
        exons = sorted(list(exons)) # collection of all exons
        if with_introns:
            introns = []
            for i in range(len(exons) - 1):
                ss5 = exons[i][1]
                ss3 = exons[i + 1][0]
                if ss5 + 1 >= ss3:
                    continue    # consecutive exons
                else:
                    assert ss5 + 1 <= ss3 - 1
                    introns.append((ss5 + 1, ss3 - 1))    # make regions non-overlapping
            exons.extend(introns)
            exons = sorted(exons)
        self.__exons = exons
        return 0

    def __chains_to_binaries(self):
        binaries = []
        for c in self.chains:
            cbinary = [1 if e in c else 0 for e in self.__exons]
            binaries.append(cbinary)
        self.__binaries = binaries
        return 0
    
    def __contruct_tree(self):
        bmatrix = self.__binaries
        treeSolver = PhylogenTreeSolver(bmatrix, self.__maxAddlNodes, formulation=self.formulation, time_limit=self.args.time_out)
        self.__is_feasible  = treeSolver.is_feasible() 
        self.__minAddlNodes = treeSolver.get_minAddlNodes()
        self.__timed_out = treeSolver.is_timed_out()
        self.__novel_binaries, self.__novel_info = treeSolver.get_novelTx_and_info()
        return 0
    
    def __random_sol(self, sol_num):
        bmatrix = self.__binaries
        self.__is_feasible  = True
        self.__minAddlNodes = sol_num
        self.__timed_out = False
        self.__set_trivial_cols()
        assert(len(self.__binaries) >= 1)
        width = len(self.__binaries[0])
        non_trivial_width = width - len(self.__trivial_cols)
        known_h = set([hash_as_str([x[i] for i in range(width) if i not in self.__trivial_cols]) for x in self.__binaries])
        
        keep = []
        counter = 0
        for r in produce_random([0,1], non_trivial_width):
            h = hash_as_str(r)
            if h in known_h:
                continue
            assert len(r) == non_trivial_width
            keep.append(r)
            counter += 1
            if counter >= sol_num:
                break
        if self.formulation == 'Random1' or self.formulation == 'RandomX':
            assert len(keep) == sol_num
        else:
            # RandomAll may output nothing
            assert len(keep) == sol_num or len(keep) == 0
        
        if len(keep) == 0:
            self.__is_feasible  = False
            self.__minAddlNodes = -1
            # print('nothing to add')
        else:
            self.__is_feasible  = True
            self.__minAddlNodes = sol_num
            # put back random positions 
            # print('trivial col', self.__trivial_cols)
            keepwtrivial = [None] * len(keep)
            for ikeep in range(len(keep)):
                l = keep[ikeep]
                # print('nontrivial random', l)
                p = [None] * width # real vector with trivial pos
                for k, v in self.__trivial_cols.items(): 
                    p[k] = v

                l_idx = 0
                for i in range(width):
                    if p[i] is None:
                        assert l_idx < non_trivial_width
                        p[i] = l[l_idx]
                        l_idx += 1
                keepwtrivial[ikeep] = p
            for k in keepwtrivial:
                assert None not in k

            self.__novel_binaries = keepwtrivial
            self.__novel_info = [dict()] * len(keepwtrivial)  # make info and binary same length
            print(f'GeneChainToTree random {sol_num} sols, Add\'l Nodes {self.__minAddlNodes}')

    # self.__trivial_cols = dict() # col: val
    def __set_trivial_cols(self):
        assert(len(self.__trivial_cols) == 0)
        assert(len(self.__binaries) >= 1)
        width = len(self.__binaries[0])
        for i in range(width):
            s = sum([x[i] for x in self.__binaries]) 
            if s == 0:
                self.__trivial_cols[i] = 0
            elif s == len(self.__binaries):
                self.__trivial_cols[i] = 1
        return 0
        


class Transcriptom():
    def __init__(self, gtf: str, statsfile='', gtfpredfile='', ignore_single_exon_isoform: bool = True, args = None):
        self.gtf_input  = gtf
        self.statsfile   = statsfile   if statsfile   != '' else basename(gtf) + '.stats'
        self.gtfpredfile = gtfpredfile if gtfpredfile != '' else basename(gtf) + '.pred.gtf'
        if os.path.exists(self.statsfile):
            print(f"File {self.statsfile} exists. Now removing it.", file=sys.stderr)
            os.remove(self.statsfile)
        if os.path.exists(self.gtfpredfile):
            print(f"File {self.gtfpredfile} exists. Now removing it.", file=sys.stderr)
            os.remove(self.gtfpredfile)
        self.args = args
        self.ignore_single_exon_isoform = ignore_single_exon_isoform
        self.maxAddlNodes           = args.max_novel_isoform if args!= None else 4
        self.countSingleIsoformGene = 0
        self.geneCount              = 0
        self.bigGene                = 0
        self.infeasiCount           = 0
        self.geneAddlnodeCounts     = [0] * (self.maxAddlNodes + 1)
        self.genes                  = list()
        self.matrix_gids            = list()
        self.IsoCountAddlnodeCounts = defaultdict(int)  #d[isoCount, AddlNum] = # occur
        self.AddlnodeIsoCountCounts = defaultdict(int)  #d[AddlNum, isoCount] = # occur
        self.isoformsDict           = defaultdict(list)
        self.genesDict              = defaultdict(list)
        self.exonMatrixCollection   = list()
        self.gene2basicinfo         = dict()            #genes2basicinfo[gid] = [fields[x] for x in ['seqname', 'strand', 'frame']]
        self.parse_gtf(gtf)    
        
    def get_trees(self, chain_type: str = 'pexon_chain', transcript_group: str = 'tsstes_level', 
                  to_save: bool = True, statsfile='', gtfpredfile='',
                  formulation: str = 'SATSimple', xi_counts=None):
        if statsfile == '':
            statsfile = self.statsfile
        if gtfpredfile == '':
            gtfpredfile = self.gtfpredfile
        chains_of_all_genes = self.get_chain_matrix(chain_type, transcript_group)
        assert len(self.matrix_gids) == len(chains_of_all_genes)
        for i, chains_1_gene in enumerate(chains_of_all_genes):
            isoN = len(chains_1_gene)
            gid = self.matrix_gids[i]
            if isoN == 0:
                continue
            self.geneCount += 1
            if isoN <= 1:
                self.countSingleIsoformGene += 1
                continue
            if isoN >= self.args.exclude_group_size:
                self.bigGene += 1
                continue
            print(f'Processing a transcript group of size {isoN} in gene {gid}')
            
            if (formulation == 'Random1') or (formulation == 'RandomX'):                
                assert xi_counts is not None
                tsstes = f"s{chains_1_gene[0][0][1]}_t{chains_1_gene[0][-1][0]}"
                txGroup = gid + "." + tsstes 
                if formulation != 'RandomAll' and ((txGroup not in xi_counts) or (xi_counts[txGroup] < 1)): continue
                randOutNum = xi_counts[txGroup] if formulation == 'RandomX' else 1
                x = GeneChainToTree(chains_1_gene, randOutNum, formulation=formulation, args=self.args)
            else:
                x = GeneChainToTree(chains_1_gene, self.maxAddlNodes, formulation=formulation, args=self.args)
            
            self.infeasiCount += 1 if not x.is_feasible() else 0
            addlN = x.get_minAddlNodes()
            if x.is_feasible():
                self.geneAddlnodeCounts[addlN] += 1
                self.save_novel_gene_in_gtf(x.get_novel_transcripts(), x.get_tx_info(), gid, chain_type, addlN, file=gtfpredfile)
            self.IsoCountAddlnodeCounts[(addlN, isoN)] += 1
            self.AddlnodeIsoCountCounts[(isoN, addlN)] += 1 

            if to_save and self.geneCount % 10 == 0:
                self.save_stats(statsfile, finished=False)
        if to_save:
            self.save_stats(statsfile, finished=True)
        return 0
    
    def save_stats(self, file='', finished = False):
        if file == '':
            file = self.statsfile
        s = []
        if not finished:
            s.append(f'Processing file {self.gtf_input}')
        else:
            d = datetime.today().strftime('%Y-%m-%d')
            t = datetime.today().strftime('%H:%M:%S')
            s.append(f'Finished processing file {self.gtf_input} at {d} {t}')
        s.append(f'# Processed transcript groups:\t{self.geneCount}')
        s.append(f'# Single isoform group:\t{self.countSingleIsoformGene}')
        s.append(f'# Excl. group w. {self.args.exclude_group_size}+ isoforms:\t{self.bigGene}')
        for i, c in enumerate(self.geneAddlnodeCounts):
            # s.append(f'# feasible counts: {str(self.geneAddlnodeCounts)}')     
            s.append(f'# T{i} groups:\t{c}')
        s.append(f'# T@ groups:\t{self.infeasiCount}')
        # s.append(f'IsoCountAddlnodeCounts: {pformat(self.IsoCountAddlnodeCounts.items())}')
        # s.append(f'AddlnodeIsoCountCounts: {pformat(self.AddlnodeIsoCountCounts.items())}')
        s.append('\n')
        with open(file, 'w') as f:
            f.write('\n'.join(s))
        return 0

    def save_novel_gene_in_gtf(self, chains, chains_info, gid, chain_type, addlN, file = ''):
        if addlN <= 0:
            return 0
        
        file = self.gtfpredfile if file == '' else file            
        assert file.endswith('.gtf')
        
        if chain_type == "exon_chain" or chain_type == "pexon_chain":
            exon_lists = self.__exon_chain_matrix_to_transcripts(chains)
        elif chain_type == "intron_chain" or chain_type == "pintron_chain":
            exon_lists = self.__intron_chain_matrix_to_exon_list(chains)
        else:
            raise RuntimeError(f"{chain_type} chain type to transript not implemented")
        
        f = open(file, 'a')
        if len(chains_info) > 0:
            assert len(chains_info) == len(exon_lists)
        for i, exon_list in enumerate(exon_lists):
            tsstes = f"s{exon_list[0][1]}_t{exon_list[-1][0]}"
            tid = gid + "." + tsstes + ".novel."+ str(i + 1)
            # filter low PctIn isoforms
            if len(chains_info) > i and 'PctIn' in chains_info[i]:
                if chains_info[i]['PctIn'] < self.args.PctIn_threshold:
                    continue
            attr_info = f'novel_transcript_num \"{addlN}\"; ' 
            if len(chains_info) > 0:
                attr_info += ' '.join([f'{k} "{v}";' for k,v in chains_info[i].items()])
            f.write(self.__exon_list_to_gtf_lines(exon_list, gid, tid, attr_info))
        f.close()
        return 0

    def __exon_list_to_gtf_lines(self, exon_list, gid, tid, attr_info="") -> str:
        info    = self.gene2basicinfo[gid]
        seqname = info[0]
        strand  = info[1]
        frame   = '.'
        source  = 'tennis'
        
        tx_start = exon_list[0][0]
        tx_end   = exon_list[-1][1]
        attr = f"gene_id \"{gid}\"; transcript_id \"{tid}\"; {attr_info}"
        line     = '\t'.join([seqname, source, 'transcript', str(tx_start), str(tx_end), str(1000), strand, frame, attr]) + '\n'
        for ex in exon_list:
            start = ex[0]
            end   = ex[1]
            line += '\t'.join([seqname, source, 'exon', str(start), str(end), str(1000), strand, frame, attr]) + '\n'
        return line
    
    def __intron_chain_matrix_to_exon_list(self, chains):
        intron_chains = self.__condense_list_of_features(chains)
        exon_lists = []
        for ichain in intron_chains:
            exon_chain = list()
            tss = max(ichain[0][0] - 100, 1) # dummy tss
            exon_chain.append((tss, ichain[0][0]))
            for i in range(len(ichain) - 1):
                exon = (ichain[i][1] + 1, ichain[i + 1][0] - 1)
                assert exon[0] <= exon[1]
                exon_chain.append(exon)
            tes = max(ichain[-1][1], ichain[-1][1] + 100) # dummy tes
            exon_chain.append((ichain[-1][1], tes))
            exon_lists.append(exon_chain)
        return exon_lists

    # condense a list of adjacent exons (can also condense introns etc.)
    def __condense_list_of_features(self, chains):
        exon_lists = []
        for chain in chains:
            exon_list = []
            assert len(chain[0]) == 2
            l = chain[0][0]
            r = chain[0][1]
            for i in range(len(chain) - 1): 
                assert len(chain[i + 1]) == 2
                r = chain[i][1]
                l2 = chain[i + 1][0]
                r2 = chain[i + 1][1]
                assert r < r2
                if r == l2 or r + 1 == l2:
                    r = r2
                else:
                    exon_list.append((l,r))
                    l = l2
                    r = r2
            exon_list.append((l,r))
            exon_lists.append(exon_list)
        return exon_lists

    def __exon_chain_matrix_to_transcripts(self, chains):
        return self.__condense_list_of_features(chains)

    def __get_exon_chain_matrix_gene_level(self):
        self.matrix_gids = self.genes
        #make sure inclusive & non-overlapping
        for geneExonMatrix in self.exonMatrixCollection:
            for isoform in geneExonMatrix:
                exon_list_valid(isoform)
        return self.exonMatrixCollection

    def __get_exon_chain_matrix_tsstes_level(self, trim_ends=True):
        matrix = self.__get_exon_chain_matrix_gene_level()
        new_matrix = []
        self.matrix_gids = [] 
        assert len(self.genes) == len(matrix)
        for i, chains_1_gene in enumerate(matrix):
            tsstes_group = defaultdict(list)
            tsstes_trimmed_ends = dict()
            for j, chain in enumerate(chains_1_gene):
                s = chain[0][1]
                t = chain[-1][0]
                tsstes = (s,t)
                tsstes_group[tsstes].append(j)
                if trim_ends:
                    if tsstes not in tsstes_trimmed_ends:
                        tsstes_trimmed_ends[tsstes] = (chain[0][0], chain[-1][1])
                    else:
                        tsstes_trimmed_ends[tsstes] = (max(tsstes_trimmed_ends[tsstes][0], chain[0][0]), min(tsstes_trimmed_ends[tsstes][1], chain[-1][1]))
            for tsstes, chain_indices in tsstes_group.items():
                s, t = tsstes
                grouped_chains = [chains_1_gene[x] for x in chain_indices]
                if trim_ends:
                    for chain in grouped_chains:
                        chain[0] = (tsstes_trimmed_ends[tsstes][0], s) 
                        chain[-1] = (t, tsstes_trimmed_ends[tsstes][1])
                new_matrix.append(grouped_chains)
                self.matrix_gids.append(self.genes[i])
        assert len(self.matrix_gids) == len(new_matrix)
        return new_matrix

    def get_chain_matrix(self, chain_type: str = 'exon_chain', transcript_group: str = 'tsstes_level'):
        if transcript_group == 'gene_level':
            matrix = self.__get_exon_chain_matrix_gene_level()
        elif transcript_group == 'tsstes_level':
            matrix = self.__get_exon_chain_matrix_tsstes_level()
            # examine same tss/tes
            for group in matrix:
                assert (len(group) >= 1)
                s = group[0][0]
                t = group[0][-1]
                for isoform in group:
                    assert isoform[0] == s
                    assert isoform[-1] == t
        else:
            assert 0

        if chain_type == 'exon_chain':
            matrix = matrix
        elif chain_type == 'intron_chain':
            matrix = self.__get_intron_chain_matrix(matrix)
        elif chain_type == 'pexon_chain' or chain_type == 'partialexon_chain':
            matrix = self.__get_partialexon_chain_matrix(matrix)
        elif chain_type == 'pintron_chain' or chain_type == 'partialintron_chain':
            matrix = self.__get_partialintron_chain_matrix(matrix)
        else:
            assert 0

        for i, isoforms in enumerate(matrix): 
            matrix[i] = self.__remove_duplicated_isoforms(isoforms, self.ignore_single_exon_isoform)
        return matrix

    def __remove_duplicated_isoforms(self, isoforms, ignore_single_exon_isoform: bool = True):
        str_hash_set = set()    
        new_isoforms = []    
        for isoform in isoforms:
            h = hash_as_str(isoform)
            if ignore_single_exon_isoform and len(isoform) <= 1 :
                continue
            if h in str_hash_set:
                continue
            str_hash_set.add(h)
            new_isoforms.append(isoform)
        return new_isoforms

    def __get_partialintron_chain_matrix(self, matrix):
        m = self.__get_intron_chain_matrix(matrix)
        return self.__get_partial_any_chain_matrix(m)

    # input matrix - inclusive
    # return matrix - inclusive
    def __get_partial_any_chain_matrix(self, original_matrix):
        M = []
        for isoforms_in_gene in original_matrix:
            all_sites = set()
            for exons_in_isoform in isoforms_in_gene:
                for ex in exons_in_isoform:
                    all_sites.add(ex[0])  
                    all_sites.add(ex[1] + 1)  
            all_sites = sorted(list(all_sites))
            idx_sites = dict()
            for idx, site in enumerate(all_sites):
                idx_sites[site] = idx

            pexonchain_in_gene = []
            for exons_in_isoform in isoforms_in_gene:
                pexon_chain = []
                for ex in exons_in_isoform:
                    idx1 = idx_sites[ex[0]]
                    idx2 = idx_sites[ex[1] + 1]
                    assert idx1 <= idx2
                    pex_in_ex = []
                    for i in range(idx1, idx2):
                        pex_in_ex.append((all_sites[i], all_sites[i + 1] - 1))
                    if idx1 == idx2:
                        pex_in_ex.append((all_sites[idx1], all_sites[idx1] - 1))
                    assert len(pex_in_ex) >= 1
                    assert pex_in_ex[0][0] == ex[0]

                    for i, pex in enumerate(pex_in_ex):
                        if i == 0:
                            assert pex_in_ex[0][0] == ex[0]
                        elif i == len(pex_in_ex) - 1:
                            assert pex_in_ex[-1][1] == ex[1]
                            break
                        if i < len(pex_in_ex) - 1:
                            assert pex_in_ex[i][1] < pex_in_ex[i + 1][0] # inclusive

                    pexon_chain.extend(pex_in_ex)
                pexonchain_in_gene.append(pexon_chain)
            M.append(pexonchain_in_gene)
        return M

    def __get_partialexon_chain_matrix(self, matrix):
        return self.__get_partial_any_chain_matrix(matrix)

    def __get_intron_chain_matrix(self, matrix):
        M = []
        for isoforms_in_gene in matrix:
            intronchain_in_gene = []
            for exons_in_isoform in isoforms_in_gene:
                intron_chain = list()
                if len(exons_in_isoform) <= 1: 
                    intron_chain.append((0,0))
                else:
                    for i in range(len(exons_in_isoform) - 1):
                        intron = (exons_in_isoform[i][1] + 1, exons_in_isoform[i + 1][0] - 1)
                        assert intron[0] <= intron[1]
                        intron_chain.append(intron)
                intronchain_in_gene.append(intron_chain)
            M.append(intronchain_in_gene)
        return M

    def parse_gtf(self, gtf):
        picklefile = ""
        if exists(f'{gtf}.pickle'):
            picklefile = f'{gtf}.pickle'
        if gtf.endswith('pickle'):
            picklefile = gtf
        if picklefile != "":
            picklehandle = open(picklefile, 'rb')
            isoform2exons, gene2isoforms, exonMatrixCollection, genes, gene2basicinfo = pickle.load(picklehandle)
            self.genes                = genes
            self.isoformsDict         = isoform2exons
            self.genesDict            = gene2isoforms
            self.exonMatrixCollection = exonMatrixCollection
            self.gene2basicinfo       = gene2basicinfo
            print(f"GTF read from pickle file {picklefile}")
            picklehandle.close()
            return 0
        
        isoform2exons = defaultdict(list)
        gene2isoforms = defaultdict(list)
        gene2basicinfo = dict()
        isoform2gene  = dict()
        genes = list()
        table = list()
        with open(gtf, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'): continue
                fields = parse_gtf_line(line)
                table.append(fields)
                
                if fields['feature'] == 'transcript':
                    tid = fields['transcript_id']
                    gid = fields['gene_id']
                    assert tid is not None
                    assert gid is not None
                    gene2isoforms[gid].append(tid)
                    if gid not in gene2basicinfo:
                        gene2basicinfo[gid] = [fields[x] for x in ['seqname', 'strand', 'frame']]
                        genes.append(gid)
                if fields['feature'] == 'exon':
                    tid = fields['transcript_id']
                    gid = fields['gene_id']
                    s = int(fields['start'])
                    t = int(fields['end'])
                    assert tid is not None
                    assert gid is not None
                    assert s is not None
                    assert t is not None
                    isoform2exons[tid].append((s,t))
                    if tid in isoform2gene:
                        assert isoform2gene[tid] == gid
                    else:
                        isoform2gene[tid] = gid
        for _, exs in isoform2exons.items():
            exs.sort()

        exonMatrixCollection = []
        assert (len(gene2isoforms) == len(genes))
        for g in genes:
            is_in_g = gene2isoforms[g]
            exonMatrix = []
            for i in is_in_g:
                exonList = isoform2exons[i]
                exon_list_valid(exonList)
                exonMatrix.append(exonList)
            exonMatrixCollection.append(exonMatrix)
        assert (len(exonMatrixCollection) == len(genes))
        
        with open(f'{gtf}.pickle', 'wb') as file:
            pickle.dump((isoform2exons, gene2isoforms, exonMatrixCollection, genes, gene2basicinfo), file)
        self.genes                = genes
        self.isoformsDict         = isoform2exons
        self.genesDict            = gene2isoforms
        self.exonMatrixCollection = exonMatrixCollection
        self.gene2basicinfo       = gene2basicinfo
        print(f"GTF read from file {gtf}")
        return 0

def parse_arguments():    
    # if argv[2] is "test", parse only argv[2:] and set test options
    is_test = True if (len(sys.argv) > 1 and sys.argv[1] == "test") else False
    args    = sys.argv[1:] if not is_test else sys.argv[2:]

    parser = argparse.ArgumentParser(description="TENNIS -- Transcript EvolutioN for New Isoform Splicing")
    
    parser.add_argument("-o", "--output_prefix",      type=str, default="tennis", help="Output prefix for stats and prediction files")
    parser.add_argument("-p", "--PctIn_threshold",    type=float, default=0.5,    help="Predicted isoforms with PctIn value lower than this threshold will be filtered out.")
    parser.add_argument("-x", "--exclude_group_size", type=int, default=100,      help="Exclude large transcript group with size greater than this value")
    parser.add_argument("-m", "--max_novel_isoform",  type=int, default=4,        help="Maximum number of allowed novel isoforms in a transcript group")
    parser.add_argument(      "--time_out",           type=int, default=900,      help="Each SAT instance time out in seconds")
    if is_test:
        parser.add_argument("-f", "--formulation", type=str, default="SATSimple", choices=["SATSimple", "Random1", "RandomX"], help="Formulation type")
        parser.add_argument("--xi_gtf_file", type=str, default=None, help="gtf file from TENNIS with Ti information")
    else:
        parser.formulation = "SATSimple"

    parser.add_argument("gtf_file", type=str, help="Input GTF file")
    args = parser.parse_args(args)
    if not is_test:
        args.formulation = "SATSimple"
    
    if args.PctIn_threshold > 1 or args.PctIn_threshold < 0:
        raise ValueError("--PctIn_threshold must be between 0 and 1")
    if args.max_novel_isoform > 10:
        print(f"Warning: --max_novel_isoform set to {args.max_novel_isoform}. Default is 4. Allowing too many isoforms will slow down the program.", file=sys.stderr)
        sleep(10)
    return args



def main():
    random.seed(2024)
    d0 = datetime.today().strftime('%Y-%m-%d')
    t0 = datetime.today().strftime('%H:%M:%S')
    
    args = parse_arguments()

    gtf_file         = args.gtf_file
    formulation      = args.formulation
    chain_type       = 'pexon_chain'
    transcript_group = 'tsstes_level' 
    save_basename = args.output_prefix 

    tsm = Transcriptom(gtf_file, f'{save_basename}.stats', f'{save_basename}.pred.gtf', args=args)

    if (formulation=='Random1') or (formulation == 'RandomX'):
        assert args.xi_gtf_file is not None 
        _, xic = get_xi_counts(args.xi_gtf_file)
        tsm.get_trees(chain_type, transcript_group, statsfile=f'{save_basename}.stats', 
                      gtfpredfile=f'{save_basename}.pred.gtf', formulation=formulation, xi_counts=xic)
    else:
        tsm.get_trees(statsfile=f'{save_basename}.stats', gtfpredfile=f'{save_basename}.pred.gtf', formulation=formulation)
    
    d = datetime.today().strftime('%Y-%m-%d')
    t = datetime.today().strftime('%H:%M:%S')
    print(f"TENNIS completed tunning! \n" + 
          f"Started  at time {d0} {t0} \n" +
          f"Finished at time {d} {t} \n" +
          f"Output files:\n" +
          " " * 8 + f"{save_basename}.stats\n"+
          " " * 8 + f"{save_basename}.pred.gtf")


if __name__ == "__main__":
    main()