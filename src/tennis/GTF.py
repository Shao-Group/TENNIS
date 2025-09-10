#!/usr/bin/env python
"""
GTF.py

Originally developed by Kamil Slowikowski (https://gist.github.com/slowkow/8101481?permalink_comment_id=321645i7)
Modified and re-distributed by Xiaofei Carl Zang
"""

"""
Original License by Kamil Slowikowski

Kamil Slowikowski
December 24, 2013

Read GFF/GTF files. Works with gzip compressed files and pandas.

    http://useast.ensembl.org/info/website/upload/gff.html

LICENSE

This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""


from collections import defaultdict
import gzip
import re


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    #return pd.DataFrame(result)
    return result


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)
    return result

def write_file(lines, dest):
    f = open(dest, 'w')
    for line in lines:
        sorted_key = ['gene_id', 'transcript_id'] + sorted([x for x in line.keys() if x not in ['gene_id', 'transcript_id'] + GTF_HEADER])
        attributes_column = ' '.join([f'{k} "{line[k]}";' for k in sorted_key])
        eight_columns = [line[col] if line[col] != None else '.' for col in GTF_HEADER]
        line_str = '\t'.join(eight_columns) + '\t' + attributes_column
        f.write(line_str + '\n')    
    f.close()
    return 0


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def get_xi_counts(gtf):
    gid2xi = dict()
    tgroup2xi = dict()
    f = open(gtf, 'r')

    for line in f.readlines():
        if line.startswith('#'): continue
        fields = parse(line)
        
        tgroup = fields['transcript_id'].split('.novel')[0]
        gid = fields['gene_id']
        xi = int(fields['novel_transcript_num'])

        if gid not in gid2xi:
            gid2xi[gid] = xi      
        if tgroup not in tgroup2xi:
            tgroup2xi[tgroup] = xi
        else:
            assert tgroup2xi[tgroup] == xi
    f.close()
    return gid2xi, tgroup2xi

def get_PctIn(gtf):
    tid2PctIn = dict()

    f = open(gtf, 'r')
    for line in f.readlines():
        if line.startswith('#'): continue
        fields = parse(line)
        gid = fields['gene_id']
        tid = fields['transcript_id']
        pctIn = float(fields['PctIn'])

        if tid not in tid2PctIn:
            tid2PctIn[tid] = pctIn
        else:
            assert tid2PctIn[tid] == pctIn
    f.close()
    return tid2PctIn