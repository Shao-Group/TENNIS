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
import math
from .GTF import parse as parse_gtf_line

def parse_psi_file(psi_file, chr_translate=None):
    """
    Parse PSI file and return a dictionary mapping (chrom, start, end) -> average PSI value.

    PSI file format (tab-separated):
    - Column 0: GENE
    - Column 1: EVENT
    - Column 2: COORD (e.g., chr3R:16935927-16936104)
    - Column 3: LENGTH
    - Column 4: FullCO
    - Column 5: COMPLEX
    - Column 6, 8, 10, ...: PSI values (can be float or "NA")
    - Column 7, 9, 11, ...: Quality values (ignored)

    Args:
        psi_file: Path to PSI file
        chr_translate: Optional dict {psi_chr_name: gtf_chr_name} to translate PSI chr names to GTF chr names

    Returns:
        dict: {(chrom, start, end): average_psi} where chrom is in GTF naming convention
    """
    coord_to_psi = dict()

    with open(psi_file, 'r') as f:
        header = f.readline()  # Skip header line
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 7:
                continue

            # Parse coordinate: chr3R:16935927-16936104
            coord_str = fields[2]
            try:
                chrom, pos = coord_str.split(':')
                start, end = pos.split('-')
                start = int(start)
                end = int(end)
            except (ValueError, IndexError):
                continue

            # Translate PSI chromosome name to GTF chromosome name if translation provided
            if chr_translate is not None and chrom in chr_translate:
                chrom = chr_translate[chrom]

            # Collect PSI values from columns 6, 8, 10, ... (0-indexed)
            psi_values = []
            for i in range(6, len(fields), 2):
                psi_str = fields[i].strip()
                if psi_str == 'NA' or psi_str == '':
                    continue
                try:
                    psi_val = float(psi_str)/100.0 # percentage to fraction
                    psi_values.append(psi_val)
                except ValueError:
                    continue

            # Calculate average PSI
            if psi_values:
                avg_psi = sum(psi_values) / len(psi_values)
                coord_to_psi[(chrom, start, end)] = avg_psi

    return coord_to_psi


def parse_chr_translate_file(translate_file):
    """
    Parse chromosome name translation file.

    File format (tab-separated, no header):
    - Column 0: chr name in PSI file (to translate)
    - Column 1: chr name in GTF file (to keep)

    Returns:
        dict: {psi_chr_name: gtf_chr_name} (translate PSI to GTF)
    """
    psi_to_gtf = dict()

    with open(translate_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 2:
                continue
            psi_chr = fields[0].strip()
            gtf_chr = fields[1].strip()
            psi_to_gtf[psi_chr] = gtf_chr

    return psi_to_gtf


def compute_psi_scores_for_gtf(input_gtf, original_gtf, psi_data, output_gtf):
    """
    Compute PSI scores for predicted isoforms and write to output GTF.
    The PSI_score is the product of all exon probabilities.
    """

    # Step 1: Parse original GTF to get exon usage statistics per gene
    # gene_id -> {exon_coord -> count of isoforms containing it}
    # gene_id -> total isoform count
    gene_exon_counts = dict()  # gene_id -> {(start, end): count}
    gene_isoform_counts = dict()  # gene_id -> total isoforms
    gene_chrom = dict()  # gene_id -> chromosome

    # First pass: count isoforms per gene
    transcript_exons = dict()  # transcript_id -> [(start, end), ...]
    transcript_gene = dict()  # transcript_id -> gene_id

    with open(original_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = parse_gtf_line(line)
            if fields['feature'] == 'transcript':
                gid = fields['gene_id']
                tid = fields['transcript_id']
                transcript_gene[tid] = gid
                if gid not in gene_isoform_counts:
                    gene_isoform_counts[gid] = 0
                    gene_exon_counts[gid] = dict()
                    gene_chrom[gid] = fields['seqname']
                gene_isoform_counts[gid] += 1
                transcript_exons[tid] = []
            elif fields['feature'] == 'exon':
                tid = fields['transcript_id']
                start = int(fields['start'])
                end = int(fields['end'])
                if tid in transcript_exons:
                    transcript_exons[tid].append((start, end))

    # Count exon usage per gene
    for tid, exons in transcript_exons.items():
        gid = transcript_gene.get(tid)
        if gid is None:
            continue
        for exon in exons:
            if exon not in gene_exon_counts[gid]:
                gene_exon_counts[gid][exon] = 0
            gene_exon_counts[gid][exon] += 1

    # Step 2: Parse predicted GTF and compute PSI scores
    predicted_transcripts = dict()  # transcript_id -> {'exons': [...], 'line_fields': [...]}
    predicted_lines = []  # All lines from predicted GTF

    with open(input_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                predicted_lines.append(line)
                continue
            fields = parse_gtf_line(line)
            predicted_lines.append(fields)

            if fields['feature'] == 'exon':
                tid = fields['transcript_id']
                gid = fields['gene_id']
                start = int(fields['start'])
                end = int(fields['end'])
                chrom = fields['seqname']

                if tid not in predicted_transcripts:
                    predicted_transcripts[tid] = {
                        'exons': [],
                        'gene_id': gid,
                        'chrom': chrom
                    }
                predicted_transcripts[tid]['exons'].append((start, end))

    # Step 3: Compute PSI score for each predicted transcript
    transcript_psi_scores = dict()

    for tid, tdata in predicted_transcripts.items():
        exons = tdata['exons']
        gid = tdata['gene_id']
        chrom = tdata['chrom']

        log_prob = 0.0
        for exon in exons:
            psi = _get_exon_psi(
                exon, chrom, gid,
                psi_data, gene_exon_counts, gene_isoform_counts
            )
            # PSI is probability of inclusion (0-100), convert to 0-1
            prob = max(psi / 100.0, 1e-10)
            log_prob += math.log(prob)

        # Convert log probability to a score (exp of log_prob, or keep as log)
        # Using exp can result in very small numbers, so we'll store the log probability
        # Alternatively, compute geometric mean or just the raw probability
        psi_score = math.exp(log_prob)
        transcript_psi_scores[tid] = psi_score

    # Step 4: Write output GTF with PSI_score attribute
    with open(output_gtf, 'w') as f:
        for line in predicted_lines:
            if isinstance(line, str):
                # Comment line
                f.write(line)
            else:
                # Add PSI_score to transcript and exon lines
                tid = line.get('transcript_id')
                if tid and tid in transcript_psi_scores:
                    line['PSI_score'] = f"{transcript_psi_scores[tid]:.6e}"

                # Write line in GTF format
                from .GTF import GTF_HEADER
                sorted_keys = ['gene_id', 'transcript_id'] + sorted(
                    [k for k in line.keys() if k not in ['gene_id', 'transcript_id'] + GTF_HEADER]
                )
                attributes = ' '.join([f'{k} "{line[k]}";' for k in sorted_keys if line.get(k) is not None])
                eight_cols = [str(line[col]) if line.get(col) is not None else '.' for col in GTF_HEADER]
                f.write('\t'.join(eight_cols) + '\t' + attributes + '\n')

    print(f"Wrote {len(transcript_psi_scores)} transcripts with PSI scores to {output_gtf}")
    return transcript_psi_scores


def _get_exon_psi(exon, chrom, gene_id, psi_data, gene_exon_counts, gene_isoform_counts):
    key = (chrom, exon[0], exon[1])
    if key in psi_data:
        return psi_data[key] * 100.0  # psi_data stores as fraction, convert to percentage

    if gene_id in gene_exon_counts:
        exon_count = gene_exon_counts[gene_id].get(exon, 0)
        total_isoforms = gene_isoform_counts.get(gene_id, 1)
        if exon_count > 0:
            return min((exon_count / total_isoforms) * 100.0, 100.0)
        else:
            return 50.0
    else:
        return 50.0
