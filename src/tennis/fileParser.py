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
