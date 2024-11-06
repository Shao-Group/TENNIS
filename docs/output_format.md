# TENNIS outputs
TENNIS outputs two files: `output_prefix.stats` and `output_prefix.pred.gtf`.

## stats file

The file `output_prefix.stats` consists of several summary statistics:

- **Processed transcript groups**:  Total number of transcript groups analyzed. This number includes all multi-exon transcript groups but excludes all single-exon transcripts.

- **Single isoform group**: Number of groups containing only one isoform. Only multi-exon isoforms are counted.

- **Excl. group w. `x+` isoforms**: Number of groups excluded due to extra large size, i.e. having over `x+` isoforms. This parameter is set by `-x` or `--exclude_group_size`.

- **T0 groups**: Number of groups that satisfy the model without additional isoforms

- **T1, T2, ..., Tm groups**: Groups requiring 1, 2, ..., m novel isoforms to satisfy the model. For example, T2 groups require 2 novel isoforms to meet the model requirements. m is set by `-m` or `--max_novel_isoform`.

- **T@ groups**: Groups that require at least `m+1` isoforms to satisfy the model; or groups that cannot be computed within the time limit set by `--time_out`.


## pred.gtf file
The `output_prefix.pred.gtf` file has a standard pdf format with additional information in the attribute column (9th column). The attributes are:

| Tag | Value |
|---|---|
| gene_id | original gene id (same as input gtf) |
| transcript_id | transcript id |
| novel_transcript_num | Number of novel isoforms required by the transcript group to satisfy the evolution model. This attr reflects the T1, T2, ..., Tm information. |
| solutions_in | This isoform is in those solutions. For example [1,2,3] means this isoform appears in solution 1, solution2, and solution 3|
| solutions_total_num | Number of optimal solutions computed. Note this entry is different from "novel_transcript_num" |
| solution_unique | Is the solution unique or not. Namely, `solutions_total_num == 1`? |
| mandatory | This isoform appears in all solutions or not. Namely, `len(solutions_in) == solutions_total_num`? |
| PctIn | Ratio of solutions that this isoform appears in. Namely, `len(solutions_in) / solutions_total_num`. Apparently, PctIn = 1 iff mandatory |
| timed_out | If True, solver did not complete within the given time limit |

An example of attr col:
> 	gene_id "Dmel_CG3151"; transcript_id "Dmel_CG3151.s2959520_t2963225.novel.4"; novel_transcript_num "2"; solutions_in "[3, 4, 5]"; solutions_total_num "6"; solution_unique "False"; mandatory "False"; PctIn "0.5"; timed_out "False";