# TENNIS 🎾: Transcript EvolutioN for New Isoform Splicing

TENNIS is an evolution-based model to predict unannotated isoforms and refine existing transcriptome annotations without requiring additional data. 



# Installation

### Prerequisites

- Python >= 3.7
- [PySAT](https://pysathq.github.io/)

### Installation 

The only dependency of TENNIS is [PySAT](https://pysathq.github.io/), which can be installed with `pip`. TENNIS can be installed by directly cloning this repository.

```sh
# install PySAT
pip install python-sat[aiger,approxmc,cryptosat,pblib]
# install TENNIS
git clone https://github.com/Shao-Group/TENNIS
cd TENNIS
chmod +x src/tennis
```

This repository also modified and re-distributes GTF.py codes (retrieved from [here](https://gist.github.com/slowkow/8101481?permalink_comment_id=321645i7)) developed by Kamil Slowikowski. Users don't have to re-download it.

### Test and Example

```sh
# display help message
./src/tennis -h
# run TENNIS on an example dataset
mkdir test
cd test
./src/tennis -o tennis_example ../example/example.gtf 
```

# Usage
If installed with conda or pip, `tennis` executable should be ready to use in `$PATH`.
If installed manually, the `tennis` executable is in `./src/` dir. 
```sh
tennis [options] -o <output_prefix> <gtf_file> 
```

The program outputs two files: `output_prefix.stats` and `output_prefix.pred.gtf`. 
More about the output format is available [here](docs/output_format.md).

### Required positional arguments:

`gtf_file` : str  
Input GTF file in standard format containing transcript annotations.

### Optional arguments:

`-h`, `--help`

`-o`, `--output_prefix` : str  
Default: "tennis"

`-p`, `--PctIn_threshold` : float    
A threshold in range [0, 1]. Predicted isoforms with *PctIn* value lower than this threshold will be filtered out. If `-p 0.0`, all isoforms are retained.    
Default: 0.5

`-x`, `--exclude_group_size` : int  
Skip analysis of transcript groups that have more isoforms than this threshold.  
Default: 100

`-m`, `--max_novel_isoform` : int  
Maximum number of novel isoforms to predict per transcript group.  
Default: 4

`--time_out` : int   
Time limit in seconds for each SAT solver instance.  
Default: 900 (15 minutes)

## Output Files

`output_prefix.stats` : 
	Statistical summary. T1 (T2, T3, ...) is the collection of transcript groups that need 1 (2, 3, ..) novel isoforms to satisfy the evolution model.

`output_prefix.pred.gtf` : 
​	GTF format file with predicted novel isoforms.

More about the output format is available [here](docs/output_format.md).

# Contributing

For bug reports or feature requests, please open an issue on the GitHub repository [here](https://github.com/Shao-Group/TENNIS/issues).



# License & Citation

TENNIS is freely available under BSD 3-Clause License. 

Copyright (c) 2024, Xiaofei Carl Zang, Ke Chen, Mingfu Shao, and The Pennsylvania State University.

The preprint of TENNIS is available on bioRxiv [here](https://doi.org/10.1101/2024.11.04.621892).

```
@article {TENNIS,
	author = {Zang, Xiaofei Carl and Chen, Ke and Khan, Irtesam Mahmud and Shao, Mingfu},
	title = {Augmenting Transcriptome Annotations through the Lens of Splicing Evolution},
	year = {2024},
	doi = {10.1101/2024.11.04.621892},
	journal = {bioRxiv}
}
```

