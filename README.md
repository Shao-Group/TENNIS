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
```

This repository also modified and re-distributes GTF.py codes (retrieved from [here](https://gist.github.com/slowkow/8101481?permalink_comment_id=321645i7)) developed by Kamil Slowikowski. Users don't have to re-download it.

### Example

```sh
# display help message
python src/tennis.py -h
# run TENNIS on an example dataset
mkdir test
cd test
python ../src/tennis.py ../example/example.gtf tennis_example
```



# **Usage** 

```sh
python src/tennis.py <gtf_file> <output_prefix>
```

Two files will be output: `output_prefix.stats` and `output_prefix.pred.gtf`



# Contributing

For bug reports or feature requests, please open an issue on the GitHub repository.



# **License**

TENNIS is freely available under BSD 3-Clause License. 

Copyright (c) 2024, Xiaofei Carl Zang, Ke Chen, Mingfu Shao, and The Pennsylvania State University.
