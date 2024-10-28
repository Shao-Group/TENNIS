# TENNIS 🎾: Transcript EvolutioN for New Isoform Splicing

TENNIS is a tool for the construction of evolutionary trajectories of transcripts with a minimal number of internal (missing) nodes.

## **Dependency**
TENNIS is implemented in Python (version >=3.7). TENNIS uses a SAT interface by [PySAT](https://pysathq.github.io/). PySAT can be installed via 

```sh
pip install python-sat[aiger,approxmc,cryptosat,pblib]
```


This repository also modified and re-distributes GTF.py codes (retrieved from [here](https://gist.github.com/slowkow/8101481?permalink_comment_id=321645i7)) developed by Kamil Slowikowski. Users don't have to re-download it.

# **Installation**

After installing the needed dependencies, TENNIS can be installed by directly cloning this repository.

```sh
git clone https://github.com/Shao-Group/tennis.git
```

# **Usage** 

```sh
python src/tennis.py <gtf_file> <output_prefix>
```

Two files will be output: output_prefix.stats and output_prefix.pred.gtf

# **License**

TENNIS is freely available under the BSD 3-Clause License.
