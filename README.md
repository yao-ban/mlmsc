# The MLMSC-II Simulator

The MLMSC-II Simulator is a program for the simulation of gene family evolution within a species tree based on the updated Multilocus Multipecies Coalescent (MLMSC-II) model. The MLMSC-II model generalises the multispecies coalescent to gene families, and is designed to capture all possible scenarios that can arise through incomplete lineage sorting, gene duplication, transfer and loss, and any interaction between these processes. The MLMSC-II combines forward- and backward-in-time modelling in order to properly account for copy number hemiplasy and linkage between loci. 

## Citation

Pending...

The simulation scripts and dataset for the paper **The effect of copy number hemiplasy on gene family evolution** can be found at https://github.com/QiuyiLi/MLMSC-II_simulation_script

## Obtaining the program

### Download

The source code of MLMSC Simulator is available in the GitHub repository:
```
https://github.com/QiuyiLi/MLMSC-II
```
You can clone the sources in your computer by executing:
```
git clone https://github.com/QiuyiLi/MLMSC-II.git
```

### Requirements
MLMSC Simulator is built under Python 3.7.2 and requirs the installation of the following packages: 
```
scikit-bio
numpy
statistics
```
You can install these packages by executing:
```
pip install -r requirements.txt
```
### Testing

The MLMSC Simulator should be ready to use now. You can test the the simulator by executing:
```
python3 MLMSC.py -i species_trees/ex_tree.newick -s 0
```
The simulator will run with the default settings and you should be able to see the following output:
```
gene tree:
                                                  /-A
                                        /--------|
                              /--------|          \-F
                             |         |
                             |          \-H
                    /--------|
                   |         |                    /-C
                   |         |          /--------|
          /--------|          \--------|          \-D
         |         |                   |
         |         |                    \-G
         |         |
         |          \-E
---------|
         |                    /-K
         |          /--------|
         |         |         |          /-J
         |         |          \--------|
          \--------|                    \-O
                   |
                   |          /-L
                    \--------|
                              \-P
finished.
```

##  Usage

### Input file
* Species tree: -i or --inputFile; the prespecified species tree written in [Newick](https://en.wikipedia.org/wiki/Newick_format) format
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick
```

### Input parameters

* Duplication rate: -d or --duplicationRate; occurrence rate of gene duplication, default: -d 0.2
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -d 0.1
```
  
* Transfer rate; -t or --transferRate; occurrence rate of horizontal gene transfer, default: -t 0.1
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -t 0.2
```
  
* Loss rate; -l or --lossRate; occurrence rate of gene loss, default: -l 0.2
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -l 0.1
```

* Coalescent rate: -c or --coalescentRate; the coalescent rate in multispecies coalescent, default: -c 1
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -c 0.5
```
  
* Recombination rate: -r or --recombinationRate; the recombiantion rate in linked coalescent, default: -r 0
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -r 0.3
```
  
* Unlink probability; -u or --unlinkProb; the probability that a duplication is unlinked, default: -u 1
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -u 0.3
```
  
* Verbose option: -v or --verbose; detailed outputs for debugging purposes, default: -v 0
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -v 1
```

* Number of repeats: -n or --numRepeats; number of gene trees simulated, defult: -n 1
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -n 10
```

* Seed number: -s or --seed; set seed for reproduciblility, default -s None (random seed)
```
e.g., python3 MLMSC.py -i species_trees/ex_tree.newick -s 0
```

### Output files
* gene_tree.newick: a multi-labelled gene tree consisting of all homologous genes
* random_tree.newick: a single-labelled gene tree in which only one gene is randomly selected for each descendant species 
* summary.txt: summary statistics including the number of surviving duplications (n_d), number of genes (n_genes), and number of species (n_species)

## Authors

* **Qiuyi Li** 
* **[Geoffrey Law](https://github.com/luojiahai)**

## License

This project is licensed under the GNU GPLv3 - see [LICENSE.md](LICENSE.md) for more details.


