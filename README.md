# PhyloGTP

## **Description**
PhyloGTP is a software tool for genome-scale microbial phylogenomics. It takes as input a collection of gene trees and estimates a species tree under a duplication-transfer-loss (DTL) model of gene family evolution. PhyloGTP uses local search heuristics to obtain a species tree topology which minimizes the global DTL reconciliation cost with the collection of input gene trees. The program can easily be run in the command line or using a python interpreter. PhyloGTP can be cited as follows:

<a href="https://compbio-engr.media.uconn.edu/wp-content/uploads/sites/2447/2024/03/MicrobialPhylogenomics_PhyloGTP_RECOMBCG_2024.pdf">Assessing the potential of gene tree parsimony for microbial phylogenomics</a><br>
Samson Weiner, Yutian Feng, J. Peter Gogartan, and Mukul S. Bansal<br>
<i>RECOMB Comparative Genomics conference (RECOMB-CG) 2024</i>; to appear.

## Installation
To download PhyloGTP, you can clone the repository with
```
git clone https://github.com/samsonweiner/PhyloGTP.git
```

### Dependencies
PhyloGTP requires two common python packages are installed in the environment:
* [Numpy](https://numpy.org/)
* [Ete3](http://etetoolkit.org/)

These packages can be easily installed using your package manager of choice. An enviroment setup and installation using `conda` following best practices are shown below.
```
conda create -n PhyloGTP python=3
conda activate PhyloGTP

conda install numpy
conda install conda-forge::ete3
```

Additionally, PhyloGTP uses the [Ranger-DTL](https://compbio.engr.uconn.edu/software/ranger-dtl/) software package. For your convenience, executables for Linux, Windows, and Mac machines are already included in the repository in the binaries folder. PhyloGTP will automatically choose the platform-appropriate executable, however you must ensure the executables are given correct permissions. To do so, nagivate into the binaries folder and change permissions as follows:
```
cd binaries
chmod u+x *
```

If you experience an issue with the provided executables, you can download and install Ranger-DTL at the link provided. Note that only the Ranger-DTL-Fast executable under SupplementaryPrograms is required to run PhyloGTP.

## Input file format
The one required input file is the file containing the gene trees. Each gene tree must appear on a separate line in the file and should be expressed in Newick format terminated by a semicolon. For each gene tree, every leaf node must be labeled with the name of the species from which that gene was sampled. If desired, the gene name can be appended to the species name separated by an underscore ‘_’ character. The gene tree may contain any number (zero, one, or more) of homologous genes from the same species.

E.g.

(((speciesA_gene1, speciesC_gene1), speciesB_geneX), speciesC_gene2); <br>
(((speciesB, speciesC), speciesB), speciesA); <br>
((speciesD_gene1, speciesC), (speciesA_gene5, speciesD)); <br>

A weight may also be assigned to each gene tree, corresponding to the contribution of each gene tree towards the global reconciliation cost. To assign weights, prepend each gene tree with ‘[w]’ where `w` is the weight.

E.g.

[1.6](((speciesB, speciesC), speciesB), speciesA); <br>
[0.7]((speciesD_gene1, speciesC), (speciesA_gene5, speciesD)); <br>

Users may also provide a starting species tree used to seed the tree search. The starting species tree should be present in its own file and must also be expressed in Newick format.


## Example run
To run PhyloGTP, cd to `main` and execute the following command:
```
python phylogtp.py -i [path-to-input-genetrees.nwk] -o [path-to-output-directory]
```

We have provided several test datasets in the `testdata` directory which you may use to validate your installation of PhyloGTP. Here is an example of running one of the test datasets:
```
python phylogtp.py -i ../testdata/genetrees1.nwk -o test/ -u
```

## Available command line options
`-i, --input` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to input file containing gene trees and optional weights. (Required)

`-o, --output` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to output directory. (Default: out) 

`-s, --start-tree-file` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to species tree file used as the starting point in the search. If unspecified, will use an additive taxon procedure. (Default: None)

`-x, --random-start` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If True, will use start tree search from a random topology. Otherwise, will use an additive taxon procedure. (Default: False)

`-u, --unrooted` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Toggle to treat gene trees as unrooted. Will do so by default if the input gene trees have no root. (Default: False) 

`-r, --ranger-path` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to Ranger-DTL-Fast executable. If unspecified, will use the executable located in the binaries folder based on the user platform automatically. (Default: None). 

`-D, --dup-cost` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cost of duplication events when computing parsimonious DTL reconciliations. (Default: 2) 

`-T, --tran-cost` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cost of transfer events when computing parsimonious DTL reconciliations. (Default: 3) 

`-L, --loss-cost` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cost of loss events when computing parsimonious DTL reconciliations. (Default: 1) 

`-t, --threshold` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Terminates loop after difference in global reconcilaition cost falls below the given threshold. (Default: 0) 

`-m, --max-iter` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Max number of iterations before termination. If set to -1, runs until convergence. (Default: -1) 

`-f, --full-SPR` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use a full SPR search search instead of a greedy search. (Default: False)

`-c, --num-cores` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of CPUs to use in parallel. (Default: 1)

## Miscelaneous
We have provided scripts used to generate the simulated datasets used in our paper. These scripts are located in the `simulations` folder. In order to use these scripts, two external software packages are required: [SaGePhy](https://compbio.engr.uconn.edu/software/sagephy/) and [IQ-Tree](http://www.iqtree.org/) (including the [AliSim](http://www.iqtree.org/doc/AliSim) package).


## Contact
If there are any questions related to PhyloGTP, please contact Samson Weiner (<samson.weiner@uconn.edu>) or Mukul Bansal (<mukul.bansal@uconn.edu>), or visit <https://compbio.engr.uconn.edu/> for a complete list of available software and datasets.
