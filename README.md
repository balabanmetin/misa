

------------------------------------
Summary
------------------------------------
MISA stands for `MIxed Sample Analysis tool` and addresses the problem of phylogenetic double placement of mix (of two) DNA sequences into an already existing reference tree.  MISA is a command-line tool and it can run on **Linux, Mac OSX, and Windows**.


------------------------------------
Publication
------------------------------------

* Balaban, Metin, and Siavash Mirarab. “Phylogenetic Double Placement of Mixed Samples.” Bioinformatics, in press, (2020).


------------------------------------
Requirements
------------------------------------
1. Python: Version >= 3.0

------------------------------------
Installation on Linux, Mac OSX, or Windows
------------------------------------

Clone the repository and change directory:

`git clone https://github.com/balabanmetin/misa.git && cd misa`


---------------------------------------------
Getting Started with MISA
---------------------------------------------

For listing all options, run the following command:

`run_misa.py -h`

---------------------------------------------
Input & Output Specification
---------------------------------------------

Input reference (backbone) tree must be in newick format. MISA can perform placements based on a distance table. Distance table between the query and reference sequences can be computed using [https://github.com/shahab-sarmashghi/Skmer](Skmer). These distances can be based on assemblies or unassembled bag of reads (genome skims).


#### Input a distance matrix
The format for distance matrix is a tab delimited csv file with column and row headers. Rows should represent query sequences and columns should represent reference sequences. You can find an example distance matrix for ten query sequences under [data/dist.mat](data/dist.mat). 
You can run MISA on the example distance matrix by running the following command:

`run_misa.py -d data/dist.mat -t data/backbone.nwk`

MISA has the following dependencies which can be retrieved from `pip`:

1. `treeswift`
2. `numpy`
3. `scipy`

To install a dependency, for example `scipy`, run `pip3 install scipy`.

#### Output
Output is a jplace file containing placement results for all queries. For more information about jplace files, please refer to Matsen et. al. (2012) [https://doi.org/10.1371/journal.pone.0031009](https://doi.org/10.1371/journal.pone.0031009). The output file can be specified using `-o` command. When output file is not specified, the result will be printed to the standard output.

### ! IMPORTANT NOTE !

Backbone tree provided to MISA has to have its branch lengths estimated using a distance based method such as minimum evolution. This is a requirement for getting good results. We recommend [FastTree2](http://www.microbesonline.org/fasttree/) for re-estimating branch lengths if the backbone tree is estimated using Maximum Likelihood based methods (e.g. RAxML, PASTA). Until we support re-estimation of branch lengths within MISA, we are expecting the user to run the following  FastTree2 command explicitly before performing any placement:

`FastTreeMP -nosupport -nt -nome -noml -log tree.log -intree backbone.nwk < ref.fa > minimum_evo_backbone.nwk`

Then perform placement on the new tree:

`run_misa.py -d dist.mat -t minimum_evo_backbone.nwk`
