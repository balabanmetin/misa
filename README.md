

------------------------------------
Summary
------------------------------------
MISA stands for `MIxed Sample Analysis tool` and addresses the problem of phylogenetic double placement of mix (of two) DNA sequences into an already existing reference tree.  MISA is a command-line tool and it can run on **Linux, Mac OSX, and Windows**.


------------------------------------
Publication
------------------------------------

* Balaban, Metin, and Siavash Mirarab. “Phylogenetic double placement of mixed samples.” Bioinformatics (Oxford, England) vol. 36,Supplement_1 (2020): i335-i343. [https://doi.org/10.1093/bioinformatics/btaa489](https://doi.org/10.1093/bioinformatics/btaa489)


------------------------------------
Requirements
------------------------------------
1. Python: Version >= 3.0

------------------------------------
Installation on Linux, Mac OSX, or Windows
------------------------------------

Install MISA using the following command in the command-line:

`pip install misa`


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
The format for distance matrix is a tab delimited csv file with column and row headers. Rows should represent query sequences and columns should represent reference sequences. You can find an example reference dataset with backbone phylogenetic tree under [data/backbone.nwk](data/backbone.nwk) and distance matrix for one query mixed sequence under [data/dist.mat](data/dist.mat).
After cloning this repository, you can run MISA on the example distance matrix and backbone tree by running the following command:

`run_misa.py -d <path_to_this_repo>/data/dist.mat -t <path_to_this_repo>/data/backbone.nwk`


#### Output
Output is a jplace file containing placement results for all queries. For more information about jplace files, please refer to Matsen et. al. (2012) [https://doi.org/10.1371/journal.pone.0031009](https://doi.org/10.1371/journal.pone.0031009). The output file can be specified using `-o` command. When output file is not specified, the result will be printed to the standard output.

### ! IMPORTANT NOTE !

Backbone tree provided to MISA has to have its branch lengths estimated using a distance based method such as minimum evolution. This is a requirement for getting good results. We recommend [FastTree2](http://www.microbesonline.org/fasttree/) for re-estimating branch lengths if the backbone tree is estimated using Maximum Likelihood based methods (e.g. RAxML, PASTA). Until we support re-estimation of branch lengths within MISA, we are expecting the user to run the following  FastTree2 command explicitly before performing any placement:

`FastTreeMP -nosupport -nt -nome -noml -log tree.log -intree backbone.nwk < ref.fa > minimum_evo_backbone.nwk`

Then perform placement on the new tree:

`run_misa.py -d dist.mat -t minimum_evo_backbone.nwk`

------------------------------------
(Alternative) Installation from github
------------------------------------

If you want to install MISA from the github repository, MISA has the following dependencies which can be retrieved from `pip`:

1. `treeswift`
2. `numpy`
3. `scipy`

To install a dependency, for example `scipy`, run `pip3 install scipy`.

Clone the repository and change directory:

`git clone https://github.com/balabanmetin/misa.git && cd misa`
