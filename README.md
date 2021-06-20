# GenericParallelFramework
Generic Parallel Framework for Inferring Large Scale Gene Regulatory Networks from Expression Profiles

This repo contains the code to implement the generic parallel framework proposed by us on the 15 inference algorithms - RELNET, ARACNE, CLR, C3NET, C3MTC, BC3NET, MutRank, MINE, PCOR, SPCOR, GeneNet, G1DBN, MRNET, MRNETB and GENIE3. You may also add inference methods of your choice to infer them parallely.

## Citation information

If you use this code in your research, please cite this repo and/or the paper (link will be added soon)

## Usage

To use this code first install R and the packages mentioned in the beginning of each file. After installation of R and the necessary packages, enter the path of the dataset and the path for the output (has been mentioned in the code) and then enter the number of divisions you want and run the entire file. 

## Dependencies
1. The different files in this code will require a different version of R or Bioconductor to execute properly.

There are chances that R will throw you an error if the version of R and the version of the package don't match so be careful while downloading the packages. Some packages might require R 4.0 or higher while others may require R 3.6 or lower. 

2. The R parallel package works properly only for LINUX version for R. Therefore it is advisable to execute these codes in LINUX R environment. However, the files- GeneNet and Genie3 will only execute in Windows R environment.
