# MSAvc: variation calling for genome-scale multiple sequence alignments

Fast variations calling from a FASTA/MAF file of genome-scale multiple sequence alignments

[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/version.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/latest_release_date.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/latest_release_relative_date.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/platforms.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/license.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/downloads.svg)](https://anaconda.org/malab/msavc)
[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/installer/conda.svg)](https://anaconda.org/malab/msavc)




## Contents
* [Introduction](#introduction)
* [Pipline](#Pipline)
* [Installation](#installation)
  * [OSX/Linux/WSL \- using conda](#osxlinuxwsl---using-conda)
  * [Windows \- from released package](#Windows---from-released-package)
* [Usage](#usage)
  * [Command\-line and options](#Command-line-and-options)
  * [Example](#example)
* [Performance](#performance)
* [Output definition](#output-definition)
* [Practical application](practical-application)
* [License](#license)
* [Feedback/Issues](#feedbackissues)
* [Citation](#citation)

## Introduction

In molecular epidemiology, the typical demand for variation calling of genome-scale multiple sequence alignment (MSA) has sharply increased. However, current tools are either difficult to interpret or omit the indels and fail to handle eukaryotic genome-scale MSA. MSAvc is a C++-based program that rapidly extracts the variations including substitution, indel and replacement from multi-FASTA (prokaryotic) and multi-MAF (eukaryotic) files of genome-scale MSA. It allows users to define reference sequences for accurate variation information and filter variations of interest. MSAvc can be easily installed via Anaconda and C++ released packages on macOS, Linux, and Windows systems and is available at https://github.com/malabz/msavc.

## Pipline
Four standard types of small variation (Danecek et al., 2011): SUB (Substitution) represents single/multiple nucleotide substitutions; INS (Insertion) represents single/multiple insertions; DEL (Deletion) represents single/multiple deletions, REP (Replacement) stands for the complex event that the co-occurrence of SUB, INS or DEL.

![VT](http://lab.malab.cn/%7Etfr/MSAvc_testdata/pipline1.svg)


## Installation

There are a few ways to install MSAvc. The simpliest way is using Conda. If you encounter an issue when installing MSAvc or encounter a bug please report it [here](https://github.com/malabz/MSAvc/issues). 
### OSX/Linux/WSL - using conda
1.Intall WSL for Windows. Instructional video [1](https://www.youtube.com/watch?v=X-DHaQLrBi8&t=5s) or [2](http://lab.malab.cn/%7Etfr/1.mp4) (Copyright belongs to the original work).

2.Download and install Anaconda. Download Anaconda versions for different systems from [here](https://www.anaconda.com/products/distribution#Downloads). Instructional video of anaconda installation [1](https://www.youtube.com/watch?v=AshsPB3KT-E) or [2](http://lab.malab.cn/%7Etfr/Install_anaconda_in_Linux.mp4) (Copyright belongs to the original work).

3.Install MSAvc.

```bash
#1 Acvtivate one of you conda environment
conda activate base

#2 Add channels to conda
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels malab

#3 Install the required package boost-cpp=1.77.0 for running msavc
conda install -c conda-forge boost-cpp=1.77.0

#4 Install msavc
conda install -c malab msavc

#5 Test msavc
msavc -h
```

### Windows - from released package
1. Download MSAvc from [relseases](https://github.com/malabz/MSAvc/releases/new).
2. Test msavc
```
.\msavc.exe -h
```

## Usage
### Command-line and options
#### Command line
```
 msavc -i <inputfile> -o <outputfile> [options]               # conda version
 
 .\msavc.exe -i <inputfile> -o <outputfile> [options]         # windows version
```
#### Options
```
   -i, --in <inputfile>              Specify the muti-FASTA/MAF input file
   -o, --out <outputfile>            Specify the output VCF file name

   -r, --reference <seqname>         Specify the reference genome during extracting 
                                     variations (default=the first sequence of the 
                                     input file)

   -g, --genotype-matrix             Output genotype matrix (default=off)

   -n, --nomerge-sub                 Do not merge the SUB variations with the same 
                                     "POS" and "REF" into one row (default=off)

   -b, --filter-begin <integer>      Filtration of the POS column by specifying an 
                                     integer such as "-b 24" in terms of the 
                                     reference genome, meaning only keep variations 
                                     POS>=24 (default=1, 1-based index)

   -e, --filter-end <integer>        Filtration of the POS column by specifying an 
                                     integer such as "-e 1000" in terms of the 
                                     reference genome, meaning only keep variations 
                                     with POS<=1000 (default=last base index)

   -c, --filter-ac <integer>         Filtration of the AC tag in the INFO column by 
                                     specifying an integer such as "-c 100", meaning 
                                     only output variations with AC>=100 (default=0)

   -t, --filter-vt <variationtype>   Filtration of the VT tag in the INFO column by 
                                     specifying one of sub/ins/del/rep (lowercase) 
                                     flags such as "-t sub", meaning only output the 
                                     substitution variations (default=off)

   -l, --filter-vl <integer>         Filtration of the VLEN tag in the INFO column 
                                     by specifying an integer such as "-l 5", 
                                     meaning only output the variations with 
                                     VLEN>=5bp (default=0)

   -s, --sub-block                   Output MSA sub-block into FASTA file, "-s" 
                                     option works only when "-b" and "-e" are both 
                                     specified, for instance "-b 24 -e 1000 -s", 
                                     meaning produce a sub MSA block, the slice 
                                     interval is 24=<POS<=1000 in terms of the 
                                     reference genome (default=off)

   -f, --force-overwrite             Overwrite existing file (default=off)
   -h, --help                        Help message
   -v, --version                     Version


```
### Example
1.Download testdata.
- FASTA format : <a href="http://lab.malab.cn/%7Etfr/MSAvc_testdata/h3_100w.tar.xz" download="h3_100w.tar.xz">h3_100w.tar.xz</a> is the alignment of 1 million respiratory syndrome coronavirus 2 (SARS‑CoV‑2) genomes via the MSA tool [HAlign 3](https://github.com/malabz/HAlign-3).
- MAF format: ... is the alignment of 3 human chromosone 1.

2.Run MSAvc.

```shell



```




## Performance


## Output definition
### Variation type
Four standard types of small variation (Danecek et al., 2011): SUB (Substitution) represents single/multiple nucleotide substitutions; INS (Insertion) represents single/multiple insertions; DEL (Deletion) represents single/multiple deletions, REP (Replacement) stands for the complex event that the co-occurrence of SUB, INS or DEL.

![VT](http://lab.malab.cn/%7Etfr/MSAvc_testdata/new/sub.svg)


### Variant Call Format

All variants will be output in standard VCF format version 4.1, which can be the input of software such as VCFtools, BCFtools for filtering or extracting variants of interest, and PLINK for GWAS analysis. The first 8 columns of the VCF file are mandatory and fixed with variation information, while the 9th to the rest columns contain the genotypes for each genome.

<table>
 <col class=xl65 width=47 style='mso-width-source:userset;mso-width-alt:1493;
 width:35pt'>
 <col class=xl65 width=149 style='mso-width-source:userset;mso-width-alt:4778;
 width:112pt'>
 <col width=800 style='mso-width-source:userset;mso-width-alt:18858;width:442pt'>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>1</td>
  <td class=xl65>#CHROM</td>
  <td>The name of user assigned reference genome</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>2</td>
  <td class=xl65>POS</td>
  <td class=xl69>The position of the first base in the reference allele according
  to the given reference genome</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>3</td>
  <td class=xl65>ID</td>
  <td>Semi-colon separated list of unique identifiers where available for variants. “.” missing value was used in MSAvc.</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>4</td>
  <td class=xl65>REF*</td>
  <td>Reference allele</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>5</td>
  <td class=xl65>ALT*</td>
  <td>Comma separated list of alternate non-reference alleles</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>6</td>
  <td class=xl65>QUAL</td>
  <td>Quality score for the assertion made in ALT for variants extracted from BAM file. “.” missing value was used in MSAvc.</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>7</td>
  <td class=xl65>FILTER</td>
  <td>“.” meaning missing value</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td rowspan=3 height=63 class=xl65 style='height:48.0pt'>8</td>
  <td rowspan=3 class=xl65>INFO</td>
  <td>AC: Alternate Allele Count</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 style='height:16.0pt'>VT: the type of variation (SUB, INS, DEL
  and REP)</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 style='height:16.0pt'>VLEN: Difference in length between REF
  and ALT alleles</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>9</td>
  <td class=xl65>FORMAT</td>
  <td>GT: genotype</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>10</td>
  <td class=xl65>Reference name</td>
  <td>The genotypes of reference genome which will always be “0”</td>
 </tr>
 <tr height=67 style='mso-height-source:userset;height:50.0pt'>
  <td height=67 class=xl65 style='height:50.0pt'>11</td>
  <td class=xl65>Sequences name**</td>
  <td class=xl66 width=800 style='width:442pt'>Given one variation, the
  genotype of sequences is assigned according the order of ALT value (if it
  matches the first alternate allele, the genotype of sequence will be assigned
  as “1” and so on)</td>
 </tr>
 </table>

 *For INS, DEL and REP variations, the REF and ALT Strings must include the base before the event (which must be reflected in the POS field), unless the event occurs at position 1 on the reference genome in which case it must include the base after the event.

 **The order of other sequences is in the same order of FASTA/MAF file after removing the chosen reference.




## Practical application



## License

MSAvc is free software, licensed under [MIT](https://github.com/malabz/MSAvc/blob/main/LICENSE).

## Feedback/Issues

MSAcv is supported by ZOU's Lab. If you have any questions and suggestions, please feel free to contact us on the [issues page](https://github.com/malabz/msavc/issues). You are also welcomed to send a copy to Furong.TANG@hotmail.com to make sure we could answer you ASAP! 

## Citation

If you use this software please cite:

MSAvc: variation calling for genome-scale multiple sequence alignments.[(2022)](http://....)


