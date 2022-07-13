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
* [Installation](#installation)
  * [OSX/Linux/WSL \- using conda](#osxlinuxwsl---using-conda)
  * [Windows \- from released package](#Windows---from-released-package)
* [Usage](#usage)
  * [Command\-line options](#Command-line-options)
  * [Example](#example)
* [Performance](#performance)
* [Output definition](#output-definition)
* [Practical application](practical-application)
* [License](#license)
* [Feedback/Issues](#feedbackissues)
* [Citation](#citation)

## Introduction

In molecular epidemiology, the typical demand for variation calling of genome-scale multiple sequence alignment (MSA) has sharply increased. However, current tools are either difficult to interpret or omit the indels and fail to handle eukaryotic genome-scale MSA. MSAvc is a C++-based program that rapidly extracts the variations including substitution, indel and replacement from multi-FASTA (prokaryotic) and multi-MAF (eukaryotic) files of genome-scale MSA. It allows users to define reference sequences for accurate variation information and filter variations of interest. MSAvc can be easily installed via Anaconda and C++ released packages on macOS, Linux, and Windows systems and is available at https://github.com/malabz/msavc.


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
### Command-line options
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
1.Download [testdata]().

2.Run MSAvc.

```shell



```




## Performance


## Output definition


* All variations are outputted in VCF formats version 4.2, which can be parsed to standard tools such as VCF/BCFtools [11, 12] for filtering and extracting interested information or PLINK [13] for GWAS analysis. The first 8 columns are fixed, whereas the genotype information from 9 to the rest columns are depend on user’s interest, which sharply increases the file size to store the genotype matrix.

<table>
 <col class=xl65 width=47 style='mso-width-source:userset;mso-width-alt:1493;
 width:35pt'>
 <col class=xl65 width=149 style='mso-width-source:userset;mso-width-alt:4778;
 width:112pt'>
 <col width=589 style='mso-width-source:userset;mso-width-alt:18858;width:442pt'>
 <tr height=28 style='mso-height-source:userset;height:21.0pt'>
  <td colspan=3 height=28 class=xl65 width=785 style='height:21.0pt;width:800pt'>Documentary
  of fixed VCF columns</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>1</td>
  <td class=xl65>#CHROM</td>
  <td>The name of user assigned reference genome</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>2</td>
  <td class=xl65>POS</td>
  <td class=xl69>The position of the first base in reference allele according
  to the given reference genome</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>3</td>
  <td class=xl65>ID</td>
  <td>“.” meaning missing value</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>4</td>
  <td class=xl65>REF</td>
  <td>Reference allele</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>5</td>
  <td class=xl65>ALT</td>
  <td>Alternate Allele</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td height=21 class=xl65 style='height:16.0pt'>6</td>
  <td class=xl65>QUAL</td>
  <td>“.” meaning missing value</td>
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
  <td class=xl65>Sequences name*</td>
  <td class=xl66 width=589 style='width:442pt'>Given one variation, the
  genotype of sequences is assigned according the order of ALT value (if it
  matches the first alternate allele, the genotype of sequence will be assigned
  as “1” and so on)</td>
 </tr>
 <tr height=21 style='height:16.0pt'>
  <td colspan=3 height=21 class=xl67 width=785 style='height:16.0pt;width:589pt'>*The
  order of other sequences is in the same order of FASTA/MAF file after removing the chosen reference.<span
  style='mso-spacerun:yes'> 
 </tr>
 <![if supportMisalignedColumns]>
 <tr height=0 style='display:none'>
  <td width=47 style='width:35pt'></td>
  <td width=149 style='width:112pt'></td>
  <td width=589 style='width:442pt'></td>
 </tr>
 <![endif]>
</table>


#### <font color=ED7D31 face="Source Code Pro">**SUB**</font>

| Input file of single or continuous substitutions:                                                                                                                                                                                                                      |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <font face="Source Code Pro">>sequence1 (as reference)</font>                                                                                                                                                                                                          |
| <font face="Source Code Pro">TCTATCTTCGCTGCTTACGGTTTCGTCC</font>                                                                                                                                                                                                       |
| <font face="Source Code Pro">>sequence2</font>                                                                                                                                                                                                                         |
| <font color=#FF0000 face="Source Code Pro">**C**</font><font face="Source Code Pro">CT</font><font color=#FF0000 face="Source Code Pro">**CC**</font><font face="Source Code Pro">CTTCGCTGCTTACGGTTTCGT</font><font color=#FF0000 face="Source Code Pro">**GT**</font> |
| <font face="Source Code Pro">>sequence3</font>                                                                                                                                                                                                                         |
| <font color=#FF0000 face="Source Code Pro">**C**</font><font face="Source Code Pro">CT</font><font color=#FF0000 face="Source Code Pro">**GG**</font><font face="Source Code Pro">CTTCGCTGCTTACGGTTTCGT</font><font color=#FF0000 face="Source Code Pro">**GT**</font> |

the output is:
|<font size=2>#CHROM</font>|<font size=2>POS</font>|<font size=2>ID</font>|<font size=2>REF</font>|<font size=2>ALT</font>|<font size=2>QUAL</font>|<font size=2>FILTER</font>|<font size=2>INFOR</font>|<font size=2>FORMAT</font>|<font size=2>sequence1</font>|<font size=2>sequence2</font>|<font size=2>sequence3</font>|
| :-----:| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>T</font>|<font size=2>C</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=SUB;VLEN=1</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|
|<font size=2>sequence1</font>|<font size=2>4</font>|<font size=2>.</font>|<font size=2>AT</font>|<font size=2>CC,GG</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1,1;VT=SUB;VLEN=2</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>2</font>|
|<font size=2>sequence1</font>|<font size=2>27</font>|<font size=2>.</font>|<font size=2>CC</font>|<font size=2>GT</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=SUB;VLEN=2</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|

#### <font color=70AD47 face="Source Code Pro">**INS**</font>

| Input file of single or continuous insertions:                                                                                                                                                                                                                                                                  |
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <font face="Source Code Pro">>sequence1 (as reference)</font>                                                                                                                                                                                                                                                   |
| <font face="Source Code Pro">--TCTATCTTCGCTGC---TTACGGTTTCGTCC---</font>                                                                                                                                                                                                                                        |
| <font face="Source Code Pro">>sequence2</font>                                                                                                                                                                                                                                                                  |
| <font color=70AD47 face="Source Code Pro">**TT**</font><font face="Source Code Pro">TCTATCTTCGCTGC</font><font color=70AD47 face="Source Code Pro">**AAA**</font><font face="Source Code Pro">TTACGG</font><font face="Source Code Pro">TTTCGTCC</font><font color=70AD47 face="Source Code Pro">**ATT**</font> |
| <font face="Source Code Pro">>sequence3</font>                                                                                                                                                                                                                                                                  |
| <font color=70AD47 face="Source Code Pro">**TT**</font><font face="Source Code Pro">TCTATCTTCGCTGC---</font><font face="Source Code Pro">TTACGG</font><font face="Source Code Pro">TTTCGTCC</font><font color=70AD47 face="Source Code Pro">**ATT**</font>                                                      |

the output is:
|<font size=2>#CHROM</font>|<font size=2>POS</font>|<font size=2>ID</font>|<font size=2>REF</font>|<font size=2>ALT</font>|<font size=2>QUAL</font>|<font size=2>FILTER</font>|<font size=2>INFOR</font>|<font size=2>FORMAT</font>|<font size=2>sequence1</font>|<font size=2>sequence2</font>|<font size=2>sequence3</font>|
| :-----:| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
|<font size=2>sequence1</font>|<font size=2>0</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>TT</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=INS;VLEN=2</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|
|<font size=2>sequence1</font>|<font size=2>14</font>|<font size=2>.</font>|<font size=2>C</font>|<font size=2>CAAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=INS;VLEN=3</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>28</font>|<font size=2>.</font>|<font size=2>C</font>|<font size=2>CATT</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=INS;VLEN=3</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|

#### <font color=5B9BD5 face="Source Code Pro">**DEL**</font>

| Input file of single or continuous deletions:                                                                                                                                                                                                                                                           |
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <font face="Source Code Pro">>sequence1 (as reference)</font>                                                                                                                                                                                                                                           |
| <font face="Source Code Pro">TCTATCTTCGCTGCTTACGGTTTCGTCC</font>                                                                                                                                                                                                                                        |
| <font face="Source Code Pro">>sequence2</font>                                                                                                                                                                                                                                                          |
| <font color=5B9BD5 face="Source Code Pro">**--**</font><font face="Source Code Pro">TATCTTCGCTG</font><font color=5B9BD5 face="Source Code Pro">**-**</font><font face="Source Code Pro">TTACGG</font><font face="Source Code Pro">TTTCG</font><font color=5B9BD5 face="Source Code Pro">**---**</font> |
| <font face="Source Code Pro">>sequence3</font>                                                                                                                                                                                                                                                          |
| <font color=5B9BD5 face="Source Code Pro">**--**</font><font face="Source Code Pro">TATCTTCGCTGCTTACGG</font><font face="Source Code Pro">TTTCG</font><font color=5B9BD5 face="Source Code Pro">**---**</font>                                                                                          |

the output is:
|<font size=2>#CHROM</font>|<font size=2>POS</font>|<font size=2>ID</font>|<font size=2>REF</font>|<font size=2>ALT</font>|<font size=2>QUAL</font>|<font size=2>FILTER</font>|<font size=2>INFOR</font>|<font size=2>FORMAT</font>|<font size=2>sequence1</font>|<font size=2>sequence2</font>|<font size=2>sequence3</font>|
| :-----:| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=DEL;VLEN=2</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|
|<font size=2>sequence1</font>|<font size=2>13</font>|<font size=2>.</font>|<font size=2>GC</font>|<font size=2>G</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=DEL;VLEN=1</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>25</font>|<font size=2>.</font>|<font size=2>GTCC</font>|<font size=2>G</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=2;VT=DEL;VLEN=3</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>1</font>|

#### <font color=FFC000 face="Source Code Pro">**REP**</font>

| Input file of mixed situations of SUB/INS/DEL:                                                                                                                                                                                                                               |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <font face="Source Code Pro">>sequence1 (as reference)</font>                                                                                                                                                                                                                |
| <font face="Source Code Pro">---TCTATCTTCGCTGC---TTACGGTTTCGTCC---</font>                                                                                                                                                                                                    |
| <font face="Source Code Pro">>sequence2</font>                                                                                                                                                                                                                               |
| <font color=FFC000 face="Source Code Pro">**TTTAA**</font><font face="Source Code Pro">TATCTTCGCTGC</font><font color=FFC000 face="Source Code Pro">**AAAAA**</font><font face="Source Code Pro">ACGGTTTCG</font><font color=FFC000 face="Source Code Pro">**A-----**</font> |
| <font face="Source Code Pro">>sequence3</font>                                                                                                                                                                                                                               |
| <font color=FFC000 face="Source Code Pro">**------T**</font><font face="Source Code Pro">TCTTCGCT</font><font color=FFC000 face="Source Code Pro">**A-**</font><font face="Source Code Pro">---TTACGGTTTCGT</font><font color=FFC000 face="Source Code Pro">**GTATT**</font> |
| <font face="Source Code Pro">>sequence4</font>                                                                                                                                                                                                                               |
| <font color=FFC000 face="Source Code Pro">**------T**</font><font face="Source Code Pro">TCTTCGCTGC</font><font color=FFC000 face="Source Code Pro">**AAAAA**</font><font face="Source Code Pro">ACGGTTTCG</font><font color=FFC000 face="Source Code Pro">**A-----**</font> |

the output is:
|<font size=2>#CHROM</font>|<font size=2>POS</font>|<font size=2>ID</font>|<font size=2>REF</font>|<font size=2>ALT</font>|<font size=2>QUAL</font>|<font size=2>FILTER</font>|<font size=2>INFOR</font>|<font size=2>FORMAT</font>|<font size=2>sequence1</font>|<font size=2>sequence2</font>|<font size=2>sequence3</font>|<font size=2>sequence4</font>|
| :-----:| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |:----: |
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|
|<font size=2>sequence1</font>|<font size=2>1</font>|<font size=2>.</font>|<font size=2>TC</font>|<font size=2>TTTAA</font>|<font size=2>.</font>|<font size=2>.</font>|<font size=2>AC=1;VT=REP;VLEN=5</font>|<font size=2>GT</font>|<font size=2>0</font>|<font size=2>1</font>|<font size=2>0</font>|<font size=2>0</font>|



## Practical application



## License

MSAvc is free software, licensed under [MIT](https://github.com/malabz/MSAvc/blob/main/LICENSE).

## Feedback/Issues

MSAcv is supported by ZOU's Lab. If you have any questions and suggestions, please feel free to contact us on the [issues page](https://github.com/malabz/msavc/issues). You are also welcomed to send a copy to Furong.TANG@hotmail.com to make sure we could answer you ASAP! 

## Citation

If you use this software please cite:

MSAvc: variation calling for genome-scale multiple sequence alignments.[(2022)](http://....)


