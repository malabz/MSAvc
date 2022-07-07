# MSAvc: variation calling for genome-scale multiple sequence alignments

Fast variations calling from a FASTA/MAF file of genome-scale multiple sequence alignments

[![Build Status](https://travis-ci.org/sanger-pathogens/snp-sites.png?branch=master)](https://travis-ci.org/sanger-pathogens/snp-sites)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/snp-sites/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/MGEN-doihhh-brightgreen)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/gubbins/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/gubbins)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![codecov](https://codecov.io/gh/sanger-pathogens/snp-sites/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/snp-sites)   

![preview badge](https://img.shields.io/docker/pulls/sangerpathogens/malaria-lftp)https://img.shields.io/docker/pulls/sangerpathogens/malaria-lftp

![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/malaria-lftp)

## Contents

* [Introduction](#introduction)

* [Installation](#installation)
  
  * [Linux/WSL \- Ubuntu/Debian](#linux---ubuntudebian)
  
  * [OSX/Linux/WSL \- using Bioconda](#osxlinux---using-bioconda)
  
  * [OSX/Linux/WSL \- from source](#osxlinux---from-source)
  
  * [OSX/Linux/WSL \- from a release tarball](#osxlinux---from-a-release-tarball)
  
  * [Windows \- Windows Subsystem for Linux (WSL)](#Windows---Windows-Subsystem-for-Linux-(WSL)])
  
  * [All platforms \- Docker](#all-platforms---docker)

* [Usage](#usage)
  
  * [Pipline](#pipline)
  
  * [Example](#example)
    
              [SUB](#sub)
              [INS](#ins)
              [DEL](#del)
              [REP](#rep)
  
  * [Example usage](#example-usage)
  
  * [Output](#output)

* [License](#license)

* [Feedback/Issues](#feedbackissues)

* [Citation](#citation)

## Introduction

﻿The precisive study of whole genome for every single individuals and cells has been made possible by the rapid development of sequencing technology. The demand for variation calling among individual whole genomes for prokaryotic and eukaryotic population studies would be sharply increased. Current tools are either time-consuming or less accurate and failed to handle eukaryotic genome-scale multiple sequence alignments. ﻿We present MSAvc, a C++ based program, which can efficiently extract variations including substitution, indel and replacement from FASTA and MAF files of multiple genome-scale sequence alignments. MSAvc allows users to set reference sequence for a ccurate variation extraction and filter variations based on position, alternate allele count and variation type. MSAvc is able to accurately extract variations from FASTA alignment of 3 million SARS-CoV-2 sequences and MAF alignment of 21 human whole genomes. It can be easily installed via Debian, Homebrew, Conda and Docker on MacOS, Linux, and Windows (32 and 64 bits) systems. MSAvc is open source code under the GNU General Public License version 3.

## Installation

There are a few ways to install MSAvc. The simpliest way is using apt (Debian/Ubuntu) or Conda. If you encounter an issue when installing MASvc or encounter a bug please report it [here](https://github.com/malabz/msavc). 

* Linux/WSL - Ubuntu/Debian
* OSX/Linux/WSL - using Bioconda
* OSX/Linux/WSL - from source
* OSX/Linux/WSL - from a release tarball
* Windows - Windows Subsystem for Linux (WSL)
* All platforms - Docker

### Linux/WSL - Ubuntu/Debian

If you have a recent version of Ubuntu or Debian then you can install it using apt.

```shell
apt-get install snp-sites
```

### OSX/Linux/WSL - using Bioconda

1. Install Anaconda and add the bioconda channels ([Instructional Video for command line installation](https://www.youtube.com/watch?v=AshsPB3KT-E)). Download Anaconda versions for different systems from [here](https://www.anaconda.com/products/individual-d): 64-Bit Command Line Installer (433 MB) for OSX; 64-Bit (x86) Installer (544 MB) for Linux/WSL.
   
   ```bash
   #1 Download installer
   wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
   ```

#2 Change the access permission
chmod +x Anaconda3-2021.05-Linux-x86_64.sh

#3 Run installer
./Anaconda3-2021.05-Linux-x86_64.sh

#4 Close bash and open a new shell
conda config --set auto_activate_base false
cd 
git
cd
mkdir
make
curl
bash

#5 Add channels
conda -h
conda config --show channels
conda config --add channels conda-forge
conda config --add channels bioconda

```
2. Install MSAvc through Anaconda.
```shell
#6 Creating an environment
conda env list
conda create -n bioinformatics

#7 Install MSAvc in the environment of bioinformatics
conda activate bioinformatics
conda install -c bioconda msavc
msavc -h
conda deactivate
conda remove --name moose --all
```

### OSX/Linux/WSL - from source

This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using Conda instead. First install a standard development environment (e.g. gcc, automake, autoconf, libtool). Download the software from [GitHub](https://github.com/sanger-pathogens/snp-sites).

```
autoreconf -i -f
./configure
make
sudo make install
```

### OSX/Linux/WSL - from a release tarball

This is a difficult method and is only suitable for someone with advanced unix skills. No support is provided with this method, since you have advanced unix skills. Please consider using Conda instead. First install a standard development environment (e.g. gcc, automa\
ke, autoconf, libtool).

```
tar xzvf snp-sites-x.y.z.tar.gz
cd snp-sites-x.y.z
./configure
make
sudo make install
```

### Windows \- Windows Subsystem for Linux (WSL)

#### Install Ubuntu on Windows 10

First, install Ubuntu on Windows 10, click [here](https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview) to access the instruction. [here](http://lab.malab.cn/%7Etfr/1.mp4)

Second, 

<video id="video" controls="" preload="none" height=450 width=800>
<source id="mp4" src="http://lab.malab.cn/%7Etfr/1.mp4" type="video/mp4">
</video

```
autoreconf -i
./configure
make
make check
```

This requires libcheck (the `check` package in Ubuntu) to be installed.

### All platforms - Docker

Bioconda produce a Docker container so you can use the software out of the box. Install Docker and then pull the container from Bioconda https://quay.io/repository/biocontainers/snp-sites

### 

## Usage

```
This program extracts small variation including substitution, indel and replacement from a multi FASTA/MAF alignment file, then outputs the SNP sites in VCF formats.

usage: msavc -i <inputfile> -o <outputfile> [options]
```

| Option &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Function                                                                                                                                                                                                           |
|:------------------------------------------------ |:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| -i, --in <inputfile>                             | Specify the muti-FASTA/MAF input file                                                                                                                                                                              |
| -o, --out <outputfile>                           | Specify the output VCF file name                                                                                                                                                                                   |
| -r, --reference <seqname>                        | Specify the reference genome during extracting variation (default=the first sequence of the input file)                                                                                                            |
| -g, --genotype-matrix                            | Output genotype matrix (default=false)                                                                                                                                                                             |
| -n, --nomerge-sub                                | Don't merge the SUB variations with the same "POS" and "REF" into one row (default=false)                                                                                                                          |
| -b, --filter-begin <integer>                     | The filtration of POS column by specifying an integer such as "-b 24" in terms of reference genome, meaning only keep variations with POS>=24 (default=0)                                                          |
| -e, --filter-end <integer>                       | The filtration of POS column by specifying an integer such as "--filter-end 1000" in terms of reference genome, meaning only keep variations with POS<1000 (default=18446744073709551615)                          |
| -c, --filter-ac <integer>                        | The filtration of AC tag in INFO column by specifying an integer such as "-c 100", meaning only output variations with AC>=100 (default=0)                                                                         |
| -t, --filter-vt <variationtype>                  | The filtration of VT tag in INFO column by specifying one of sub/ins/del/rep (lowercase) flag such as "-t sub", meaning only output the substitution variations (default=false)                                    |
| -l, --filter-vl <integer>                        | The filtration of VLEN tag in INFO column by specifying an integer such as "-l 25", meaning only output the variations with VLEN>=25bp (default=0)                                                                 |
| -s, --sub-block                                  | Output MSA sub-block into FASTA file after the filtration of POS column, for instance "-b 24 -e 100 -s", meaning produce a sub MSA block, the slice interval is 24=<POS<100 in terms of reference  (default=false) |
| -h, --help                                       | Help message                                                                                                                                                                                                       |
| -v, --version                                    | Version                                                                                                                                                                                                            |

### Example

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

### Example usage

```
msacv -fasta input_FASTA_filename -o output_filename.vcf

msacv -fasta input_FASTA_filename -gt -derep -subcompress -o output_filename.vcf 

msacv -fasta input_FASTA_filename -filPOS 24:296 -filAC 100 -filVT SUB -o output_filename.vcf 


msacv -maf input_MAF_filename -o output_filename.vcf

msacv -maf input_MAF_filename -gt -derep -subcompress -o output_filename.vcf

msacv -maf input_MAF_filename -filPOS 24:296 -filAC 100 -filVT SUB -o output_filename.vcf
```

### Output

* All variations are outputted in VCF formats version 4.2, which can be parsed to standard tools such as VCF/BCFtools [11, 12] for filtering and extracting interested information or PLINK [13] for GWAS analysis. The first 8 columns are fixed, whereas the genotype information from 9 to the rest columns are depend on user’s interest, which sharply increases the file size to store the genotype matrix.

<table>
 <col class=xl65 width=47 style='mso-width-source:userset;mso-width-alt:1493;
 width:35pt'>
 <col class=xl65 width=149 style='mso-width-source:userset;mso-width-alt:4778;
 width:112pt'>
 <col width=589 style='mso-width-source:userset;mso-width-alt:18858;width:442pt'>
 <tr height=28 style='mso-height-source:userset;height:21.0pt'>
  <td colspan=3 height=28 class=xl65 width=785 style='height:21.0pt;width:589pt'>Documentary
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

## License

MSAvc is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/snp-sites/blob/master/LICENSE).

## Feedback/Issues

MSAcv is supported by ZOU's Lab. If you have any questions and suggestions, please feel free to contact us on the [issues page](https://github.com/malabz/msavc/issues). You are also welcomed to send a copy to Furong.TANG@hotmail.com to make sure we could answer you ASAP! 

## Citation

If you use this software please cite:

MSAvc: variation calling for genome-scale multiple sequence alignments, Furong TANG[^#], Jiannan CHAO[^#], Fenglong Yang, Lei Xu[^*] and Quan Zou[^*], [?????(?), (2021)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000056)

[^#]: First author
[^*]: Corresponding author

