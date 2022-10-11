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
* [Practical application](practical-application)
* [Output definition](#output-definition)
* [License](#license)
* [Feedback/Issues](#feedbackissues)
* [Citation](#citation)

## Introduction

In molecular epidemiology, the typical demand for variation calling of genome-scale multiple sequence alignment (MSA) has sharply increased. However, current tools are either difficult to interpret or omit the indels and fail to handle eukaryotic genome-scale MSA. MSAvc is a C++-based program that rapidly extracts the variations including substitution, indel and replacement from multi-FASTA (prokaryotic) and multi-MAF (eukaryotic) files of genome-scale MSA. It allows users to define reference sequences for accurate variation information and filter variations of interest. MSAvc can be easily installed via Anaconda and C++ released packages on macOS, Linux, and Windows systems and is available at: https://github.com/malabz/msavc.

## Pipline
Four standard types of small variation (Danecek et al., 2011): SUB (Substitution) represents single/multiple nucleotide substitutions; INS (Insertion) represents single/multiple insertions; DEL (Deletion) represents single/multiple deletions, REP (Replacement) stands for the complex event that the co-occurrence of SUB, INS or DEL.

![VT](http://lab.malab.cn/%7Etfr/MSAvc_testdata/pipline2.svg)


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
   -i, --in <inputfile>              Specify the multi-FASTA/MAF input file
   -o, --out <outputfile>            Specify the output VCF file name

   -r, --reference <seqname>         Specify the reference genome during extracting variations 
                                     (default=the first sequence of the input file)

   -g, --genotype-matrix             Output genotype matrix (default=off)

   -n, --nomerge-sub                 Do not merge the SUB variations with the same "POS" and "REF" 
                                     into one row (default=off)

   -b, --filter-begin <integer>      Filtration via the POS column by specifying an integer in terms 
                                     of the reference genome: "-b 24" keep the variations with POS>=24 
                                     (default=1, 1-based index)

   -e, --filter-end <integer>        Filtration via the POS column by specifying an integer: "-e 1000"
                                     keep the variations with POS<=1000 (default=the last base index, 
                                     1-based index)

   -c, --ac-greater <integer>        Filtration via the AC tag in the INFO column by specifying an 
                                     integer: "-c 10" output the variations with AC>=10 (default=0)

   -d, --ac-less <integer>           Filtration via the AC tag: "-d 100" output the variations with 
                                     AC<=100 (default=the total number of sequences)

   -t, --filter-vt <variationtype>   Filtration via the VT tag in the INFO column by specifying one of 
                                     sub/ins/del/rep (lowercase) flags: "-t sub" only output the 
                                     substitution variations (default=off)

   -l, --vl-greater <integer>        Filtration via the VLEN tag in the INFO column by specifying an 
                                     integer: "-l 5" output the variations with VLEN>=5bp (default=0)

   -m, --vl-less <integer>           Filtration via the VLEN tag: "-m 10" output the variations with 
                                     VLEN<=10bp (default=length of reference genome)

   -s, --sub-block                   Output MSA sub-block into FASTA file, "-s" option works only when 
                                     "-b" and "-e" are both specified: "-b 24 -e 1000 -s" produce a 
                                     sub MSA block; the slice interval is 24=<POS<=1000 in terms of 
                                     the reference genome (default=off)

   -C, --compress-bgz                Compress the VCF output file. As the number of sequences and
                                     variations increases, the VCF file with "-g" becomes super large.

   -B, --buffer-size                 String buffer size: For large files, set this size smaller in 
                                     order to lower the memory usage, only used in FASTA file 
                                     (default=1048576)

   -N, --duplicate-name              Check duplicate name in file (default=off)
   -f, --force-overwrite             Overwrite existing file (default=off)
   -h, --help                        Help message
   -v, --version                     Version


```
### Example
1.Download testdata.
- FASTA format : <a href="http://lab.malab.cn/%7Etfr/MSAvc_testdata/halign3_sars_cov_2_10kseq.tar.xz" download="halign3_sars_cov_2_10kseq.tar.xz">halign3_sars_cov_2_10kseq.tar.xz</a> is the alignment of 10 thousand respiratory syndrome coronavirus 2 (SARS‑CoV‑2) genomes via the MSA tool [HAlign 3](https://github.com/malabz/HAlign-3). The reference genome (EPI_ISL_402124) is the first one. 
```
file                           format  type  num_seqs      sum_len  min_len  avg_len  max_len
halign3_sars_cov_2_10kseq.fas  FASTA   DNA     10,000  431,940,000   43,194   43,194   43,194
```
- MAF format: <a href="http://lab.malab.cn/%7Etfr/MSAvc_testdata/parsnp_human_chr1_21seq.tar.xz" download="parsnp_human_chr1_21seq.tar.xz">parsnp_human_chr1_21seq.tar.xz</a> is the alignment of 21 human chromosone 1 via the MSA tool [Parsnp](https://github.com/marbl/parsnp). The reference genome (GRCh38.p13) is the first one. 

2.Run MSAvc.

```shell



```




## Performance



## Practical application


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





## License

MSAvc is free software, licensed under [MIT](https://github.com/malabz/MSAvc/blob/main/LICENSE).

## Feedback/Issues

MSAcv is supported by ZOU's Lab. If you have any questions and suggestions, please feel free to contact us on the [issues page](https://github.com/malabz/msavc/issues). You are also welcomed to send a copy to Furong.TANG@hotmail.com to make sure we could answer you ASAP! 

## Citation

If you use this software please cite:

MSAvc: variation calling for genome-scale multiple sequence alignments.[(2022)](http://....)


