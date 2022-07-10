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
./MSAvc.exe -h
```

