# MSAvc: variation calling for genome-scale multiple sequence alignments

Fast variations calling from a FASTA/MAF file of genome-scale multiple sequence alignments



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

