# MSAvc: variation calling for genome-scale multiple sequence alignments

Fast variations calling from a FASTA/MAF file of genome-scale multiple sequence alignments

[![Anaconda-Server Badge](https://anaconda.org/malab/msavc/badges/version.svg)](https://anaconda.org/malab/msavc)



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

In molecular epidemiology, the typical demand for variation calling of genome-scale multiple sequence alignment (MSA) has sharply increased. However, current tools are either difficult to interpret or omit the indels and fail to handle eukaryotic genome-scale MSA. MSAvc is a C++-based program that rapidly extracts the variations including substitution, indel and replacement from multi-FASTA (prokaryotic) and multi-MAF (eukaryotic) files of genome-scale MSA. It allows users to define reference sequences for accurate variation information and filter variations of interest. MSAvc can be easily installed via Anaconda and C++ released packages on macOS, Linux, and Windows systems and is available at https://github.com/malabz/msavc.


## Installation

There are a few ways to install MSAvc. The simpliest way is using Conda. If you encounter an issue when installing MASvc or encounter a bug please report it [here](https://github.com/malabz/msavc). 
* Linux/WSL - Ubuntu/Debian
* OSX/Linux/WSL - using Bioconda
* OSX/Linux/WSL - from source
* OSX/Linux/WSL - from a release tarball
* Windows - Windows Subsystem for Linux (WSL)
* All platforms - Docker

