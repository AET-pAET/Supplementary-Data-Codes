# Supplementary-Data-Codes


- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

# System Requirements

## Hardware Requirements

We recommend a computer with 32G DRAM, AMD Ryzen 7 6800H with Radeon Graphics 3.20 GHz CPU, and a NVIDIA GeForce RTX 3070 Ti GPU to run most data analysis source codes. But for the atomic tracing, we recommend a computer with 64G DRAM.

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: CentOS 6 2.6.32
Windows: Windows 11 Version 24H2 
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.   

### Matlab Version Requirements

This package has been tested with `Matlab` R2022b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2021a or higher to test the data and source codes.

# Repositary Contents


### 1. NP1_aSi_D2nm

Folder: [Si Nanopartile 2nm](./NP1_aSi_D2nm)

This folder contains the simulation data of a Si nanoparticle with size of 2nm ... (add other content to describe the folder)


#### 1.1 projections

Folder: [Projections](./NP1_aSi_D2nm/1.1_projections)

This folder contains the projections after alignment as well as their corresponding angles.

#### 1.2 Reconstruction

Folder: [Projections](./NP1_aSi_D2nm/1.2_reconstructions)

Run the tomography reconstruction code to obtain the 3D reconstruction of the nanoparticles.

#### 1.3 Atom Tracing

Folder: [Projections](./NP1_aSi_D2nm/1.3_tracing)

Run the polynomial atom tracing code to obtain the 3D coordinates of atoms in the nanoparticles.























