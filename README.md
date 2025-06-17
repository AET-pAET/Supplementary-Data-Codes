# Supplementary-Data-Codes

**Methods for accurate determination of the 3D atomic structure of amorphous materials**

Yuxuan Liao<sup>1</sup>, Haozhi Sha<sup>1</sup>, Colum M. O’Leary<sup>1</sup>, Hanfeng Zhong<sup>1</sup>, Yao Yang<sup>2</sup>, and Jianwei Miao<sup>1*</sup>

*<sup>1</sup>Department of Physics and Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.*                     
*<sup>2</sup>School of Engineering, Westlake University, Hangzhou, China.*  
<sup>*</sup>*Correspondence and requests for materials should be addressed to J.M. (j.miao@ucla.edu)*   


## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

This supplementary repository provides all data and MATLAB source codes associated with our paper, “Methods for accurate determination of the 3D atomic structure of amorphous materials.” It includes complete and reproducible workflows for robust image preprocessing, precise projection alignment, advanced tomographic reconstruction using the RESIRE algorithm, and systematic postprocessing for atomic tracing and elemental classification. These tools enable accurate 3D atomic structure determination across a range of disordered material systems. The repository supports both ADF-STEM-based atomic electron tomography (AET) and ptychographic AET (pAET), allowing researchers to reproduce key results, benchmark against published simulations, and apply the methodology to new experimental or simulated datasets. Code and data are modularly organized to facilitate adaptation to other atomic-scale imaging problems.

Due to GitHub’s file size limitations, some data file is not included directly in this repository.
Please download the large data (>100 MB) files from **[Large data file link](https://github.com/AET-pAET/Supplementary-Data-Codes/commits/v1)**.
For superlarge data files (>1 GB), please generate by running the code.  


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

# Repository Contents


### 1. NP1_aSi_208

Folder: [Si Nanoparticle 208 atoms](./NP1_aSi_208)

This folder contains the simulated data and source codes for the application of the AET workflow to a 208-atom amorphous Si nanoparticle.

### 1.1. NP1_aSi_208_From_Busch

Folder: [Si Nanoparticle 208 from Busch](./NP1_aSi_208_From_Busch)

This folder contains the simulated data and source codes used for the reconstruction of a 208-atom amorphous Si nanoparticle by Busch et al.

### 2. NP2_aSi_19993

Folder: [Si Nanoparticle 19993](./NP2_aSi_19993)

This folder contains the simulated data and source codes for the application of the AET workflow to a 19993-atom amorphous Si nanoparticle.


### 3. NP3_aSiGeSn

Folder: [SiGeSn Nanoparticle](./NP3_aSiGeSn)

This folder contains the simulated data and source codes for the application of the AET workflow to an amorphous SiGeSn nanoparticle, containing 4,610 Si, 4,709 Ge, and 4,663 Sn atoms.

### 3.1. NP3_aSiGeSn_Projection_from_Busch

Folder: [SiGeSn Nanoparticle from Busch](./NP3_aSiGeSn_Projection_from_Busch)

This folder contains the simulated data and source codes used for the reconstruction of an amorphous SiGeSn nanoparticle by Busch et al.

### 4. NP4_aCoPdPt

Folder: [CoPdPt Nanoparticle](./NP4_aCoPdPt)

This folder contains the simulated data and source codes for the application of the AET workflow to an amorphous CoPdPt nanoparticle, containing 8,322 Co, 6,896 Pd, and 3,138 Pt atoms. 

### 4.1. NP4_aCoPdPt_Projection_from_Busch

Folder: [CoPdPt Nanoparticle from Busch](./NP4_aCoPdPt_Projection_from_Busch)

This folder contains the simulated data and source codes used for the reconstruction of an amorphous CoPdPt nanoparticle by Busch et al.

### 5. NP5_aSiO2

Folder: [SiO2 Nanoparticle](./NP5_aSiO2)

This folder contains the simulated data and source codes for the application of the pAET workflow to an amorphous SiO2 nanoparticle, containing 2,603 Si and 5,101 O atoms. 

### 6. BM3D_Codes

Folder: [BM3D](./BM3D_Codes)

This folder contains the BM3D source code used for image denoising.

### 7. SIRT_Code_from_Busch

Folder: [SIRT](./SIRT_Code_from_Busch)

This folder contains the SIRT source code as used by Busch et al. in their reconstructions. 











