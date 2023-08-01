# SCOS data processing

Communications Biology 2023 - Measuring human cerebral blood flow and brain function with fiber-based speckle contrast optical spectroscopy system

We created a data analysis pipeline for SCOS under photon starved conditions. This repository contains:

1) Sample speckle images and dark camera image from Hamamatsu ORCA-Fusion camera (./sample data)
2) Camera parameters such as gain and mask (./camera calibration)
3) Functions and script to run the noise correction pipeline

To run the code:

1) Download all the files and subfolders in one folder
2) Run HamamatsuPipeline.m

The HamamatsuPipeline.m outputs: 

1) Processed data files (img_file_mean_img.mat, analyzed_K2_img_file_YYYY_MM_DD_HH_mm.mat) in ./sample data.
2) Raw standard deviation and mean intensity time course   
   <img width="537" alt="Screenshot 2023-08-01 at 1 01 13 PM" src="https://github.com/BUNPC/2023-SCOS/assets/55467463/4ddceb44-6f8d-41f6-8597-b59867c97461">
3) Linear estimation of window I(t) in the contrast measurement
   <img width="540" alt="Screenshot 2023-08-01 at 1 01 35 PM" src="https://github.com/BUNPC/2023-SCOS/assets/55467463/6443c906-7d1e-40f4-bada-4362fe3dd914">
4) Contribution of noise sources to the measured speckle contrast
   <img width="539" alt="Screenshot 2023-08-01 at 1 01 53 PM" src="https://github.com/BUNPC/2023-SCOS/assets/55467463/d70cc54a-b65f-4537-9e13-c04fab0c9bad">
   
5) Fundamental speckle contrast after noise correction and the corresponding blood flow index

   <img width="543" alt="Screenshot 2023-08-01 at 1 02 06 PM" src="https://github.com/BUNPC/2023-SCOS/assets/55467463/578c7c4f-e610-4b15-a806-f895a126be78">

For any issue reporting or suggestions, please contact Byungchan Kim, kennykim@bu.edu


[![DOI](https://zenodo.org/badge/672077193.svg)](https://zenodo.org/badge/latestdoi/672077193)

