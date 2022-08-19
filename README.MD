# Thesis Code

This repository contains code written as part of an ESAC Traineeship project in 2022. There are two repositories below. **Note:** The code below was written for specific file formats and directories and will not run otherwise. Thus it is included as an example. 

## VOTable Tools

This module contains a set of functions that can unpack **Gaia** photometry data from a RAW data link table. The functions can unpack specific sources based on their Gaia Source ID and photometric band. The functions can then plot one light curve or all the sources in an input table. The main enabler of this module is Astropy's XML parsing capabilities. 

## Light Curves

This repo contains four modules devoted to the fitting of epoch photometry data with a periodic light curve function. 

*symmetric_lc.py* contains a functions used to fit a regular symmetric sinusoidal function to a photometry time series. The code is adapted from an IDL script orignally written by F. Jimenez-Esteban from a previous 2006 study. The fitting function takes input photometry data from an ascii file with the following structure

| Time (JD)          | m (G-band)        | $\Delta m$           |
|--------------------|-------------------|----------------------|
| 2456902.119789502  | 7.58108803057825  | 0.02532920080357525  |
| 2456902.2959800074 | 7.7289054278422   | 0.006525154125789667 |
| 2456941.7547726203 | 8.097049542037922 | 0.01587221099975183  |
|         ...        |        ...        |          ...         |

