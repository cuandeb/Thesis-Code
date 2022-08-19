# Thesis Code

This repository contains code written as part of an ESAC Traineeship project in 2022. There are two repositories below. **Note:** The code below was written for specific file formats and directories and will not run otherwise. Thus it is included as an example. 

## VOTable Tools

This module contains a set of functions that can unpack **Gaia** photometry data from a RAW data link table. The functions can unpack specific sources based on their Gaia Source ID and photometric band. The functions can then plot one light curve or all the sources in an input table. The main enabler of this module is Astropy's XML parsing capabilities. 

## Light Curves

This repo contains four modules devoted to the fitting of Long Period Variable epoch photometry data with a periodic light curve function. 

*symmetric_lc.py* contains a functions used to fit a regular symmetric sinusoidal function to a photometry time series. The code is adapted from an IDL script orignally written by F. Jimenez-Esteban from a previous 2006 study. The fitting function takes input photometry data from an ascii file with the following structure

| Time (JD)          | m (G-band)        | $\Delta m$           |
|--------------------|-------------------|----------------------|
| 2456902.119789502  | 7.58108803057825  | 0.02532920080357525  |
| 2456902.2959800074 | 7.7289054278422   | 0.006525154125789667 |
| 2456941.7547726203 | 8.097049542037922 | 0.01587221099975183  |
|         ...        |        ...        |          ...         |


*asymmetric_lc.py* performs the same task but fits the data with an asymmetric cosine function, typical of Mira variable stars. 

*plotting_lc.py* plots the result of the fitting performed by the two previous modules. Four plots are generated: A plot of the epoch photometry, a plot of the fitted curve, a periodogram ($$\chi^2$$ vs. period for the test periods), and a folded light curve constructed from the period of best fit. An example of such a plot is shown below

![image](https://user-images.githubusercontent.com/90631697/185609355-5ef14008-2a97-4e48-84e2-a1d0457d06c2.png)

Finally, *ALC_fitting.py* combines all of the three previous modules. Given a list of ascii files with the structure shown above it will fit a light curve, pass the fit parameters to the plotting function and outputs the 4 plots. The function includes an option to first fit a sinusoidal light curve to the data and pass these parameters to the asymmetric fitting function which will then perform an abridged fit using the input parameters to constrain the new fit. 
