# img_process
This repo contains python code for pre-processing of image data. 

Modules contain functions which extract, filter and register calcium imaging data. 
Accompanying ipynotebooks demonstrate how to use the modules. 

admin_functions - simple functions ordering and naming data.
plot - plotting function for data processing visualisation
cellGUI - widget to visualise individual neurons in TIFF recording - work in progress.=
order_functions - functions for reading and splitting large tiffs into single tiffs
suite - running suite2p from python
registration - functions for applying image registration to images and points, using ants library. 
deeplabcut - runs convnet DLC to label tail data. 
extract - extracts suite2p data, and runs pre-processing on it - this includes extraction and denoising of ca2+ data, binarisation with hidden markov model, and normalisation. 





