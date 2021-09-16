# img_process


## What is this repo for?
* the processing of preprocessed brain image data files for image registration
* the extraction and processing of segmented suite2p cell traces and spatial coordinates
* the extraction and processing of behavioural image data


## What does this repo contain?
* Modules contain functions which for image processing and cell trace processing
* Accompanying ipynotebooks demonstrate how to use the modules


### Modules
'admin_functions.py' - useful administrative functions useful 

'extract.py' - extracts suite2p segmented cells and runs pre-processing on it - this includes extraction and denoising of ca2+ data, binarisation with hidden markov model, and normalisation

'order_functions.py' - functions for reading and splitting large tiffs into single tiffs

'plot.py' - plotting functions for the visualisation of 

'registration.py' - functions for applying image registration to images and points, using ants library


### Notebooks

'extract.ipynb' - extracting segmented cells and processing fluorescence traces

'registration.ipynb' - register segmented cells to brain atlas and label neurons

'suite.ipynb' - running suite2p from python

'extract_tail_angle.ipynb' - extraction of tail segment positions and plotting 






