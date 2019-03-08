# GDGT_Models

The file AllFunctions.R contains a variety of functions to implement the model described in Eley, Y. et al. A Gaussian Process calibration for GDGT-based palaeothermometry.

The rest of the files contain either pre-trained models or modern calibration data.

You will need to ensure that you have a version of Python installed (the code is based on Python 3.6) along with the GPy library (https://sheffieldml.github.io/GPy/).

The code uses the reticulate R package to make use of GPy's extensive Gaussian Process framework within R, so you will also need to install reticulate and set it up by pointing to your Python distribution via 
Sys.setenv(RETICULATE_PYTHON = path_to_your_python_distribution)
