# Focused-Ultrasound-Simulations

## Authors 
This repository was created from the work of Dr. Charles F. Caskey's lab group at Vanderbilt University with significant contributions from Dr. Michelle K. Sigona and Dr. Thomas J. Manuel. 

## Overview 
This repository contains scripts we use to simulate focused ultrasound. Scripts included: 
- Slicer2Kwave_freeField.m : free field simulation

- Slicer2Kwave_TranscranialSteering.m : transcranial simulation with steering of a multielement array

Functions included: 
- getInputAmplitude.m : called to calculate the input amplitude used at the element surface based on a linear scaling of estimated free-field pressure. Uses 'InputAmplitudeCalibration.mat' Output is then fed into 'getFocusedElementVals.m' 

- getFocusedElementVals.m : called to calculate the amplitudes and phases on each element of a steered array.  Code based on: Ebbini and Cain, IEEE Trans. Ultrason.
Ferroelectr. Freq. Control, vol. 36, no. 5, pp. 540-8, Jan. 1989. 

- makeIGTNiftiTransducer.m : used to create a NIFTI file of your transducer, requires the coordinate locations (in mm) of your transducer. For our IGT transducer, this is the file 'IGTelemLocs.mat' 

## Requirements 
- MATLAB [Getting Started with MATLAB] (https://www.mathworks.com/help/matlab/getting-started-with-matlab.html)
- [k-wave] (http://www.k-wave.org/)
- [3D Slicer]  (https://download.slicer.org/) 


