PVOCEX sources v1.0.


These have the option to use the libraries from FFTW (http://fftw.org). Anyone wanting to build from the
sources will need to download and build the FFTW libraries, according to the instructions in the FFTW documentation. 
The current code is for FFTW v 2.1.5.

Note that PVOCEX requires the 'floats' version of the FFTW libraries. 
The path to the libraries in the Makefiles will then need to be edited to reflect where you have installed FFTW. 

The FFTW libraries are used when USE_FFTW is defined (see the Makefiles for other necessary changes). 

The default is the CDP fft functions (converted from Fortran) in mxfft.c.
These are still fast (for example, used in 'pvplay'), and may well be more than sufficient for most purposes. 


Richard Dobson 2000, 2022
