# 				CDP PVXTOOLS

This is a first release of a small toolkit demonstrating the PVOCEX file format for phase vocoder (pvoc) analysis files. It is intended to develop this resource further over time.

The PVOCEX file format (extension .pvx) was developed in 2000 and presented at the ICMC conference that year. It was added to Csound soon afterwards, in which it is now in very general use, superceding the original non-portable "pv" file format.

It is now also fully supported in the new version of the CDP system ("Release 8" - see the repository "CDP8"), alongside the long established ".ana" format based somewhat too closely on the standard WAVE format.

## The programs

The programs (and libraries) are provided as source files in the usual way of github, but to encourage immediate explorations a few precomiled binaries for PC and Mac (64bit) are also provided, in appropriate zip files. Run any program without arguments to see the usage message.


*  ANA2PVX

	A simple program to convert a CDP .ana file into a mono .pvx file.

* PVOCEX

	A fairly straight port of the original CARL phase vocoder (Mark Dolson) - takes a mono soundfile (.wav, .aiff) and outputs a mono .pvx file. See the supplied documentation for details on usage. Long-standing CDP users will recognise this as the original "pvoc" program as supplied with the first Atari-based systems.
	
* PVOCEX2

	Similar to PVOCEX, but with a slightly different arrangement of flag options - these include a -I flag to print the analysis properties of a .pvx file. This version of the CARL phase vocoder can analyse a stereo input file to a stereo PVOCEX file (which **pvplay** can play). When processing a .pvx file back to audio, there is an option to use an oscillator bank for resynthesis instead of the usual inverse FFT. This program is written in elementary C++ wrapping the original C code.

* PVPLAY

	A multi-option command line player program for both analysis files and (multichannel) soundfiles.  Type the name without arguments to see all options.	
* DIRSF

	This is CDP's workhorse directory listing program. Used without arguments, it lists all CDP-recognised soundfiles in the current directory, with information on format, duration, sample type, etc. This now includes information on .pvx analysis files. Use with a ? argument to see options.

* SPECVU

	A simple program to print the amp/freq values of a nominated window of an analysis file to a text file, or (Mode 2)  max/min amp and freq of the whole file. Type the program name without arguments to see a usage message. Note that a wildcard facility (on filenames) is available for a selective listing.
	
* PVXIO

	The basic library to read and write .pvx files.
	
* PORTSF

	The library to read/write soundfiles (.wav, .aiff), from the CDP M/C Toolkit and the MIT Audio Programming Book.
	

### Examples

A few short example PVOCEX analysis files (mono and stereo) are available in the **examples** directory. We may update/replace these, for fun, at irregular intervals. Download ad lib., and play with **pvplay**; convert back to soundfiles with **pvocex**/**pvocex2**.

### To Come:

all sorts of things, including conversions to use CMake as per the CDP system.

Check the directory **misc** periodically for assorted pvoc-related code (old and new).

-------------------------------------------------------------------------------------
Richard Dobson 

last revision: 02/01/2023