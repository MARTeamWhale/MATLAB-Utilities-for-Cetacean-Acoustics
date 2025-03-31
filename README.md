# MATLAB Utilities for Cetacean Acoustics (MUCA)


## Overview
This library contains various MATLAB functions and classes that are helpful for analyzing cetacean acoustic data. It is organized as a monolithic MATLAB package divided into subpackages. The subpackages currently consist of the following:

- ***algorithms*** - functions that perform miscellaneous data manipulation tasks, such as breaking up a vector into smaller pieces
- ***audio*** - tools for working with and extracting data from audio files
- ***dcs_analysis*** - tools for evaluating the performance of signal detection and/or classification systems
- ***dsp*** - functions and classes to perform digital signal processing tasks, such as creating and manipulating spectrograms
- ***filepaths*** - tools for locating files and manipulating file and folder path strings
- ***io*** - tools for reading or writing common files (i.e., text or image files)
- ***time*** - functions for working with dates and times

Note that some of the code within the MUCA library is interdependent; that is, some functions may depend on other functions within MUCA. Therefore, it is not a good idea to download individual files within this library. The library should be cloned in its entirety.


## Installation
- Clone or download the entire repository *MATLAB-Utilities-for-Cetacean-Acoustics* and save it somewhere on your computer
- Add the location of the cloned repository folder to your MATLAB path (there is no need to add subfolders)


## Usage
The MUCA library uses MATLAB's packaging system, which consists of folders prefixed with the "+" symbol.

To access the MUCA code within MATLAB, use either dot notation or the "import" statement. For example:

	Using dot notation:
		dt = MUCA.time.readDateTime(fileNames)
	
	Using an import statement:
		import MUCA.time.readDateTime
		dt = readDateTime(fileNames)
	
See the MATLAB documentation on packages and *import* for more information.

Refer to the headers within each M-file to learn what they do and how to use them (tip: it's also possible to see headers without opening files by using the *help* command in MATLAB followed by the M-file name).
