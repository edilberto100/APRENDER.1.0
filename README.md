# APRENDER.1.0

**Requirements for the installation of APRENDER**
The APRENDER pipeline is open source and can be downloaded from the project website and only is available to be compiled on Linux operating systems. For installation it is important to have a compiled and installed version of HDF5 and Cfitsio libraries, make sure you have the addresses of the libraries' location in the user's .bashrc. For the following example, the installation libraries are located in /usr/local/hdf5 and /usr/local/cfitsio.

export LD_LIBRARY_PATH="/usr/local/hdf5/lib:/usr/local/cfitsio/lib:$LD_LIBRARY_PATH"
export INCLUDE="/usr/local/hdf5/include:/usr/local/cfitsio/include:$INCLUDE"

If you do not have a version of these libraries, a version of HDF5 and Cfitiso is also made available to you on this field.

**APRENDER installation** 

1.- Download the package APRENDER.1.0.tar.gz.
2.- Unzip the package APRENDER.1.0.tar.gz
	 tar -xvf LEARN.1.0.tar.gz
3.- Enter the APRENDER folder
	 cd APRENDER/
4.- Compile the code with Make
	 make
5.- There should be no problems if the HDF5 and Cfitsio libraries are installed correctly. An executable file named APRENDER.1.0 is generated


**Execute APRENDER**

It is recommended to copy the executable file to the /usr/bin directory. APRENDER runs on a command line and to execute it is necessary to have a list of the images that are going to be processed in four different text files (bias, darks flats and Scientific Data), a list with the positions in pixels of the objects to which the photometry will be calculated in x,y format and format and you must define the diameters in pixels corresponding to the length of the annulus and dannulus.

The input parameters available for the APRENDER pipeline are as shown in the figure.
 as follows. (You can see this list do APRENDER.1.0 --help)
Usage: APRENDER.1.0 [OPTION...]  
APRENDER: Automated Pipeline for the REduction of astroNomical images and
asteroiD photomEtRy.-

  -b, --BiasList=STRING          File list that include the Bias Files
  -d, --DarkList=STRING          File list that include the Dark Files
  -f, --FlatList=STRING             File list that include the Flat Files
  -l, --PositionList=STRING      File list that include the star position x,y
  -s, --ScienceList=STRING     File list that include the Scientific Data
  -t, --TimeInt=int                      Integration time in seconds. Only if the image
                                                header does not include the integration time.
  -x, --D1=STRING            Diameters of annulus in pixels[d1,d2,d3, ... ,dn].
  -y, --D2=int                      Smallest diameter of dannulus
  -z, --D3=int                      Largest diameter of dannulus
  -?, --help                          Give this help list
      --usage                        Give a short usage message
  -V, --version                     Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

a pipeline for  TAOSII Proyect http://taos2.astrosen.unam.mx

Report bugs to <edilberto@astro.unam.mx>.
