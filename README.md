This repository provides a description of the data analysis tools used for a Flash X-ray Imaging (FXI) experiment which was performed at 
the Linac Coherent Light Source (LCLS) and is described in 

**Daurer B.J., Okamoto K., et al. TITLE. Manuscript in preparation.**

The data has been deposited in the Coherent X-ray Imaging Data Base (CXIDB) with ID XX and can be downloaded from here: 
http://cxidb.org/id-XX.html

### List of available files: ###
File name                                           | Name     | Description
--------------------------------------------------- | -------- | ----------------------------------
http://cxidb.org/data/XX/cxidb_XX_hits.tar.gz       | **HITS** | Diffraction hits saved as CXI files.
http://cxidb.org/data/XX/cxidb_XX_background.tar.gz | **BKGR** | Diffraction background saved as CXI files.
http://cxidb.org/data/XX/cxidb_XX_metadata.tar.gz   | **META** | Auxiliary files.

### Inspecting CXI files ###
The easiest way to inspect CXI files is to use the viewing tool Owl 
(http://github.com/FXIhub/owl), but any inspection tool for HDF5 files can be used.

### Requirements ###
In order to be able to run all the provided scripts and jupyter notebooks, the following has to be installed:

* python 2.7
* numpy
* scipy
* h5py
* matplotlib
* jupyter-notebook
* libspimage (http://github.com/FXIhub/libspimage)
* condor (http://github.com/FXIhub/condor)

### Contact information
For questions about this repository and the data analysis tools, feel free to contact the author: 
* Benedikt Daurer (benedikt@xray.bmc.uu.se)

________________________________________________

## Step-by-step instructions on data analysis ##
The following sections provide a detailed description of the individual data analysis steps performed for this work. The names **HITS**, **BKGR** and **META** are used as reference to the data files which are available for download using the links above. The entire processing pipeline for this experiment looks like this:

![Overview](overview.png?raw=true)
