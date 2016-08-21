This is a repository providing data analysis tools used for the Flash X-ray Imaging (FXI) experiment performed 
at the Linac Coherent Light Source (LCLS) which is described in 

Daurer B.J., Larsson D.S.D., et al. Advances in imaging of small virus particles at an X-ray laser. Submitted.

Data collected in this experiment has been deposited on CXIDB [http://cxidb.org](http://cxidb.org) with ID XX and can be 
downloaded from here: 

# Requirements
The data (CXI files) needs to be downloaded from the link above and a soft link has to be created inside this repository

```bash
ln -s /path/to/downloaded/data/ data
```

pointing to a directory containing all CXI files.

For the scripts to run properly, the following has to be installed
* python 2.7
* libspimage (http://github.com/FilipeMaia/libspimage)
* matplotlib
* numpy
* scipy

# Viewing the data
The easiest way to look at the data is to use the viewing tool Owl ([http://github.com/FilipeMaia/owl](http://github.com/FilipeMaia/owl)) 
which allows for detailed inspection of files in the CXI format.

Another option is to inspect data directly in an IPython session:

```ipython
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Loading 9 strongest patterns (back detector)
cxifile = 'data/cxic9714-r0163.cxi'
with h5py.File(cxifile, 'r') as f:
    hitscores = f['cheetah/event_data/peakNpix'][:]
    strongest10 = hitscores > np.sort(hitscores)[-11]
    back_detector = f['entry_1/image_1/data'][strongest10,:,:]

# Plotting
fig, arr = plt.subplots(3,3)
for i in range(9):
    im = arr[i/3,i%3].imshow(back_detector[i],vmin=15, vmax=1000)    
plt.show()
```

The most interesting data entries are summarized in the following table

