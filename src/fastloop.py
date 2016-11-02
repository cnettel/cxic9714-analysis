"""
A module for a fast looping through a dataset (using python's multiprocessing)

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se) and Max Hantke
Last update:  March 9, 2014
"""
import multiprocessing
import numpy as np
import h5py, os, sys, time
import utils
__all__ = ['FastLoop']


def _worker_call(worker, pipe):
    work_package = pipe.recv()
    while work_package != "fika":
        pipe.send(worker(work_package))
        work_package = pipe.recv()

class FastLoop(object):
    def __init__(self, inputfile, outputfile, Nprocesses, Njobs, loopfunc, updaterate, *args, **kwargs):
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.index = 0
        self.Nprocesses = Nprocesses
        self.Njobs = Njobs
        self.framerate_update = updaterate
        self.loopfunc = loopfunc
        self.args = args
        self.kwargs = kwargs
        self.load()

    def load(self):
        with h5py.File(self.inputfile, 'r') as datafile: self.input = {item:datafile[item][:] for item in datafile.keys()}
        with h5py.File(self.outputfile, 'r') as datafile: self.output = {item:datafile[item][:] for item in datafile.keys()}
        #pass

    #@profile
    def next(self):
        self.index += 1
        #with h5py.File(self.inputfile, 'r') as datafile: data = datafile['data'][self.index-1]
        data = {item:self.input[item][self.index-1] for item in self.input}
        input = {'i':self.index-1, 'data':data, 'args':self.args, 'kwargs':self.kwargs}
        return input
    
    #@profile
    def save(self, res):
        for k in res:
            #with h5py.File(self.outputfile, 'a') as file:
            #    if not k == 'i':
            #        file[k][res['i']] = res[k]
            if not k == 'i': 
                self.output[k][res['i']] = res[k]
                #if res[k].ndim == 0:
                #    self.output[k][res['i']] = res[k]
                #else:
                #    self.output[k][res['i']][...] = res[k]

    #@profile
    def write(self, progress = 1.):
        utils.progressbar(progress, 'Saving results to file')
        with h5py.File(self.outputfile, 'a') as datafile:
            for item in self.output: datafile[item][...] = self.output[item]

    #@profile
    def start(self):        
        
        # Decide how many parallel processes to use
        available_cpus = multiprocessing.cpu_count()
        if self.Nprocesses == 'max': self.Nprocesses = available_cpus
        elif self.Nprocesses > available_cpus: self.Nprocesses = available_cpus
    
        # Initialize lists keeping track of read/write pipes
        pipes_end_host = list(np.zeros(self.Nprocesses))
        pipes_end_worker = list(np.zeros(self.Nprocesses))
        processes = list(np.zeros(self.Nprocesses))

        # Start multiple processes
        for i in range(self.Nprocesses):
            pipes_end_host[i], pipes_end_worker[i] = multiprocessing.Pipe()
            processes[i] = multiprocessing.Process(target=_worker_call, args=(self.loopfunc, pipes_end_worker[i],) )
            processes[i].start()

        # Variables to keep track of jobs started/done
        Njobs_done = 0            
        Njobs_started = 0

        # Send initial jobs
        for r in pipes_end_host: 
            r.send(self.next())
            Njobs_started += 1

        # Some parameters for diagnostics
        t_start = time.time()
        message = 'Datarate %.2f Hz; job %i/%i; process'        

        # This is the main loop, it waits for new jobs and sends the input
        while Njobs_done < self.Njobs:
            for r in pipes_end_host:
                if r.poll():
                    result = r.recv()
                    if Njobs_started < self.Njobs: 
                        r.send(self.next())
                        Njobs_started += 1
                    self.save(result)
                    Njobs_done += 1                    

                    # Give some feedback to the command line once in a while
                    if ((Njobs_done + 1) % self.framerate_update) == 0:
                        progress = float(Njobs_done) / self.Njobs
                        datarate = (Njobs_done + 1) / (time.time() - t_start)
                        utils.progressbar(progress, message %(datarate, Njobs_done + 1, self.Njobs), t_start)

        # Close all processes
        for i in range(self.Nprocesses):
            pipes_end_host[i].send('fika')
            processes[i].join()
            pipes_end_host[i].close()

# ===================================
# Some code for testing/profiling
# ===================================
if __name__ == '__main__':

    from fitting import fit_photon_histograms
    import logging
    logging.captureWarnings(True)

    histfile = h5py.File('testing/histogram.h5', 'r')
    Hmap = histfile['data/histogram']
    Hbinsize = histfile['data/histogramBinsize'][0]
    Hcount = histfile['data/histogramCount'][0]
    Hmin = histfile['data/histogramMin'][0]
    Hnbins = histfile['data/histogramNbins'][0]
    Hbins = np.linspace(Hmin, Hnbins + Hmin - 1, Hnbins)
    NY = Hmap.shape[0]
    NX = Hmap.shape[1]
    SH = (NY, NX)
    NXY = NX * NY

    infile = h5py.File('testing/tmpin.h5', 'w')
    infile['histogram'] = Hmap[:].reshape(NXY, Hnbins)
    infile['mask'] = np.array(NXY * [False])
    infile.close()
    histfile.close()
    
    outfile = h5py.File('testing/tmpout.h5', 'w')
    outfile.create_dataset('bg_offset', (NXY,))
    outfile.create_dataset('bg_amp', (NXY,))
    outfile.create_dataset('bg_sigma', (NXY,))
    outfile.create_dataset('photon_offset', (NXY,))
    outfile.create_dataset('photon_amp', (NXY,))
    outfile.create_dataset('photon_sigma', (NXY,))
    outfile.create_dataset('status', (NXY,), dtype=h5py.special_dtype(vlen=str))
    outfile.close()
    
    fastloop = FastLoop('testing/tmpin.h5', 'testing/tmpout.h5', 30, NXY, fit_photon_histograms, Hbins)
    fastloop.start()
    fastloop.write()
    
    os.remove('testing/tmpin.h5')
    os.remove('testing/tmpout.h5')



