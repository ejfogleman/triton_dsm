# triton_dsm.py
# matlab plotting scripts translated to python
# Eric Fogleman -- 8/2022
#
# usage:  import triton_dsm
# to re-import
# import importlib
# importlib.reload(triton_dsm)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

plt.ioff()  # disable interactive mode

def plot_dsm_psd_linear_freq( datfile, delta, Fs, d, OSR, pltres ):
    """ 
    triton_dsm.plot_dsm_psd_linear_freq
	Plot psd in dB; scale based on spectrum; linear x axis
       Provides 2x scaling to give V^2/Hz for + freq

	Usage:
	plot_dsm_psd_linear_freq( datfile, delta, Fs, d, OSR, pltres )

	Where:
        datfile = name of file containing output from dsm.c
        datfile(:,1) = normalized frequency vector (cycles/sample)
	datfile(:,2) = psd (codes^2/(cycles/sample))
	datfile(:,3) = res psd (codes^2/(cycles/sample))
 	delta = lsb size (V)
 	Fs = sample rate (Hz)
        OSR = oversampling ratio 
        d = dynamic range of vertical axis
	pltres:  0=>plot both, 1=>plot sig only, 2=>plot res only
    """ 

    try:
        datvar = np.loadtxt( datfile )
    except:
        print( 'Cannot open' + datfile )
        return()

    # column 0 is frequency
    if Fs < 10e6:
        fa = datvar[:,0] * Fs/1000  # make x axis scale kHz
    else:
        fa = datvar[:,1] * Fs/1e6  # make x axis scale MHz 

    # add a tiny non-zero component to datvar cols 1, 2 to 
    # suppress the log(0) errors
    e = 1.0e-99*np.ones(datvar[:,1].shape)
    datvar[:,1] += e
    datvar[:,2] += e

    # column 1 is psd of signal+noise
    dBSxx = 10*np.log10( 2.0*(delta**2/Fs)*datvar[:,1] )  # scale to V^2/Hz

    # column 2 is psd of noise only 
    dBSrr = 10*np.log10( 2.0*(delta**2/Fs)*datvar[:,2] )  # scale to V^2/Hz 

    fig, ax = plt.subplots()

    # select what to plot -- 0: S+N, 1: S only, 2: N only
    if pltres==0:
        ax.plot( fa, dBSxx, label='S+N', color='blue' )
        ax.plot( fa, dBSrr, label='N', color='red' )
    elif pltres==1:
        ax.plot( fa, dBSxx, label='S+N', color='blue' )
    elif pltres==2:
        ax.plot( fa, dBSrr, label='N', color='red' )
    else:
        print( 'pltres == 0, 1, or 2 only' )
        return

    # find max over band of interest for y axis scaling
    dBSxx_max = np.amax(dBSxx[0:int(len(dBSxx)/2/OSR)])
    dBSxx_max = 10*np.ceil(dBSxx_max/10)

    if Fs < 10e6:
        ax.axis( [0, 0.5*(Fs/1000)/OSR, -d+dBSxx_max, dBSxx_max] )
        ax.set_xlabel('Frequency (kHz)')
    else:
        ax.axis( [ 0, 0.5*Fs/1e6/OSR, -d+dBSxx_max, dBSxx_max] )
        ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('PSD (dBV/Hz)')
    ax.set_title( datfile )
    ax.legend()
    plt.grid(visible=True)
    # plt.show()
    print("Saving to plot_dsm_psd_linear_freq.png")
    plt.savefig("plot_dsm_psd_linear_freq.png")
    return
            
 
def plot_dsm_psd_log_freq( datfile, delta, Fs, d, pltres ):
    """ 
    triton_dsm.plot_dsm_psd_log_freq
	Plot psd in dB; scale based on spectrum; log x axis
        Provides 2x scaling to give V^2/Hz for + freq.

	Usage:
	plot_dsm_psd_log_freq( datvar, delta, Fs, d, pltres )

	Where:
        datfile = name of file containing output from dsm.c
        datvar(:,1) = normalized frequency vector (cycles/sample)
	datvar(:,2) = psd (codes^2/(cycles/sample))
	datvar(:,3) = res psd (codes^2/(cycles/sample))
 	delta = lsb size (V)
 	Fs = sample rate (Hz)
        d = dynamic range of vertical axis
 	pltres:  0=>plot both, 1=>plot sig only, 2=>plot res onlya

    """ 
    try:
        datvar = np.loadtxt( datfile )
    except:
        print( 'Cannot open' + datfile )
        return()

    # column 0 is frequency
    datlen_by_2 = int(len(datvar[:,0])/2)  # take positive frequencies only

    if Fs < 10e6:
        fa = datvar[0:datlen_by_2,0] * Fs/1000  # make x axis scale kHz
    else:
        fa = datvar[0:datlen_by_2,1] * Fs/1e6  # make x axis scale MHz 

    # add a tiny non-zero component to datvar cols 1, 2 to 
    # suppress the log(0) errors
    e = 1.0e-99*np.ones(datvar[:,1].shape)
    datvar[:,1] += e
    datvar[:,2] += e

    # column 1 is psd of signal+noise
    dBSxx = 10*np.log10( 2.0*(delta**2/Fs)*datvar[0:datlen_by_2,1] )  # scale to V^2/Hz

    # column 2 is psd of noise only 
    dBSrr = 10*np.log10( 2.0*(delta**2/Fs)*datvar[0:datlen_by_2,2] )  # scale to V^2/Hz 

    fig, ax = plt.subplots()

    # select what to plot -- 0: S+N, 1: S only, 2: N only
    if pltres==0:
        ax.semilogx( fa, dBSxx, label='S+N', color='blue' )
        ax.semilogx( fa, dBSrr, label='N', color='red' )
    elif pltres==1:
        ax.semilogx( fa, dBSxx, label='S+N', color='blue' )
    elif pltres==2:
        ax.semilogx( fa, dBSrr, label='N', color='red' )
    else:
        print( 'pltres == 0, 1, or 2 only' )
        return

    # find max over band of interest for y axis scaling
    dBSxx_max = np.amax(dBSxx[0:datlen_by_2])
    dBSxx_max = 10*np.ceil(dBSxx_max/10)

    if Fs < 10e6:
        ax.axis( [0.5e-5*(Fs/1000), 0.5*(Fs/1000), -d+dBSxx_max, dBSxx_max] )
        ax.set_xlabel('Frequency (kHz)')
    else:
        ax.axis( [ 0.5e-5*(Fs/1e6), 0.5*(Fs/1e6), -d+dBSxx_max, dBSxx_max] )
        ax.set_xlabel('Frequency (MHz)')
    
    ax.set_ylabel('PSD (dBV/Hz)')
    ax.set_title( datfile )
    ax.legend()
    plt.grid(visible=True)
    plt.show()
    print("Saving to plot_dsm_psd_log_freq.png")
    plt.savefig("plot_dsm_psd_log_freq.png")
    return               
                    
                        
                            
                                
