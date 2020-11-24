import cochlea as co
import thorns as th
import numpy as np
import matplotlib.pyplot as plt

def sig_to_train(sig, fs, db=80.0):
    sig = co.set_dbspl(sig, db)
    spikes = co.run_zilany2014(sig, fs, anf_num=(10, 7, 3), cf = (
        125.0, 6400.0, 29), species='human', seed=1)
    return spikes 


def static_complex(f_vals, dur, fs, magnitudes=None, phase_vals=None,
        truncate=False):
    '''
    '''
    ## Deal with parameters
    if phase_vals==None:
        phase_vals = np.zeros_like(f_vals)
    assert len(phase_vals) == len(f_vals), \
            "Incompatible number of frequencies and phase offsets!"
    if magnitudes==None:
        magnitudes = np.ones_like(f_vals)
    assert len(magnitudes)==len(f_vals), \
            "Incompatible number of frequencies and magnitudes!"
    dt = 1./fs
    if truncate:
        T = 1/f_vals[0]
        num_periods = np.floor(dur/T)
        dur = T*num_periods
    t = np.arange(0, dur, dt)
    ## Generate Signal
    sig = np.zeros_like(t)
    for p in range(1, len(f_vals)+1):
        sig += magnitudes[p-1]*np.cos(2*np.pi*f_vals[p-1]*p*t + phase_vals[p-1])
    return sig

def harmonic_complex(f_0, dur, fs, num_h):
    return static_complex([f_0*k for k in range(1, num_h+1)],
                          dur,
                          fs,
                          truncate=True)
if __name__=="__main__":
    dur = 0.1
    fs = 100e3
    f0 = 155.0
    sig = harmonic_complex(f0, dur, fs, 6)
    spikes = sig_to_train(sig, fs)
    spikes.to_pickle(r"/home/dahlbom/research/rinzelmodel/spikes.pkl")
    anf_acc = th.accumulate(spikes, keep=['cf', 'duration'])
    anf_acc.sort_values('cf', ascending=False, inplace=True)
    # isih, bins = th.isih(spikes, 5e-05)

    fig, ax = plt.subplots(3,1)
    th.plot_signal(
        signal=sig,
        fs=fs,
        ax=ax[0]
    )
    th.plot_neurogram(
        anf_acc,
        fs,
        ax=ax[1]
    )
    th.plot_isih(
        spikes,
        5e-05,
        ax=ax[2]
    )
    plt.show()
