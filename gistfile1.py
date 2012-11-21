import math

def goertzel(samples, sample_rate, f_start, f_end):
    """
    Implementation of the Goertzel algorithm, useful for calculating individual
    terms of a discrete Fourier transform.
    """
    window_size = len(samples)
    f_step = sample_rate / float(window_size)
    # Calculate which DFT bins we'll have to compute
    k_start = int(math.ceil(f_start / f_step))
    k_end = int(math.floor(f_end / f_step))
    
    if k_end > window_size - 1: raise ValueError('frequency out of range %s' % k_end)

    # For all the bins between `f_start` and `f_end`, calculate the DFT
    # term
    n_range = range(0, window_size)
    freqs = []
    results = []
    for k in range(k_start, k_end + 1):
        # Bin frequency and coefficients for the computation
        f = k * f_step
        w_real = 2.0 * math.cos(2.0 * math.pi * f)
        w_imag = math.sin(2.0 * math.pi * f)

        # Doing the calculation on the whole sample
        d1, d2 = 0.0, 0.0
        for n in n_range:
            y  = samples[n] + w_real * d1 - d2
            d2, d1 = d1, y

        # Storing results `(real part, imag part, power)`
        results.append((
            0.5 * w_real * d1 - d2, w_imag * d1,
            d2**2 + d1**2 - 2 * w_real * d1 * d2)
        )
        freqs.append(f)
    return freqs, results

if __name__ == '__main__':
    # quick test
    import numpy as np
    import pylab

    t = np.linspace(0, 1, 44100)
    sine_wave = np.sin(2*np.pi*441*t)[:1024]
    freqs, results = goertzel(sine_wave, 44100, 0, 22049)
    print np.array(results)
    pylab.plot(freqs, np.array(results)[:,2])
    pylab.show()