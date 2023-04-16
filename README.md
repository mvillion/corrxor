# Benchmark for correlation using xor for low quantification signals (typically 1 bit)


## Overview

Correlatation of a received signal with a reference signal can be computed using np.convolve(sig, ref, "valid") on signed integer values.
Faster implementation can be searched for unsigned integer values with only a few bits for both received and reference signals.

In the special case with only 1 (signed) bit for both, products can be implemented using xor.

On 32-bit machines, 32 multiplications can be performed in a single cycles but accumulations can be longer if popcount instruction is not available.

### Introductory results

Goertzel is indeed faster than FFT.

![Alt text](media/intro_cost_db.png?raw=true "Goertzel vs DFT vs FFT (cost)")

Cost is represented in measured time with a logarithmic scale.


## Environments
### Machine
* OS: Ubuntu 22.04
* CPU: AMD A4-5000 APU with Radeon(TM) HD Graphics
* RAM: 16.00 GB Dual-Channel DDR3 @ 667 MT/s

### Version of Python and packages
* Python: Python 3.10.7
* Numpy: 1.23.4

## Installation

* To uninstall this package:

  ```bash
  $ ./setup build_ext --inplace --cpu-baseline="avx"
  ```

## Usage
* Evaluate the power of a single DFT term by `goertzel`

  ```python
  import gofft
  import numpy as np

  fs = 1000   # sampling frequency
  ft = 60     # target frequency to be evaluated (60 Hz)
  dur = 2     # duration of signal
  num = fs*dur  # sampling points
  t = np.linspace(0, dur, num)  # time series
  data = np.sin(2*np.pi*ft*t)   # signal to be evaluated (60 Hz)

  mag = gofft.alg.goertzel(data, fs, ft, fs)
  print(mag)  # 0.4969141358692001
  ```

## Implemented algorithms

1. `corrxor.alg.goertzel`: Normal Goertzel algorithm.

## Algorithm verification and benchmark

* Run all benchmark cases and plot result

  ```bash
  $ python3 bench_and_test.py
  ```

## Reference
[wikipedia - Goertzel](https://en.wikipedia.org/wiki/Goertzel_algorithm)
[stackoverflow - Implementation of Goertzel algorithm in C](http://stackoverflow.com/questions/11579367)

[STFT]: https://en.wikipedia.org/wiki/Short-time_Fourier_transform
[launch_on_binder]: https://mybinder.org/v2/gh/NaleRaphael/goertzel-fft/master?filepath=doc%2Fipynb%2Fdemo_simple_example.ipynb
