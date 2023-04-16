# Benchmark for correlation using xor for low quantification signals (typically 1 bit)


## Overview

Correlation of a received signal with a reference signal can be computed using np.convolve(sig, ref, "valid") on signed integer values.
Faster implementation can be searched for unsigned integer values with only a few bits for both received and reference signals.

In the special case with only 1 (signed) bit for both, products can be implemented using xor.

On 32-bit machines, 32 multiplications can be performed in a single cycle but accumulations may require many more cycles if popcount instruction is not available.

### Introductory results

FFS

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


## Implemented algorithms

1. `corrxor.alg.goertzel`: Normal Goertzel algorithm.

## Algorithm verification and benchmark

* Run all benchmark cases and plot result

  ```bash
  $ python3 bench_and_test.py
  ```

## Reference

