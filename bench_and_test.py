#!/usr/bin/python3
import numpy as np
import sys

from enum import Enum
from corrxor_directory import corrxor
from pathlib import Path
from pyfftw.interfaces import numpy_fft
from time import time


def bench_goertzel(BenchType, data_len, cx=False, n_test=10000):
    rng = np.random.Generator(np.random.SFC64())
    dtype = np.float32
    in_data = rng.random((n_test, data_len), dtype)
    if cx:
        in_data = in_data+1j*rng.random((n_test, data_len), dtype)
    in_data = in_data.astype(np.float64)

    cost = np.empty(len(BenchType))
    error = np.empty(len(BenchType))

    # fft
    for name in ["fft", "fftf"]:
        try:
            etype = BenchType[name]
        except KeyError:
            continue
        is_f32 = name in ["fftf"]
        if is_f32:
            in_data2 = in_data.astype(np.float32)
            t0 = time()
            out = numpy_fft.fft(in_data2)
            cost[etype.value] = time()-t0
        else:
            t0 = time()
            out = np.fft.fft(in_data)
            cost[etype.value] = time()-t0
            out_fft = out

        if is_f32:
            error[etype.value] = (out-out_fft).std(axis=-1).max()
        else:
            error[etype.value] = np.nan

    # dft, protection against memory errors
    try:
        if data_len < 1024:
            x = np.arange(data_len)
            F = np.exp(-2j*np.pi/data_len*np.outer(x, x))
            t0 = time()
            out = in_data @ F
            cost[BenchType.dft.value] = time()-t0
        else:
            out = np.nan
            cost[BenchType.dft.value] = np.nan
        error[BenchType.dft.value] = (out-out_fft).std(axis=-1).max()
    except AttributeError:
        pass

    for etype in BenchType:
        if "_dft" not in etype.name:
            continue
        try:
            fun = getattr(corrxor, etype.name)
        except AttributeError:
            continue
        t0 = time()
        out = fun(in_data, np.nan)
        cost[etype.value] = time()-t0
        error[etype.value] = (out-out_fft).std(axis=-1).max()

    for etype in BenchType:
        name = etype.name
        if "_dft" in name or name in ["dft", "fft", "fftf"]:
            continue
        is_f32 = "goertzelf" in name
        if is_f32:
            name = name.replace("goertzelf", "goertzel")
            in_data2 = in_data.astype(np.float32)
        else:
            in_data2 = in_data
        try:
            fun = getattr(corrxor, name)
        except AttributeError:
            fun = globals()[name]
            # try:
            # except KeyError:
            #     continue
        out = np.empty_like(out_fft)
        cost_goertzel = 0
        for k in range(data_len):
            t0 = time()
            out_k = fun(in_data2, k)
            cost_goertzel += time()-t0
            out[:, k] = out_k
        cost[etype.value] = cost_goertzel/data_len
        error[etype.value] = (out-out_fft).std(axis=-1).max()

    return cost, error


def bench_range(BenchType, len_range, cx=False, n_test=10000):
    cost = np.empty((len(BenchType), len(len_range)))
    error = np.empty((len(BenchType), len(len_range)))
    for i, data_len in enumerate(len_range):
        print("data_len %d" % data_len)
        cost[:, i], error[:, i] = bench_goertzel(
            BenchType, data_len, cx=cx, n_test=n_test)

    return cost, error


def plot_bench(
        BenchType, len_range, cost, error, bench_list, media_path, plot_prefix,
        title_str):

    prefix = "/".join([str(media_path), plot_prefix])
    from matplotlib import pyplot as plt
    plt.close("all")
    fig_size = (12, 9)
    plt.figure(1, figsize=fig_size)
    for name in bench_list:
        k = BenchType[name]
        plt.plot(len_range, 10*np.log10(cost[k.value, :]), label=k.name)
    plt.legend()
    plt.ylabel("time (dBs)")
    plt.xlabel("length (samples)")
    plt.title(title_str)
    plt.savefig("%s_cost_db.png" % prefix, bbox_inches="tight")

    plt.figure(2, figsize=fig_size)
    for name in bench_list:
        k = BenchType[name]
        plt.plot(len_range, cost[k.value, :], label=k.name)
    plt.legend()
    plt.ylabel("time (s)")
    plt.xlabel("length (samples)")
    plt.title(title_str)
    plt.savefig("%s_cost.png" % prefix, bbox_inches="tight")

    plt.figure(3, figsize=fig_size)
    for name in bench_list:
        k = BenchType[name]
        if np.isnan(error[k.value, :]).all():
            plt.plot(np.arange(1), np.arange(1), label="%s vs fft" % k.name)
            continue
        error_k = error[k.value, :]
        error_k = 10*np.log10(error_k)
        plt.plot(len_range, error_k, label="%s vs fft" % k.name)
    plt.legend()
    plt.ylabel("error")
    plt.xlabel("length (samples)")
    plt.title(title_str)
    plt.savefig("%s_error.png" % prefix, bbox_inches="tight")
    # plt.show()
    plt.pause(1)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print("dummy exit")
        sys.exit()

    media_path = Path(".") / "media_tmp"
    media_path.mkdir(exist_ok=True)
    n_test = 100

    len_range = np.concatenate((
        np.arange(1, 64), np.arange(64, 1024, 64), 2**np.arange(10, 13)))

    # real-value tests---------------------------------------------------------
    bench_list = [
        "dft",
        "fft",
        "fftf",
        "goertzel",
        "goertzelf",
    ]
    BenchType = Enum("BenchType", bench_list, start=0)
    cost, error = bench_range(BenchType, len_range, n_test=n_test)

    title_str = "Faster Goertzel w/ radix"
    bench_list = [
        "fft",
        "goertzel",
    ]
    plot_bench(
        BenchType, len_range, cost, error, bench_list, media_path, "radix",
        title_str)
