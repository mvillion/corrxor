#!/usr/bin/python3
import numpy as np
import sys

from enum import Enum
from corrxor_directory import corrxor
from pathlib import Path
from time import time


def corrxor_py(in_sig, n_sig, in_ref, n_ref):
    in_sig = np.unpackbits(in_sig.view(np.uint8), axis=-1, bitorder="little")
    in_sig = in_sig[:, :n_sig]
    in_ref = np.unpackbits(in_ref.view(np.uint8), axis=-1, bitorder="little")
    in_ref = in_ref[:, :n_ref]

    n_out = n_sig-n_ref+1
    out = np.empty(in_sig.shape[:-1]+(n_out,), dtype=np.uint32)
    for i_out in range(n_out):
        sig_k = in_sig[:, i_out:i_out+n_ref]
        out[:, i_out] = (in_ref*sig_k+(1-in_ref)*(1-sig_k)).sum(axis=-1)

    return out


def bench_corrxor(BenchType, n_ref, n_out=3, n_test=10000):
    rng = np.random.Generator(np.random.SFC64(seed=12345))
    in_ref = rng.random((n_test, n_ref), dtype=np.float32) < 0.5
    n_sig = n_ref+n_out-1
    in_sig = rng.random((n_test, n_sig), dtype=np.float32) < 0.5

    in_ref = in_ref.astype(np.uint8)
    n_pad = in_ref.shape[-1]
    n_pad = (n_pad+31)//32*32-n_pad
    in_ref = np.column_stack((in_ref, np.zeros((n_test, n_pad), np.uint8)))
    in_ref = np.packbits(in_ref, axis=-1, bitorder="little").view(np.uint32)

    in_sig = in_sig.astype(np.uint8)
    n_pad = in_sig.shape[-1]
    n_pad = (n_pad+31)//32*32-n_pad
    in_sig = np.column_stack((in_sig, np.zeros((n_test, n_pad), np.uint8)))
    in_sig = np.packbits(in_sig, axis=-1, bitorder="little").view(np.uint32)

    cost = np.empty(len(BenchType))
    error = np.empty(len(BenchType))

    # reference test
    for name in ["corrxor_py"]:
        try:
            etype = BenchType[name]
        except KeyError:
            continue
        t0 = time()
        out_ref = corrxor_py(in_sig, n_sig, in_ref, n_ref)
        cost[etype.value] = time()-t0
        error[etype.value] = np.nan

    for etype in BenchType:
        if etype.name in ["corrxor_py"]:
            continue
        try:
            fun = getattr(corrxor, etype.name)
        except AttributeError:
            continue
        t0 = time()
        out = fun(in_sig, n_sig, in_ref, n_ref)
        cost[etype.value] = time()-t0
        error[etype.value] = (out-out_ref).std(axis=-1).max()

    return cost, error


def bench_range(BenchType, len_range, n_test=10000):
    cost = np.empty((len(BenchType), len(len_range)))
    error = np.empty((len(BenchType), len(len_range)))
    for i, data_len in enumerate(len_range):
        print("data_len %d" % data_len)
        cost[:, i], error[:, i] = bench_corrxor(
            BenchType, data_len, n_test=n_test)

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
            plt.plot(np.arange(1), np.arange(1), label="%s vs ref" % k.name)
            continue
        error_k = error[k.value, :]
        # error_k = 10*np.log10(error_k)
        plt.plot(len_range, error_k, label="%s vs ref" % k.name)
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

    # len_range = np.concatenate((
    #     np.arange(1, 64), np.arange(64, 1024, 64), 2**np.arange(10, 13)))

    len_range = 2**np.arange(10, 13)

    # real-value tests---------------------------------------------------------
    bench_list = [
        "corrxor_py",
        "corrxor",
        "corrxor_popcount",
        "corrxor_nopop",
        "corrxor_popcount_3quad",
        "corrxor_popcount_8octet",
    ]
    BenchType = Enum("BenchType", bench_list, start=0)
    cost, error = bench_range(BenchType, len_range, n_test=n_test)

    title_str = "corrxor comparison"
    bench_list = [
        "corrxor_py",
        "corrxor",
        "corrxor_popcount",
        "corrxor_nopop",
        "corrxor_popcount_3quad",
        "corrxor_popcount_8octet",
    ]
    plot_bench(
        BenchType, len_range, cost, error, bench_list, media_path, "radix",
        title_str)
    from matplotlib import pyplot as plt
    plt.show()
