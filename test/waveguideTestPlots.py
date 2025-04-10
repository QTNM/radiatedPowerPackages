#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py


def parse_arguments():
    parser = argparse.ArgumentParser(description="""Cross section plotter""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input",
                        default="wg_data.h5",
                        type=str,
                        help="Input HDF5 file with data",
                        nargs="?")
    parser.add_argument("--output",
                        default="collectedpowers.pdf",
                        type=str,
                        help="Output file string",
                        nargs="?")
    return parser.parse_args()


def make_plots(fname_in="wg_data.h5", string_out="collectedpowers.pdf"):
    # Set LaTeX font for matplotlib
    plt.rc('font', family='serif')

    with h5py.File(fname_in, 'r') as f:
        # Assuming the data is stored in a dataset named 'data'
        circWgPos = f['Circular']['collectedPower'][:, 0]
        circWgPowerFrac = f['Circular']['collectedPower'][:, 1]
        wr42Pos = f['WR42']['collectedPower'][:, 0]
        wr42PowerFrac = f['WR42']['collectedPower'][:, 1]

    plt.figure(figsize=(12, 5))
    plt.subplot(121)
    plt.plot(circWgPos * 1e3, circWgPowerFrac)
    plt.xlabel("Position [mm]")
    plt.ylabel("Collected power fraction")
    plt.title("Power fraction in circular waveguide")
    plt.xlim(-5, 5)
    plt.ylim(0, 1.1 * np.max(circWgPowerFrac))

    plt.subplot(122)
    plt.plot(wr42Pos * 1e3, wr42PowerFrac)
    plt.xlabel("x [mm]")
    plt.ylabel("Colled power fraction")
    plt.title("Power fraction in WR42 waveguide")
    plt.xlim(-5, 5)
    plt.ylim(0, 1.1 * np.max(wr42PowerFrac))
    plt.savefig(string_out)


if __name__ == "__main__":
    options = parse_arguments()
    make_plots(fname_in=options.input, string_out=options.output)
