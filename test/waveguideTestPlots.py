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
    return parser.parse_args()


def make_plots(fname_in="wg_data.h5"):
    # Set LaTeX font for matplotlib
    plt.rc('font', family='serif')

    with h5py.File(fname_in, 'r') as f:
        # Get the collected power data. Factor of 4 accounts for downmixing
        circWgPos = f['Circular']['collectedPower'][:, 0]
        circWgPowerFrac = f['Circular']['collectedPower'][:, 1] * 4.0
        wr42Pos = f['WR42']['collectedPower'][:, 0]
        wr42PowerFrac = f['WR42']['collectedPower'][:, 1] * 4.0
        # Get the field data
        circWgField_P = f['FieldData']['circFieldDataTE11_P'][:]
        circWgField_M = f['FieldData']['circFieldDataTE11_M'][:]
        wr42Field = f['FieldData']['rectFieldDataTE10'][:]

    plt.figure(figsize=(12, 5))
    plt.subplot(121)
    plt.plot(circWgPos * 1e3, circWgPowerFrac)
    plt.xlabel("Position [mm]")
    plt.ylabel("Collected power fraction")
    plt.title("Power fraction in circular waveguide")
    plt.xlim(-5, 5)
    plt.ylim(0, 1.1)

    plt.subplot(122)
    plt.plot(wr42Pos * 1e3, wr42PowerFrac)
    plt.xlabel("x [mm]")
    plt.ylabel("Colled power fraction")
    plt.title("Power fraction in WR42 waveguide")
    plt.xlim(-5, 5)
    plt.ylim(0, 1.1)
    plt.savefig("collectedpowers.pdf", bbox_inches='tight')

    # Now make a 2D colour plot of the field magnitudes
    # The first dimension is the x-axis, the second dimension is the y-axis and the value is the field magnitude
    plt.figure(figsize=(12, 5))
    plt.subplot(121)
    x = np.linspace(-5, 5, 41)
    y = np.linspace(-5, 5, 41)
    X, Y = np.meshgrid(x, y)
    plt.subplot(121)
    plt.contourf(X, Y, circWgField_P, levels=100)
    plt.title("TE11+ field in circular waveguide")
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.subplot(122)
    plt.contourf(X, Y, circWgField_M, levels=100)
    plt.title("TE11- field in circular waveguide")
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.savefig("TE11Fields.pdf", bbox_inches='tight')

    plt.figure(figsize=(6, 5))
    plt.contourf(X, Y, wr42Field, levels=100)
    plt.title("TE10 field in WR42 waveguide")
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.savefig("TE10Fields.pdf", bbox_inches='tight')


if __name__ == "__main__":
    options = parse_arguments()
    make_plots(fname_in=options.input)
