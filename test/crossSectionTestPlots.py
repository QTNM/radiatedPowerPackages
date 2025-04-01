#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT


def parse_arguments():
    parser = argparse.ArgumentParser(description="""Cross section plotter""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input",
                        default="crosssections.root",
                        type=str,
                        help="Input ROOT file with TGraphs",
                        nargs="?")
    parser.add_argument("--output",
                        default="crosssections_H.png",
                        type=str,
                        help="Output PNG file",
                        nargs="?")
    return parser.parse_args()


def make_plots(fname_in="crosssections.root", fname_out="crosssections_H.png"):
    input_file = ROOT.TFile(fname_in)
    graphHElastic = input_file.Get("grHElastic")
    graphHIonisation = input_file.Get("grH")
    graphHShah = input_file.Get("grHShah")

    plt.rc('font', family='serif')

    # Get the x and y values from the TGraphs
    x_elastic = np.array(graphHElastic.GetX())
    y_elastic = np.array(graphHElastic.GetY())
    x_ionisation = np.array(graphHIonisation.GetX())
    y_ionisation = np.array(graphHIonisation.GetY())
    x_shah = np.array(graphHShah.GetX())
    y_shah = np.array(graphHShah.GetY())

    plt.plot(x_elastic, y_elastic, label="Elastic", color="blue")
    plt.plot(x_ionisation, y_ionisation, label="Ionisation", color="red")
    plt.plot(x_shah, y_shah, '.', label="Shah (1987)", color="green")
    plt.xlabel("Energy [eV]")
    plt.ylabel(r"Cross section [$10^{-20}$ m$^2$]")
    plt.xscale("log")
    plt.legend()
    plt.savefig(fname_out)


if __name__ == "__main__":
    options = parse_arguments()
    make_plots(fname_in=options.input, fname_out=options.output)
