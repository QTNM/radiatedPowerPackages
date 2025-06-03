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
                        default="crosssections",
                        type=str,
                        help="Output start string",
                        nargs="?")
    return parser.parse_args()


def make_plots(species="H",
               fname_in="crosssections.root", string_out="crosssections"):
    input_file = ROOT.TFile(fname_in)
    graphElastic = input_file.Get("gr"+species+"Elastic")
    graphIonisation = input_file.Get("gr"+species)
    graphShah = input_file.Get("gr"+species+"Shah")

    plt.rc('font', family='serif')

    # Get the x and y values from the TGraphs
    x_elastic = np.array(graphElastic.GetX())
    y_elastic = np.array(graphElastic.GetY())
    x_ionisation = np.array(graphIonisation.GetX())
    y_ionisation = np.array(graphIonisation.GetY())
    x_shah = np.array(graphShah.GetX())
    y_shah = np.array(graphShah.GetY())

    plt.figure(figsize=(8, 6))
    plt.plot(x_elastic, y_elastic, label="Elastic", color="blue")
    plt.plot(x_ionisation, y_ionisation, label="Ionisation", color="red")
    plt.plot(x_shah, y_shah, '.', label="Shah data", color="green")
    plt.xlabel("Energy [eV]")
    plt.ylabel(r"Cross section [$10^{-20}$ m$^2$]")
    plt.title("Cross sections for "+species)
    plt.xscale("log")
    plt.ylim(0.0, 0.7)
    plt.xlim(10.0, 20e3)
    plt.legend()
    fname_out = string_out + "_" + species + ".pdf"
    print('Creating output plot file:', fname_out)
    plt.savefig(fname_out)


if __name__ == "__main__":
    options = parse_arguments()
    make_plots("H", fname_in=options.input, string_out=options.output)
    make_plots("He", fname_in=options.input, string_out=options.output)
