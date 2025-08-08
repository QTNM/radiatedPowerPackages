/// GraphUtils.cxx - ROOT TGraph utility functions
#include "utilities/ROOTUtils/GraphUtils.h"

#include <iostream>
#include <cfloat>

#include "utilities/BasicCore/Constants.h"
#include "utilities/ROOTUtils/HistogramUtils.h"
#include "utilities/ROOTUtils/FFTAnalysis.h"
#include <cmath>
#include "TH1D.h"

using namespace rad;

void rad::setGraphAttr(TGraph *gr) {
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}

void rad::setGraphAttr(std::unique_ptr<TGraph> &gr) {
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetLineWidth(2);
}

void rad::setGraphAttr(TGraph &gr) {
  gr.GetXaxis()->SetTitleSize(0.05);
  gr.GetYaxis()->SetTitleSize(0.05);
  gr.GetXaxis()->SetLabelSize(0.05);
  gr.GetYaxis()->SetLabelSize(0.05);
  gr.SetLineWidth(2);
  gr.SetLineColor(kBlack);
  gr.SetMarkerStyle(20);
  gr.SetMarkerColor(kBlack);
}

TGraph *rad::DownmixInPhase(TGraph *grInput, const double freq) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(
        i, time,
        grInput->GetPointY(i) * std::cos(2 * PI * freq * time));
  }
  return grOut;
}

TGraph *rad::DownmixQuadrature(TGraph *grInput, const double freq) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);
  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    grOut->SetPoint(
        i, time,
        grInput->GetPointY(i) * std::sin(2 * PI * freq * time));
  }
  return grOut;
}

void rad::ScaleGraph(TGraph *grInput, const double scale) {
  for (int i = 0; i < grInput->GetN(); i++) {
    grInput->SetPointY(i, grInput->GetPointY(i) * scale);
  }
}

TGraph *rad::SumGraphs(std::vector<TGraph *> grInput) {
  TGraph *grOut = new TGraph();
  setGraphAttr(grOut);

  // First of all check that our graphs have the same x (time) spacing
  double testSpacing = grInput[0]->GetPointX(1) - grInput[0]->GetPointX(0);
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    double thisSpacing =
        grInput[iGr]->GetPointX(1) - grInput[iGr]->GetPointX(0);
    // Return empty graph if not
    if ((thisSpacing - testSpacing) / testSpacing > 1e-10) {
      std::cout << "Graphs do not have equivalent time spacing! -- returning "
                   "empty graph."
                << std::endl;
      return grOut;
    }
  }

  // The time series may be of different lengths
  // Need to determine the overlapping times between all the graphs
  double latestStart = -DBL_MAX;
  double earliestEnd = DBL_MAX;
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    if (grInput[iGr]->GetPointX(0) > latestStart)
      latestStart = grInput[iGr]->GetPointX(0);

    if (grInput[iGr]->GetPointX(grInput[iGr]->GetN() - 1) < earliestEnd)
      earliestEnd = grInput[iGr]->GetPointX(grInput[iGr]->GetN() - 1);
  }

  std::vector<int> startIndices;
  std::vector<int> endIndices;
  // Find the indices where these start and end times are satisfied
  for (int iGr = 0; iGr < grInput.size(); iGr++) {
    // Loop through points from start
    for (int iPnt = 0; iPnt < grInput[iGr]->GetN(); iPnt++) {
      if (grInput[iGr]->GetPointX(iPnt) == latestStart) {
        startIndices.push_back(iPnt);
        break;
      }
    }
    // Now loop through the points from the end to get the end
    for (int iPnt = grInput[iGr]->GetN() - 1; iPnt >= 0; iPnt--) {
      if (grInput[iGr]->GetPointX(iPnt) == earliestEnd) {
        endIndices.push_back(iPnt);
        break;
      }
    }  // Loop through points
  }

  // Now sum the graphs between the appropriate ranges
  int nPointsToSum = endIndices[0] - startIndices[0];
  for (int iPnt = 0; iPnt < nPointsToSum; iPnt++) {
    double xVal = grInput[0]->GetPointX(startIndices[0] + iPnt);
    double yVal = 0;
    for (int iGr = 0; iGr < grInput.size(); iGr++) {
      yVal += grInput[iGr]->GetPointY(startIndices[iGr] + iPnt);
    }
    grOut->SetPoint(grOut->GetN(), xVal, yVal);
  }

  return grOut;
}

TGraph *rad::SampleWaveform(TGraph *grInput, const double sRate) {
  TGraph *grOut = new TGraph();
  const double sampleSpacing = 1.0 / sRate;
  double sampleTime = grInput->GetPointX(0);

  for (int i = 0; i < grInput->GetN(); i++) {
    double time = grInput->GetPointX(i);
    if (time < sampleTime)
      continue;
    else if (i == 0) {
      double calcV = grInput->GetPointY(0);
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    } else {
      // Sample the distribution using linear interpolation
      double calcV = grInput->GetPointY(i - 1) +
                     (sampleTime - grInput->GetPointX(i - 1)) *
                         (grInput->GetPointY(i) - grInput->GetPointY(i - 1)) /
                         (time - grInput->GetPointX(i - 1));
      grOut->SetPoint(grOut->GetN(), sampleTime, calcV);
      sampleTime += sampleSpacing;
    }
  }

  return grOut;
}

TH1D *rad::GraphToHistogram(TGraph *grInput) {
  double binWidth = grInput->GetPointX(1) - grInput->GetPointX(0);
  double firstPoint = grInput->GetPointX(0);
  double lastPoint = grInput->GetPointX(grInput->GetN() - 1);
  TH1D *h = new TH1D("h", "", grInput->GetN(), firstPoint - binWidth / 2,
                     lastPoint + binWidth / 2);
  SetHistAttr(h);
  for (int i = 0; i < grInput->GetN(); i++) {
    h->SetBinContent(i + 1, grInput->GetPointY(i));
  }

  return h;
}

TGraph *rad::SignalProcessGraph(TGraph *grInput, const double downmixFreq,
                                const double sampleRate) {
  TGraph *grDM = DownmixInPhase(grInput, downmixFreq);
  TGraph *grS1 = SampleWaveform(grDM, 10 * sampleRate);
  delete grDM;
  TGraph *grF = BandPassFilter(grS1, 0.0, sampleRate / 2.0);
  delete grS1;
  TGraph *grS2 = SampleWaveform(grF, sampleRate);
  delete grF;

  return grS2;
}

double rad::SumPower(const TGraph *gr, int firstBin, int lastBin) {
  double integral{0};
  double deltaF = gr->GetPointX(1) - gr->GetPointX(0);
  if (firstBin < 0) firstBin = 0;
  if (lastBin < 0) lastBin = gr->GetN() - 1;
  for (int i = firstBin; i <= lastBin; i++) {
    integral += gr->GetPointY(i);
  }
  return integral;
}

double rad::SumVoltageSquared(const TGraph *gr, int firstBin, int lastBin) {
  double integral{0};
  if (firstBin < 0) firstBin = 0;
  if (lastBin < 0) lastBin = gr->GetN() - 1;
  for (int i = firstBin; i <= lastBin; i++) {
    integral += gr->GetPointY(i) * gr->GetPointY(i);
  }
  return integral;
}