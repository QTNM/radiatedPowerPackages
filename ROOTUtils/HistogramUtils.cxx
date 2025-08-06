/// HistogramUtils.cxx - ROOT histogram utility functions
#include "ROOTUtils/HistogramUtils.h"

void rad::SetHistAttr(TH1 *h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH1 &h) {
  h.GetXaxis()->SetTitleSize(0.05);
  h.GetYaxis()->SetTitleSize(0.05);
  h.GetXaxis()->SetLabelSize(0.05);
  h.GetYaxis()->SetLabelSize(0.05);
  h.SetLineWidth(2);
}

void rad::SetHistAttr(std::unique_ptr<TH1D> &h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->SetLineWidth(2);
}

void rad::SetHistAttr(TH2 *h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

void rad::SetHistAttr(std::unique_ptr<TH2> &h) {
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetZaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}