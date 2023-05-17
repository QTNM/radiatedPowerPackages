// rayleighDistributionChecks.cxx

#include "BasicFunctions/BasicFunctions.h"
#include "Antennas/HalfWaveDipole.h"
#include "SignalProcessing/Signal.h"
#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"

#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

using namespace rad;

double rayleighNorm(double *x, double *par)
{
  double xx = sqrt(x[0] * 18.5e3 * 750e6 / (0.5 * 20270));
  double retVal = RayleighPDF(xx, par[0]);
  return retVal;
}

int main(int argc, char** argv)
{
  std::cout<<"Have "<<argc<<" arguments"<<std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<std::endl;
  }
  TString outputFile = argv[1];

  const double TEff = 4.0;
  const double tAcq = 54e-6;
  const double bandwidth = 1 / tAcq; // Hz
  std::cout<<"Bandwidth "<<bandwidth<<std::endl;
  const double sampleRate = 750e6; // Hz
  const double deltaT = 1 / sampleRate;
  const int N = sampleRate / bandwidth;
  int newLength = (N/2) + 1;
  double* x1 = new double[newLength];
  double* y1 = new double[newLength];
  double tempF = 0;
  
  for (int i = 0; i < newLength; i++) {
    x1[i] = tempF;
    tempF += bandwidth;
    y1[i] = 0;
  }
  
  TGraph* grNoise = new TGraph(newLength, x1, y1);
  AddWhiteNoiseFrequencyDomainPowerNorm(grNoise, TEff);
  setGraphAttr(grNoise);
  TH1D* hNoiseVals = new TH1D("hNoiseVals", "", 100, 0, 5e-9);
  SetHistAttr(hNoiseVals);
  for (int i = 0; i < grNoise->GetN(); i++) {
    double noiseRoot = sqrt(grNoise->GetPointY(i)*bandwidth/(deltaT*0.5*N));
    hNoiseVals->Fill(noiseRoot);
  }
		    
  const double sigma = TMath::Sqrt(TMath::K() * TEff * bandwidth);
  TF1* fCDF = new TF1("fCDF", RayleighCDFFunc, 0, 4*sigma, 1);
  fCDF->SetParameter(0, sigma);
  fCDF->SetTitle("CDF: Rayleigh distribution with T = 4 K, B = 18.5 kHz; x; f(x; #sigma)");
  fCDF->GetHistogram()->GetXaxis()->SetTitle("#sqrt{k_{B}TB} [W^{0.5}]");

  TF1* fPDF = new TF1("fPDF", RayleighPDFFunc, 0, 4*sigma, 1);//rayleighNorm, 0, 25e-27, 1);
  fPDF->SetParameter(0, sigma*sqrt(0.5));
  fPDF->SetTitle("PDF: Rayleigh distribution with T = 4 K, B = 18.5 kHz; x; f(x; #sigma)");
  fPDF->GetHistogram()->GetXaxis()->SetTitle("#sqrt{k_{B}TB} [W^{0.5}]");

  TH1D* hRayleighSq = new TH1D("hRayleighSq", "", 200, 0, 1e-17);
  for (int ii = 0; ii < 500000; ii++) {
    hRayleighSq->Fill(pow(fPDF->GetRandom(), 2));
  }

  TVector3 antennaPoint1(0.05, 0, 0); 
  TVector3 antennaDirZ1(0, 1, 0);
  TVector3 antennaDirX1(1, 0, 0);
  HalfWaveDipole* antenna1 = new HalfWaveDipole(antennaPoint1, antennaDirX1, antennaDirZ1, 27.01e9);

  const double loadResistance = 70; // Ohms
  const double downmixFreq = 26.75e9; // Hertz 
  
  InducedVoltage iv1("/home/sjones/work/qtnm/trajectories/electronTraj1ms90Deg.root", antenna1);
  LocalOscillator myLO(downmixFreq * 2 * TMath::Pi());
  GaussianNoise noise1(TEff, loadResistance);
  Signal sig1(iv1, myLO, sampleRate, {noise1}, tAcq);
  TGraph* grVI = sig1.GetVITimeDomain();
  TGraph* grVQ = sig1.GetVQTimeDomain();
  std::vector<TGraph*> vecV = {grVI, grVQ};
  TGraph* grVIPgram = MakePowerSpectrumPeriodogram(grVI);
  TGraph* grVQPgram = MakePowerSpectrumPeriodogram(grVQ);
  grVIPgram->GetYaxis()->SetTitle("Power [W]");
  grVQPgram->GetYaxis()->SetTitle("Power [W]");
  ScaleGraph(grVIPgram, 1/loadResistance);
  ScaleGraph(grVQPgram, 1/loadResistance);
  double totalPowerVI = FFTtools::sumPower(grVIPgram) * 1e15;
  double totalPowerVQ = FFTtools::sumPower(grVQPgram) * 1e15;  
  grVIPgram->SetTitle(Form("Total power %.4f fW", totalPowerVI));
  grVQPgram->SetTitle(Form("Total power %.4f fW", totalPowerVQ));
  TGraph* grSum = SumGraphs(vecV);
  TGraph* grSumPgram = MakePowerSpectrumPeriodogram(grSum);
  double totalPowerSum = FFTtools::sumPower(grSumPgram) * 1e15;
  ScaleGraph(grSumPgram, 1/loadResistance);
  grSumPgram->SetTitle(Form("Total power %.4f fW", totalPowerSum));

  auto hNoiseSqrtPower{std::make_unique<TH1D>("hNoiseSqrtPower", "#sqrt{P} in 18.5 kHz bins; #sqrt{P} [W^{0.5}]; N_{bins}", 70, 0, 5e-9)};
  SetHistAttr(hNoiseSqrtPower);
  auto hNoisePower{std::make_unique<TH1D>("hNoisePower", "P in 18.5 kHz bins; P [W]; N_{bins}", 70, 0, 2.5e-17)};
  SetHistAttr(hNoisePower);
  
  for (int i = 0; i < grSumPgram->GetN(); i++) {
    hNoiseSqrtPower->Fill(sqrt(grVIPgram->GetPointY(i)));
    hNoisePower->Fill(grVIPgram->GetPointY(i));
  }
  
  auto fout{std::make_unique<TFile>(outputFile, "RECREATE")};
  fout->cd();
  grVIPgram->Write("grVIPgram");
  grVQPgram->Write("grVQPgram");
  grSumPgram->Write("grSumPgram");
  fCDF->Write();
  fPDF->Write();
  grNoise->Write("grNoise");
  hRayleighSq->Write();  
  hNoiseVals->Write();
  hNoisePower->Write();
  hNoiseSqrtPower->Write();
  
  fout->Close();  
  return 0;
}
