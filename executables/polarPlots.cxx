// polarPlots.cxx



#include "TFile.h"
#include "TH2.h"
#include "TEllipse.h"
#include "TTree.h"

using namespace rad;

int main()
{
  const int nBinsX = 200;
  const int nBinsY = 200;
  const int nBinsZ = 200;

  TH2D* h2FieldXY = new TH2D("h2FieldXY", "|S|; x [m]; y [m]; |S| [W m^{-2}]", nBinsX, -0.05, 0.05, nBinsY, -0.05, 0.05);
  SetHistAttr(h2FieldXY);
  TH2D* h2FieldXZ = new TH2D("h2FieldXZ", "|S|; x [m]; z [m]; |S| [W m^{-2}]", nBinsX, -0.05, 0.05, nBinsZ, -0.05, 0.05);
  SetHistAttr(h2FieldXZ);

  TFile* fin = new TFile("/home/sjones/work/qtnm/trajectories/90DegOnAxis.root", "read");
  TTree* tree = (TTree*)fin->Get("tree");
  double time;
  double xPos, yPos, zPos;
  double xVel, yVel, zVel;
  double xAcc, yAcc, zAcc;
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("xPos", &xPos);
  tree->SetBranchAddress("yPos", &yPos);
  tree->SetBranchAddress("zPos", &zPos);
  tree->SetBranchAddress("xVel", &xVel);
  tree->SetBranchAddress("yVel", &yVel);
  tree->SetBranchAddress("zVel", &zVel);
  tree->SetBranchAddress("xAcc", &xAcc);
  tree->SetBranchAddress("yAcc", &yAcc);
  tree->SetBranchAddress("zAcc", &zAcc);

  tree->GetEntry(1);
  ROOT::Math::XYZPoint pos(xPos, yPos, zPos);
  ROOT::Math::XYZVector vel(xVel, yVel, zVel);
  ROOT::Math::XYZVector acc(xAcc, yAcc, zAcc);

  std::cout<<"x, y, z = "<<xPos<<", "<<yPos<<", "<<zPos<<std::endl;
  std::cout<<"xVel, yVel, zVel = "<<xVel<<", "<<yVel<<", "<<zVel<<std::endl;
  std::cout<<"xAcc, yAcc, zAcc = "<<xAcc<<", "<<yAcc<<", "<<zAcc<<std::endl;
  
  for (int x = 1; x <= nBinsX; x++) {
    for (int y = 1; y <= nBinsY; y++) {
      ROOT::Math::XYZPoint detPos(h2FieldXY->GetXaxis()->GetBinCenter(x),
		      h2FieldXY->GetYaxis()->GetBinCenter(y), 0.0);
      ROOT::Math::XYZVector poyntingVec = CalcPoyntingVec(detPos, pos, vel, acc);
      h2FieldXY->SetBinContent(x, y, sqrt(poyntingVec.Mag2()));
    }
    
    for (int z = 1; z <= nBinsZ; z++) {
      ROOT::Math::XYZPoint detPos(h2FieldXZ->GetXaxis()->GetBinCenter(x), 0.0,
		      h2FieldXZ->GetYaxis()->GetBinCenter(z));
      ROOT::Math::XYZVector poyntingVec = CalcPoyntingVec(detPos, pos, vel, acc);
      h2FieldXZ->SetBinContent(x, z, sqrt(poyntingVec.Mag2()));
    }
  }

  fin->Close();
  delete fin;
  
  TFile* fout = new TFile("polarPlotsOut.root", "recreate");
  fout->cd();
  h2FieldXY->Write("h2FieldXY");
  h2FieldXZ->Write("h2FieldXZ");
  TEllipse *el = new TEllipse(0, 0, 0.05, 0.05);
  el->SetLineWidth(2);
  el->SetLineColor(kRed);
  el->Write("el");
  
  fout->Close();
  delete fout;
  
  return 0;
}
