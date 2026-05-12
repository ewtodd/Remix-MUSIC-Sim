#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


// Reference
// http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
class CSVRow
{
public:
  std::string const& operator[](std::size_t index) const
  {
    return m_data[index];
  }
  std::size_t size() const
  {
    return m_data.size();
  }
  void readNextRow(std::istream& str)
  {
    std::string         line;
    std::getline(str, line);

    std::stringstream   lineStream(line);
    std::string         cell;

    m_data.clear();
    while(std::getline(lineStream, cell, ','))
      {
	m_data.push_back(cell);
      }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
      {
	// If there was a trailing comma then add an empty element.
	m_data.push_back("");
      }
  }
private:
  std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
  data.readNextRow(str);
  return str;
}   


using namespace std;


int ShowTraces(){

  int cluster = 1;
  int select_strip = 4;  // Select strip number to be analyzed (0-9), -1 to skip
  bool savehisto = false;
  
  TH1F* hstp = new TH1F("hstp","",18,-0.5,17.5);
  TFile* fileStp = 0;
  TH2F* hbk = 0;
  ifstream datafile;

  hbk = new TH2F("hbk","",18,-0.5,17.5,500,0,5);
  hbk->SetStats(0);
  hbk->GetXaxis()->SetTitle("strip number");
  hbk->GetXaxis()->CenterTitle();
  hbk->GetYaxis()->SetTitle("Energy deposited in left+right segments (MeV)");
  hbk->GetYaxis()->CenterTitle();
  TH2F* htraces = new TH2F("htraces","traces",18,-0.5,17.5,500,0,5);
  htraces->SetStats(0);
  
  datafile.open(Form("data_cluster_%d%d.0.csv",cluster, select_strip));
  fileStp = new TFile(Form("../Results_Oct/alpha_p_stp%d.root", select_strip));
  if (fileStp)
    hbk->Add((TH2F*)fileStp->Get("htraces"));


  TGraph* gr = new TGraph();
  //  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);

  TCanvas* can = new TCanvas("can","traces",1200,900);
  hbk->Draw("colz");
  
  CSVRow row;
  while (datafile >> row) {
    float eid = atof(row[0].c_str());
    float cla = atof(row[19].c_str());
    float stp = atof(row[20].c_str());
    hstp->Fill(stp);
    if (select_strip>-1 && select_strip==stp) {
      for (auto i=0; i<18; i++) {
	gr->SetPoint(i, (float)i, atof(row[i+1].c_str()));
	htraces->Fill((float)i,  atof(row[i+1].c_str()));
      }
      hbk->SetTitle(Form("clu%d event %.0f stp %.0f class %.0f", cluster, eid, stp, cla));
      gr->Draw("lp same");
      can->Update();
      can->WaitPrimitive();
    }
  }
  
  hstp->Draw();

  if (savehisto) {
    TFile* RF = new TFile("hist_traces.root", "recreate");
    htraces->Write();
    RF->Close();
  }
    
  
  return 0;
}
