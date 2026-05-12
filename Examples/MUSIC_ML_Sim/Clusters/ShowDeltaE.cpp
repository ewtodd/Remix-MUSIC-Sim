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


int ShowDeltaE(){

  int select_strip = 4;  // Select strip number to be analyzed (0-9), -1 to skip
  int devi = select_strip + 1;
  int devf = select_strip + 5;
  int dehi = select_strip;
  int dehf = 16;

  bool savehisto = false;

  auto can = new TCanvas("can1","traces",1200,900);
  
  auto hbk  = new TH2F("hbk","",500,0,60,500,0,30);
  hbk->SetStats(0);
  hbk->GetXaxis()->SetTitle("#DeltaE1 (MeV)");
  hbk->GetXaxis()->CenterTitle();
  hbk->GetYaxis()->SetTitle("#DeltaE0 (MeV)");
  hbk->GetYaxis()->CenterTitle();
  hbk->Draw("colz");

  auto htr = new TH2F("htr","",18,-0.5,17.5,500,0,5);
  htr->SetStats(0);
  htr->GetXaxis()->SetTitle("strip number");
  htr->GetXaxis()->CenterTitle();
  htr->GetYaxis()->SetTitle("Energy deposited in left+right segments (MeV)");
  htr->GetYaxis()->CenterTitle();
  htr->Draw();

  const int maxClu = 10;
  auto hde = new TH2F*[maxClu];
  for (int clu=0; clu<maxClu; clu++) {
    hde[clu] = new TH2F(Form("hde%d",clu),"",500,0,60,500,0,30);
    if (clu+1==10)
      hde[clu]->SetMarkerColor(11);
    else
      hde[clu]->SetMarkerColor(clu+1);
    hde[clu]->SetMarkerStyle(20+clu);
    hde[clu]->SetMarkerSize(3);
  }
  
  ifstream datafile;
  CSVRow row;
  TGraph* gr = new TGraph();
  //  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  for (int clu=0; clu<maxClu; clu++) {
    datafile.open(Form("data_cluster_%d%d.0.csv",clu, select_strip));
    if (datafile.is_open()) {      
      cout << Form("data_cluster_%d%d.0.csv",clu, select_strip) << endl;
      while (datafile >> row) {     
	float eid = atof(row[0].c_str());
	float cla = atof(row[19].c_str());
	float stp = atof(row[20].c_str());
	float de0 = 0;
	float de1 = 0;

	for (int i=0; i<21; i++) 
	  cout << "[" << i << "]="<< atof(row[i].c_str()) << " ";
	cout << endl;
	
	for (int stp=0; stp<17; stp++) {
	  gr->SetPoint(stp, (float)stp, atof(row[stp+1].c_str()));
	  if (stp>=devi && stp<=devf)
	    de0 += atof(row[stp+1].c_str());
	  if (stp>=dehi && stp<=dehf)
	    de1 += atof(row[stp+1].c_str());
	}

	if (de1<30) {
	  htr->SetTitle(Form("event %.0f stp %.0f", eid, stp));
	  gr->DrawClone("lp");
	  // can->Update();
	  // can->WaitPrimitive();
	}
	cout << "  " << de0 << " " << de1 << endl;
	hde[clu]->Fill(de1,de0);
	hbk->Fill(de1,de0);
      }
      datafile.close();
    }
  }

  can = new TCanvas("can2","deltaE",1200,900);
  hbk->Draw("colz");
  for (int clu=0; clu<maxClu; clu++) {
    if (hde[clu]->GetEntries()>0)
      hde[clu]->Draw("same");
  }

  // if (savehisto) {
  //   TFile* RF = new TFile("hist_traces.root", "recreate");
  //   htraces->Write();
  //   RF->Close();
  // }
    
  
  return 0;
}
