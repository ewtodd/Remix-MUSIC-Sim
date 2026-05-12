#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>

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


int CompareTraces(){

  int select_strip = 2;  // Select strip number to be analyzed (0-14), -1 to skip

  TH1F* hstpaa = new TH1F("hstpaa","Distribution of event classification",18,-0.5,17.5);
  hstpaa->SetStats(0);
  hstpaa->SetLineWidth(2);
  hstpaa->SetLineColor(kGray+1);
  hstpaa->GetXaxis()->SetTitle("strip");
  hstpaa->GetXaxis()->CenterTitle();
  TH1F* hstpap = new TH1F("hstpap","Distribution of event classification",18,-0.5,17.5);
  hstpap->SetStats(0);
  hstpap->SetLineWidth(2);
  hstpap->SetLineColor(kRed);
  hstpap->GetXaxis()->SetTitle("strip");
  hstpap->GetXaxis()->CenterTitle();
  
  ifstream datafile;


  
  //  datafile.open("alpha_alpha_0.csv");
  
  TH2F* hbk = new TH2F("hbk","MUSIC data",18,-0.5,17.5,500,0,5);
  hbk->SetStats(0);
  hbk->GetXaxis()->SetTitle("strip");
  hbk->GetXaxis()->CenterTitle();
  hbk->GetYaxis()->SetTitle("Energy [MeV]");
  hbk->GetYaxis()->CenterTitle();

  TMultiGraph* grap = new TMultiGraph();
  TMultiGraph* graa = new TMultiGraph();
  TGraph* gr;


  TCanvas* can = new TCanvas("can","traces",1200,900);
  can->Divide(2);
 
  
  CSVRow row;

  // Load (a,p) data
  datafile.open("alpha_p_with_index.csv");
  while (datafile >> row) {
    int eid = atoi(row[0].c_str());
    int stp = atoi(row[19].c_str());
    hstpap->Fill(stp);
    if (select_strip>-1 && select_strip==stp) {
      //      gr = new TGraph
      gr = new TGraph();
      gr->SetLineWidth(2);
      gr->SetLineColor(kRed);
      for (auto i=0; i<18; i++)
	gr->SetPoint(i, (float)i, atof(row[i+1].c_str()));
      gr->SetTitle(Form("(a,p), s%d, e%d", stp, eid));
      gr->SetName(Form("grap_s%d_e%d", stp, eid));
      grap->Add(gr);
    }
  }
  datafile.close();
  
  // Load (a,a) data
  datafile.open("alpha_alpha_with_index.csv");
  while (datafile >> row) {
    int eid = atoi(row[0].c_str());
    int stp = atoi(row[19].c_str());
    hstpaa->Fill(stp);
    if (select_strip>-1 && select_strip==stp) {
      gr = new TGraph();
      gr->SetLineWidth(2);
      gr->SetLineColor(kGray+1);
      for (auto i=0; i<18; i++)
	gr->SetPoint(i, (float)i, atof(row[i+1].c_str()));
      gr->SetTitle(Form("(a,a), s%d, e%d", stp, eid));
      gr->SetName(Form("graa_s%d_e%d", stp, eid));
      graa->Add(gr);
    }
  }
  

  can->cd(1);
  hbk->Draw();
  graa->Draw("l");
  grap->Draw("l");
  auto leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(hstpap, "(#alpha,p)", "l");
  leg->AddEntry(hstpaa, "(#alpha,#alpha)", "l");
  leg->Draw();

  can->cd(2)->SetLogy();
  hstpap->Draw();
  hstpaa->Draw("same");
  
  return 0;
}
