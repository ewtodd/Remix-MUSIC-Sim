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

  string reac = "(a,p)";  // Select reaction: (a,p), (a,a), ...
  int select_strip = 2;  // Select strip number to be analyzed (0-14), -1 to skip

  TH1F* hstp = new TH1F("hstp","",18,-0.5,17.5);
  
  ifstream datafile;

  if (reac=="(a,p)")
    datafile.open("alpha_p_0.csv");
  else if (reac=="(a,a)")
    datafile.open("alpha_alpha_0.csv");
  
  TH2F* hbk = new TH2F("hbk","",18,-0.5,17.5,500,0,5);
  hbk->SetStats(0);
  
  TGraph* gr = new TGraph();
  //  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);

  TCanvas* can = new TCanvas("can","traces",1200,900);
  hbk->Draw();
  
  CSVRow row;
  while (datafile >> row) {
    int eid = atoi(row[0].c_str());
    int stp = atoi(row[19].c_str());
    hstp->Fill(stp);
    if (select_strip>-1 && select_strip==stp) {
      for (auto i=0; i<18; i++)
	gr->SetPoint(i, (float)i, atof(row[i+1].c_str()));
      hbk->SetTitle(Form("%s event %d", reac.c_str(), eid));
      gr->Draw("lp same");
      can->Update();
      can->WaitPrimitive();
    }
  }

  hstp->Draw();
  
  return 0;
}
