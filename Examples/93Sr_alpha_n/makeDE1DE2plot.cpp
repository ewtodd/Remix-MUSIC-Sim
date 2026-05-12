{
  //  gSystem->SetOptStats(0);
  // gSystem->SetObjectStat(0)
  TFile f1("traces_an_stp6_600Torr_lise.root");
  TTree* simt1 = (TTree*)f1.Get("simt");
  
  TFile f2("traces_ap_stp6_600Torr_lise.root");
  TTree* simt2 = (TTree*)f2.Get("simt");

  TFile f3("traces_aa_stp6_600Torr_lise.root");
  TTree* simt3 = (TTree*)f3.Get("simt");

  TFile f4("unreacted_93Sr_600Torr_lise.root");
  TTree* simt4 = (TTree*)f4.Get("simt");

  TFile f5("unreacted_93Y_600Torr_lise.root");
  TTree* simt5 = (TTree*)f5.Get("simt");

#if 1
  simt1->SetAlias("de1","de_l[6]+de_r[6]+de_l[7]+de_r[7]+de_l[8]+de_r[8]+de_l[9]+de_r[9]");
  simt1->SetAlias("de2","de1+de_l[10]+de_r[10]+de_l[11]+de_r[11]+de_l[12]+de_r[12]+de_l[13]+de_r[13]+de_l[14]+de_r[14]+de_l[15]+de_r[15]");
  simt2->SetAlias("de1","de_l[6]+de_r[6]+de_l[7]+de_r[7]+de_l[8]+de_r[8]+de_l[9]+de_r[9]");
  simt2->SetAlias("de2","de1+de_l[10]+de_r[10]+de_l[11]+de_r[11]+de_l[12]+de_r[12]+de_l[13]+de_r[13]+de_l[14]+de_r[14]+de_l[15]+de_r[15]");
  simt3->SetAlias("de1","de_l[6]+de_r[6]+de_l[7]+de_r[7]+de_l[8]+de_r[8]+de_l[9]+de_r[9]");
  simt3->SetAlias("de2","de1+de_l[10]+de_r[10]+de_l[11]+de_r[11]+de_l[12]+de_r[12]+de_l[13]+de_r[13]+de_l[14]+de_r[14]+de_l[15]+de_r[15]");
  simt4->SetAlias("de1","de_l[6]+de_r[6]+de_l[7]+de_r[7]+de_l[8]+de_r[8]+de_l[9]+de_r[9]");
  simt4->SetAlias("de2","de1+de_l[10]+de_r[10]+de_l[11]+de_r[11]+de_l[12]+de_r[12]+de_l[13]+de_r[13]+de_l[14]+de_r[14]+de_l[15]+de_r[15]");
  simt5->SetAlias("de1","de_l[6]+de_r[6]+de_l[7]+de_r[7]+de_l[8]+de_r[8]+de_l[9]+de_r[9]");
  simt5->SetAlias("de2","de1+de_l[10]+de_r[10]+de_l[11]+de_r[11]+de_l[12]+de_r[12]+de_l[13]+de_r[13]+de_l[14]+de_r[14]+de_l[15]+de_r[15]");
  
  gROOT->ProcessLine("simt1->Draw(\"de1:de2>>han(200,145,165,200,55,65)\")");
  han->SetMarkerStyle(7);
  han->SetMarkerColor(4);
  han->SetTitle("");
  han->GetYaxis()->SetTitle("#DeltaE_{7-10} [MeV]");
  han->GetYaxis()->CenterTitle();
  han->GetXaxis()->SetTitle("#DeltaE_{7-18} [MeV]");
  han->GetXaxis()->CenterTitle();
  
  gROOT->ProcessLine("simt2->Draw(\"de1:de2>>hap(200,145,165,200,55,65)\")");
  hap->SetMarkerStyle(7);
  hap->SetMarkerColor(2);

  gROOT->ProcessLine("simt3->Draw(\"de1:de2>>haa(200,145,165,200,55,65)\")");
  haa->SetMarkerStyle(7);
  haa->SetMarkerColor(3);

  gROOT->ProcessLine("simt4->Draw(\"de1:de2>>hub(200,145,165,200,55,65)\")");
  hub->SetMarkerStyle(7);
  hub->SetMarkerColor(1);
  
  gROOT->ProcessLine("simt5->Draw(\"de1:de2>>hubY(200,145,165,200,55,65)\")");
  hubY->SetMarkerStyle(7);
  hubY->SetMarkerColor(921);
#endif

  han->Draw();
  hap->Draw("same");
  haa->Draw("same");
  hub->Draw("same");
  hubY->Draw("same");

  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(han, "(#alpha,n)", "p");
  Leg->AddEntry(hap, "(#alpha,p)", "p");
  Leg->AddEntry(haa, "(#alpha,#alpha)", "p");
  Leg->AddEntry(hub, "93Sr beam", "p");
  Leg->AddEntry(hubY, "93Y beam", "p");
  Leg->Draw();

}
