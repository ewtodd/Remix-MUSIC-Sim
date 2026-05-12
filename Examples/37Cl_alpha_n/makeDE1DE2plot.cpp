{
  //  disesteem->SetOptStats(0);
  // gSystem->SetObjectStat(0)
  TFile f1("traces_an_stp6_220Torr_lise_79MeV.root");
  TTree* simt1 = (TTree*)f1.Get("simt");
  
  TFile f2("traces_ap_stp6_220Torr_lise_79MeV.root");
  TTree* simt2 = (TTree*)f2.Get("simt");

  TFile f3("traces_aa_stp6_220Torr_lise_79MeV.root");
  TTree* simt3 = (TTree*)f3.Get("simt");

  TFile f4("unreacted_37Cl_220Torr_lise_79MeV.root");
  TTree* simt4 = (TTree*)f4.Get("simt");

#if 0
  gROOT->ProcessLine("simt1->Draw(\"(de_l[6]+de_r[6]):(de_l[7]+de_r[7])>>han(200,0,5,200,0,5)\")");
  han->SetMarkerStyle(7);
  han->SetMarkerColor(4);
  han->SetTitle("");
  han->GetXaxis()->SetTitle("#DeltaE_7+#DeltaE_8 [MeV]");
  
  gROOT->ProcessLine("simt2->Draw(\"(de_l[6]+de_r[6]):(de_l[7]+de_r[7])>>hap(200,0,5,200,0,5)\")");
  hap->SetMarkerStyle(7);
  hap->SetMarkerColor(2);

  gROOT->ProcessLine("simt3->Draw(\"(de_l[6]+de_r[6]):(de_l[7]+de_r[7])>>haa(200,0,5,200,0,5)\")");
  haa->SetMarkerStyle(7);
  haa->SetMarkerColor(921);

  gROOT->ProcessLine("simt4->Draw(\"(de_l[6]+de_r[6]):(de_l[7]+de_r[7])>>hub(200,0,5,200,0,5)\")");
  hub->SetMarkerStyle(7);
  hub->SetMarkerColor(1);
#endif

#if 1
    gROOT->ProcessLine("simt1->Draw(\"(de_l[6]+de_r[6]+de_l[7]+de_r[7]):(de_l[8]+de_r[8]+de_l[9]+de_r[9])>>han(200,4,6.5,200,4,6.5)\")");
  han->SetMarkerStyle(20);
  han->SetMarkerColor(4);
  han->SetTitle("");
  han->GetYaxis()->SetTitle("#DeltaE_{8}+#DeltaE_{9} [MeV]");
  han->GetYaxis()->CenterTitle();
  han->GetXaxis()->SetTitle("#DeltaE_{10}+#DeltaE_{11} [MeV]");
  han->GetXaxis()->CenterTitle();
  
  gROOT->ProcessLine("simt2->Draw(\"(de_l[6]+de_r[6]+de_l[7]+de_r[7]):(de_l[8]+de_r[8]+de_l[9]+de_r[9])>>hap(200,4,6.5,200,4,6.5)\")");
  hap->SetMarkerStyle(20);
  hap->SetMarkerColor(2);

  gROOT->ProcessLine("simt3->Draw(\"(de_l[6]+de_r[6]+de_l[7]+de_r[7]):(de_l[8]+de_r[8]+de_l[9]+de_r[9])>>haa(200,4,6.5,200,4,6.5)\")");
  haa->SetMarkerStyle(20);
  haa->SetMarkerColor(921);

  gROOT->ProcessLine("simt4->Draw(\"(de_l[6]+de_r[6]+de_l[7]+de_r[7]):(de_l[8]+de_r[8]+de_l[9]+de_r[9])>>hub(200,4,6.5,200,4,6.5)\")");
  hub->SetMarkerStyle(20);
  hub->SetMarkerColor(1);
#endif

  han->Draw();
  hap->Draw("same");
  haa->Draw("same");
  hub->Draw("same");

  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(han, "(#alpha,n)", "p");
  Leg->AddEntry(hap, "(#alpha,p)", "p");
  Leg->AddEntry(haa, "(#alpha,#alpha)", "p");
  Leg->AddEntry(hub, "beam", "p");
  Leg->Draw();

}
