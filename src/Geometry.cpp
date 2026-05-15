#include "Simulator.hpp"

// Anode geometry, hardcoded from agfiles/AnodeGeometry (2019/11/04).
// 20 rows × 2 columns. Strip IDs 0 and 17 are single-column full-width strips;
// 1..16 are split into left/right halves; rows 0 and 19 are dead-layer caps.
void Simulator::LoadHardcodedAnodeGeometry() {
  AnodeRows = 20;
  AnodeCols = 2;
  AnodeDX = new Double_t *[AnodeRows];
  AnodeDY = new Double_t *[AnodeRows];
  AnodeDZ = new Double_t *[AnodeRows];
  AnodeColor = new Short_t *[AnodeRows];
  AnodeSegName = new TString *[AnodeRows];
  AnodeStpID = new Int_t *[AnodeRows];
  for (Int_t r = 0; r < AnodeRows; ++r) {
    AnodeDX[r] = new Double_t[AnodeCols];
    AnodeDY[r] = new Double_t[AnodeCols];
    AnodeDZ[r] = new Double_t[AnodeCols];
    AnodeColor[r] = new Short_t[AnodeCols];
    AnodeSegName[r] = new TString[AnodeCols];
    AnodeStpID[r] = new Int_t[AnodeCols];
    for (Int_t c = 0; c < AnodeCols; ++c) {
      AnodeDX[r][c] = AnodeDY[r][c] = AnodeDZ[r][c] = 0;
      AnodeColor[r][c] = kWhite;
      AnodeSegName[r][c] = "";
      AnodeStpID[r][c] = -1;
    }
  }

  auto put = [&](Int_t r, Int_t c, Int_t id, const char *nm, Double_t dx,
                 Double_t dy, Double_t dz, Short_t color) {
    AnodeStpID[r][c] = id;
    AnodeSegName[r][c] = nm;
    AnodeDX[r][c] = dx;
    AnodeDY[r][c] = dy;
    AnodeDZ[r][c] = dz;
    AnodeColor[r][c] = color;
  };

  // Upstream dead layer (full width).
  put(0, 0, -1, "DeadUS", 9, 10, 3.590, 920);
  // Strip 0 (full width).
  put(1, 0, 0, "S0", 9, 10, 1.578, 416);
  // Strips 1..16 split L/R. Halves sum to 9 cm (not 10). Short half
  // alternates: odd s → R is short (4 cm), L is long (5 cm); even s →
  // L is short (4 cm), R is long (5 cm).
  for (Int_t s = 1; s <= 16; ++s) {
    Int_t row = 1 + s;
    Double_t dxR = (s % 2 == 1) ? 4 : 5;
    Double_t dxL = 9 - dxR;
    char nameR[16], nameL[16];
    std::snprintf(nameR, sizeof(nameR), "S%dC0", s);
    std::snprintf(nameL, sizeof(nameL), "S%dC1", s);
    put(row, 0, s, nameR, dxR, 10, 1.578, 425);
    put(row, 1, s, nameL, dxL, 10, 1.578, 609);
  }
  // Strip 17 (full width).
  put(18, 0, 17, "S17", 9, 10, 1.578, 416);
  // Downstream dead layer. stpid -2 distinguishes from the upstream -1.
  put(19, 0, -2, "DeadDS", 9, 10, 3.522, 920);
}

Int_t Simulator::SetAnode(Short_t Trans, Int_t ELossBins, Float_t MaxELoss) {
  AnodeDepth = AnodeLength = AnodeHeight = 0;
  AnodeRows = AnodeCols = 0;

  if (PrintLevel > 0) {
    Log << "\nmusicsim::SetAnode ************************************************"
        << std::endl;
  }

  LoadHardcodedAnodeGeometry();
  if (verbose_) {
    std::cout << "Anode strips: " << AnodeRows << std::endl;
    std::cout << "Anode columns: " << AnodeCols << std::endl;
  }
  if (AnodeRows <= 0 || AnodeCols <= 0)
    return 0;

  for (Int_t anodeRow = 0; anodeRow < AnodeRows; anodeRow++)
    AnodeDepth += AnodeDZ[anodeRow][0];
  for (Int_t col = 0; col < AnodeCols; col++)
    AnodeLength += AnodeDX[0][col];
  AnodeHeight = AnodeDY[0][0];

  // VolAnode is allocated in all modes — downstream code uses (VolAnode != 0)
  // as a "geometry is configured" sentinel. Workers (Geo==nullptr) leave the
  // entries null; PropagateParticle finds cells via AnodeDX/AnodeDZ.
  VolAnode = new TGeoVolume **[AnodeRows];
  for (Int_t row = 0; row < AnodeRows; ++row) {
    VolAnode[row] = new TGeoVolume *[AnodeCols];
    for (Int_t col = 0; col < AnodeCols; ++col)
      VolAnode[row][col] = nullptr;
  }
  if (Geo) {
    Double_t z0 = 0;
    for (Int_t row = 0; row < AnodeRows; row++) {
      z0 += AnodeDZ[row][0] / 2;
      Double_t x0 = -AnodeLength / 2;
      for (Int_t col = 0; col < AnodeCols; col++) {
        if (AnodeDX[row][col] > 0) {
          VolAnode[row][col] = Geo->MakeBox(
              Form("VolAnode%d%d", row, col), Vacuum, AnodeDX[row][col] / 2,
              AnodeDY[row][col] / 2, AnodeDZ[row][col] / 2);
          VolAnode[row][col]->SetLineColor(AnodeColor[row][col]);
          VolAnode[row][col]->SetTransparency(Trans);
          x0 += AnodeDX[row][col] / 2;
          VolTop->AddNode(VolAnode[row][col], 1,
                          new TGeoTranslation(x0, 0, z0));
          x0 += AnodeDX[row][col] / 2;
        }
      }
      z0 += AnodeDZ[row][0] / 2;
    }
  }

  // Per-strip energy-deposit arrays. The last column (AnodeCols) holds the
  // sum across columns for that strip.
  DeltaEB_ave = new Double_t *[AnodeRows];
  DeltaEB = new Double_t *[AnodeRows];
  DeltaEL = new Double_t *[AnodeRows];
  DeltaEH = new Double_t *[AnodeRows];
  DeltaED1 = new Double_t *[AnodeRows];
  DeltaED2 = new Double_t *[AnodeRows];
  for (Int_t row = 0; row < AnodeRows; row++) {
    DeltaEB_ave[row] = new Double_t[AnodeCols + 1];
    DeltaEB[row] = new Double_t[AnodeCols + 1];
    DeltaEL[row] = new Double_t[AnodeCols + 1];
    DeltaEH[row] = new Double_t[AnodeCols + 1];
    DeltaED1[row] = new Double_t[AnodeCols + 1];
    DeltaED2[row] = new Double_t[AnodeCols + 1];
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      DeltaEB_ave[row][col] = DeltaEB[row][col] = 0;
      DeltaEL[row][col] = DeltaEH[row][col] = 0;
      DeltaED1[row][col] = DeltaED2[row][col] = 0;
    }
  }
  DeltaE_EvaR = new Double_t **[maxEvaporations];
  DeltaE_EvaP = new Double_t **[maxEvaporations];
  for (Int_t er = 0; er < maxEvaporations; er++) {
    DeltaE_EvaR[er] = new Double_t *[AnodeRows];
    DeltaE_EvaP[er] = new Double_t *[AnodeRows];
    for (Int_t row = 0; row < AnodeRows; row++) {
      DeltaE_EvaR[er][row] = new Double_t[AnodeCols + 1];
      DeltaE_EvaP[er][row] = new Double_t[AnodeCols + 1];
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        DeltaE_EvaR[er][row][col] = DeltaE_EvaP[er][row][col] = 0;
    }
  }

  // Trace canvases show the 18 readout strips only; dead-layer dE lives on
  // the MC tree (DeadUS_dE / DeadDS_dE), not the trace. Only allocated for
  // the interactive visualizer — workers run with ctf.Update=0 and would
  // otherwise collide on gROOT's name registry.
  if (ctf.Update) {
    const Int_t NReadout = 18;
    HCTB = new TH2F("HCTB", "Beam", NReadout, -0.5, NReadout - 0.5, ELossBins,
                    0, MaxELoss);
    HCTB->GetXaxis()->SetTitle("Strip ID");
    HCTB->GetXaxis()->CenterTitle();
    HCTB->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCTB->GetYaxis()->CenterTitle();
    HCTB->SetStats(0);

    HCT = new TH2F("HCT", "Column traces", NReadout, -0.5, NReadout - 0.5,
                   ELossBins, 0, MaxELoss);
    HCT->GetXaxis()->SetTitle("Strip ID");
    HCT->GetXaxis()->CenterTitle();
    HCT->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCT->GetYaxis()->CenterTitle();
    HCT->SetStats(0);

    HPT = new TH2F("HPT", "Particle traces", NReadout, -0.5, NReadout - 0.5,
                   ELossBins, 0, MaxELoss);
    HPT->GetXaxis()->SetTitle("Strip ID");
    HPT->GetXaxis()->CenterTitle();
    HPT->GetYaxis()->SetTitle("Energy loss [MeV]");
    HPT->GetYaxis()->CenterTitle();
    HPT->SetStats(0);
  }

  if (Geo)
    Geo->CloseGeometry();
  if (ctf.Update)
    DrawMUSIC(Eve, 85);

  if (verbose_)
    std::cout << "Anode dimensions: " << AnodeLength << "x" << AnodeHeight
              << "x" << AnodeDepth << "cm^3" << std::endl;
  return 1;
}

void Simulator::DrawMUSIC(TEveManager *gEve, Short_t Transparency) {
  if (VolAnode == 0)
    return;
  if (Transparency < 0 || Transparency > 100) {
    std::cout << "Warning: Transparency level must be from 0 to 100."
              << std::endl;
    Transparency = 0;
  }

  TopNode = new TEveGeoTopNode(Geo, Geo->GetTopNode());
  gEve->AddGlobalElement(TopNode);

  TEveArrow *Xaxis = new TEveArrow(20, 0, 0, -10, 0, 0);
  Xaxis->SetName("x axis");
  Xaxis->SetMainColor(kGray);
  Xaxis->SetMainTransparency(65);
  Xaxis->SetTubeR(0.01);
  gEve->AddElement(Xaxis);

  TEveArrow *Yaxis = new TEveArrow(0, 20, 0, 0, -10, 0);
  Yaxis->SetName("y axis");
  Yaxis->SetMainColor(kYellow);
  Yaxis->SetMainTransparency(65);
  Yaxis->SetTubeR(0.01);
  gEve->AddElement(Yaxis);

  TEveArrow *Zaxis = new TEveArrow(0, 0, 1, 0, 0, 0);
  Zaxis->SetName("z axis");
  Zaxis->SetMainColor(kWhite);
  Zaxis->SetTubeR(0.1);
  Zaxis->SetMainTransparency(65);
  gEve->AddElement(Zaxis);

  gEve->Redraw3D(kTRUE);
}
