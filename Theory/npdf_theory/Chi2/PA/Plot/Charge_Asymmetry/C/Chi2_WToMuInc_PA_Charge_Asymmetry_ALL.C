void Chi2_WToMuInc_PA_Charge_Asymmetry_ALL()
{
//=========Macro generated from canvas: c/c
//=========  (Thu Mar 28 20:02:04 2019) by ROOT version 6.12/07
   TCanvas *c = new TCanvas("c", "c",0,0,4000,1000);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c->Range(-20.02326,-1.5,105.1221,11);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLeftMargin(0.16);
   c->SetRightMargin(0.04);
   c->SetTopMargin(0.08);
   c->SetBottomMargin(0.12);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   Double_t _fx1[388] = {
   0,
   0.25,
   0.5,
   0.75,
   1,
   1.25,
   1.5,
   1.75,
   2,
   2.25,
   2.5,
   2.75,
   3,
   3.25,
   3.5,
   3.75,
   4,
   4.25,
   4.5,
   4.75,
   5,
   5.25,
   5.5,
   5.75,
   6,
   6.25,
   6.5,
   6.75,
   7,
   7.25,
   7.5,
   7.75,
   8,
   8.25,
   8.5,
   8.75,
   9,
   9.25,
   9.5,
   9.75,
   10,
   10.25,
   10.5,
   10.75,
   11,
   11.25,
   11.5,
   11.75,
   12,
   12.25,
   12.5,
   12.75,
   13,
   13.25,
   13.5,
   13.75,
   14,
   14.25,
   14.5,
   14.75,
   15,
   15.25,
   15.5,
   15.75,
   16,
   16.25,
   16.5,
   16.75,
   17,
   17.25,
   17.5,
   17.75,
   18,
   18.25,
   18.5,
   18.75,
   19,
   19.25,
   19.5,
   19.75,
   20,
   20.25,
   20.5,
   20.75,
   21,
   21.25,
   21.5,
   21.75,
   22,
   22.25,
   22.5,
   22.75,
   23,
   23.25,
   23.5,
   23.75,
   24,
   24.25,
   24.5,
   24.75,
   25,
   25.25,
   25.5,
   25.75,
   26,
   26.25,
   26.5,
   26.75,
   27,
   27.25,
   27.5,
   27.75,
   28,
   28.25,
   28.5,
   28.75,
   29,
   29.25,
   29.5,
   29.75,
   30,
   30.25,
   30.5,
   30.75,
   31,
   31.25,
   31.5,
   31.75,
   32,
   32.25,
   32.5,
   32.75,
   33,
   33.25,
   33.5,
   33.75,
   34,
   34.25,
   34.5,
   34.75,
   35,
   35.25,
   35.5,
   35.75,
   36,
   36.25,
   36.5,
   36.75,
   37,
   37.25,
   37.5,
   37.75,
   38,
   38.25,
   38.5,
   38.75,
   39,
   39.25,
   39.5,
   39.75,
   40,
   40.25,
   40.5,
   40.75,
   41,
   41.25,
   41.5,
   41.75,
   42,
   42.25,
   42.5,
   42.75,
   43,
   43.25,
   43.5,
   43.75,
   44,
   44.25,
   44.5,
   44.75,
   45,
   45.25,
   45.5,
   45.75,
   46,
   46.25,
   46.5,
   46.75,
   47,
   47.25,
   47.5,
   47.75,
   48,
   48.25,
   48.5,
   48.75,
   49,
   49.25,
   49.5,
   49.75,
   50,
   50.25,
   50.5,
   50.75,
   51,
   51.25,
   51.5,
   51.75,
   52,
   52.25,
   52.5,
   52.75,
   53,
   53.25,
   53.5,
   53.75,
   54,
   54.25,
   54.5,
   54.75,
   55,
   55.25,
   55.5,
   55.75,
   56,
   56.25,
   56.5,
   56.75,
   57,
   57.25,
   57.5,
   57.75,
   58,
   58.25,
   58.5,
   58.75,
   59,
   59.25,
   59.5,
   59.75,
   60,
   60.25,
   60.5,
   60.75,
   61,
   61.25,
   61.5,
   61.75,
   62,
   62.25,
   62.5,
   62.75,
   63,
   63.25,
   63.5,
   63.75,
   64,
   64.25,
   64.5,
   64.75,
   65,
   65.25,
   65.5,
   65.75,
   66,
   66.25,
   66.5,
   66.75,
   67,
   67.25,
   67.5,
   67.75,
   68,
   68.25,
   68.5,
   68.75,
   69,
   69.25,
   69.5,
   69.75,
   70,
   70.25,
   70.5,
   70.75,
   71,
   71.25,
   71.5,
   71.75,
   72,
   72.25,
   72.5,
   72.75,
   73,
   73.25,
   73.5,
   73.75,
   74,
   74.25,
   74.5,
   74.75,
   75,
   75.25,
   75.5,
   75.75,
   76,
   76.25,
   76.5,
   76.75,
   77,
   77.25,
   77.5,
   77.75,
   78,
   78.25,
   78.5,
   78.75,
   79,
   79.25,
   79.5,
   79.75,
   80,
   80.25,
   80.5,
   80.75,
   81,
   81.25,
   81.5,
   81.75,
   82,
   82.25,
   82.5,
   82.75,
   83,
   83.25,
   83.5,
   83.75,
   84,
   84.25,
   84.5,
   84.75,
   85,
   85.25,
   85.5,
   85.75,
   86,
   86.25,
   86.5,
   86.75,
   87,
   87.25,
   87.5,
   87.75,
   88,
   88.25,
   88.5,
   88.75,
   89,
   89.25,
   89.5,
   89.75,
   90,
   90.25,
   90.5,
   90.75,
   91,
   91.25,
   91.5,
   91.75,
   92,
   92.25,
   92.5,
   92.75,
   93,
   93.25,
   93.5,
   93.75,
   94,
   94.25,
   94.5,
   94.75,
   95,
   95.25,
   95.5,
   95.75,
   96,
   96.25,
   96.5,
   96.75};
   Double_t _fy1[388] = {
   1.034392,
   0,
   0,
   0,
   0.8934321,
   0,
   0,
   0,
   1.309509,
   0,
   0,
   0,
   1.008663,
   0,
   0,
   0,
   1.010455,
   0,
   0,
   0,
   1.041761,
   0,
   0,
   0,
   1.094561,
   0,
   0,
   0,
   1.001026,
   0,
   0,
   0,
   1.133942,
   0,
   0,
   0,
   0.9878416,
   0,
   0,
   0,
   1.30377,
   0,
   0,
   0,
   1.146221,
   0,
   0,
   0,
   1.083601,
   0,
   0,
   0,
   1.192699,
   0,
   0,
   0,
   0.9567562,
   0,
   0,
   0,
   1.032435,
   0,
   0,
   0,
   1.02154,
   0,
   0,
   0,
   0.968696,
   0,
   0,
   0,
   1.059217,
   0,
   0,
   0,
   1.043311,
   0,
   0,
   0,
   1.025057,
   0,
   0,
   0,
   1.015425,
   0,
   0,
   0,
   0.9766455,
   0,
   0,
   0,
   0.9126156,
   0,
   0,
   0,
   1.128356,
   0,
   0,
   0,
   1.015267,
   0,
   0,
   0,
   1.037929,
   0,
   0,
   0,
   1.076193,
   0,
   0,
   0,
   1.012974,
   0,
   0,
   0,
   1.059634,
   0,
   0,
   0,
   1.001882,
   0,
   0,
   0,
   1.032306,
   0,
   0,
   0,
   0.9888612,
   0,
   0,
   0,
   1.121118,
   0,
   0,
   0,
   0.9386541,
   0,
   0,
   0,
   1.042965,
   0,
   0,
   0,
   1.209362,
   0,
   0,
   0,
   1.022243,
   0,
   0,
   0,
   1.016072,
   0,
   0,
   0,
   0.9816091,
   0,
   0,
   0,
   1.061901,
   0,
   0,
   0,
   1.051206,
   0,
   0,
   0,
   1.060435,
   0,
   0,
   0,
   1.110148,
   0,
   0,
   0,
   0.9852564,
   0,
   0,
   0,
   0.9577576,
   0,
   0,
   0,
   1.004293,
   0,
   0,
   0,
   1.047434,
   0,
   0,
   0,
   0.9721477,
   0,
   0,
   0,
   1.034026,
   0,
   0,
   0,
   0.9721989,
   0,
   0,
   0,
   1.004122,
   0,
   0,
   0,
   1.063626,
   0,
   0,
   0,
   1.071541,
   0,
   0,
   0,
   1.001822,
   0,
   0,
   0,
   1.022833,
   0,
   0,
   0,
   1.091822,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0};
   TGraph *graph = new TGraph(388,_fx1,_fy1);
   graph->SetName("");
   graph->SetTitle("Correlation Matrix;PDF set number;#Chi^{2}/ndf");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetFillColor(ci);
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Correlation Matrix",388,0,106.425);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(10);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineStyle(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetTitle("PDF set number");
   Graph_Graph1->GetXaxis()->SetRange(1,365);
   Graph_Graph1->GetXaxis()->CenterTitle(true);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("#Chi^{2}/ndf");
   Graph_Graph1->GetYaxis()->CenterTitle(true);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.7);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ab");
   
   Double_t _fx2[388] = {
   0,
   0.25,
   0.5,
   0.75,
   1,
   1.25,
   1.5,
   1.75,
   2,
   2.25,
   2.5,
   2.75,
   3,
   3.25,
   3.5,
   3.75,
   4,
   4.25,
   4.5,
   4.75,
   5,
   5.25,
   5.5,
   5.75,
   6,
   6.25,
   6.5,
   6.75,
   7,
   7.25,
   7.5,
   7.75,
   8,
   8.25,
   8.5,
   8.75,
   9,
   9.25,
   9.5,
   9.75,
   10,
   10.25,
   10.5,
   10.75,
   11,
   11.25,
   11.5,
   11.75,
   12,
   12.25,
   12.5,
   12.75,
   13,
   13.25,
   13.5,
   13.75,
   14,
   14.25,
   14.5,
   14.75,
   15,
   15.25,
   15.5,
   15.75,
   16,
   16.25,
   16.5,
   16.75,
   17,
   17.25,
   17.5,
   17.75,
   18,
   18.25,
   18.5,
   18.75,
   19,
   19.25,
   19.5,
   19.75,
   20,
   20.25,
   20.5,
   20.75,
   21,
   21.25,
   21.5,
   21.75,
   22,
   22.25,
   22.5,
   22.75,
   23,
   23.25,
   23.5,
   23.75,
   24,
   24.25,
   24.5,
   24.75,
   25,
   25.25,
   25.5,
   25.75,
   26,
   26.25,
   26.5,
   26.75,
   27,
   27.25,
   27.5,
   27.75,
   28,
   28.25,
   28.5,
   28.75,
   29,
   29.25,
   29.5,
   29.75,
   30,
   30.25,
   30.5,
   30.75,
   31,
   31.25,
   31.5,
   31.75,
   32,
   32.25,
   32.5,
   32.75,
   33,
   33.25,
   33.5,
   33.75,
   34,
   34.25,
   34.5,
   34.75,
   35,
   35.25,
   35.5,
   35.75,
   36,
   36.25,
   36.5,
   36.75,
   37,
   37.25,
   37.5,
   37.75,
   38,
   38.25,
   38.5,
   38.75,
   39,
   39.25,
   39.5,
   39.75,
   40,
   40.25,
   40.5,
   40.75,
   41,
   41.25,
   41.5,
   41.75,
   42,
   42.25,
   42.5,
   42.75,
   43,
   43.25,
   43.5,
   43.75,
   44,
   44.25,
   44.5,
   44.75,
   45,
   45.25,
   45.5,
   45.75,
   46,
   46.25,
   46.5,
   46.75,
   47,
   47.25,
   47.5,
   47.75,
   48,
   48.25,
   48.5,
   48.75,
   49,
   49.25,
   49.5,
   49.75,
   50,
   50.25,
   50.5,
   50.75,
   51,
   51.25,
   51.5,
   51.75,
   52,
   52.25,
   52.5,
   52.75,
   53,
   53.25,
   53.5,
   53.75,
   54,
   54.25,
   54.5,
   54.75,
   55,
   55.25,
   55.5,
   55.75,
   56,
   56.25,
   56.5,
   56.75,
   57,
   57.25,
   57.5,
   57.75,
   58,
   58.25,
   58.5,
   58.75,
   59,
   59.25,
   59.5,
   59.75,
   60,
   60.25,
   60.5,
   60.75,
   61,
   61.25,
   61.5,
   61.75,
   62,
   62.25,
   62.5,
   62.75,
   63,
   63.25,
   63.5,
   63.75,
   64,
   64.25,
   64.5,
   64.75,
   65,
   65.25,
   65.5,
   65.75,
   66,
   66.25,
   66.5,
   66.75,
   67,
   67.25,
   67.5,
   67.75,
   68,
   68.25,
   68.5,
   68.75,
   69,
   69.25,
   69.5,
   69.75,
   70,
   70.25,
   70.5,
   70.75,
   71,
   71.25,
   71.5,
   71.75,
   72,
   72.25,
   72.5,
   72.75,
   73,
   73.25,
   73.5,
   73.75,
   74,
   74.25,
   74.5,
   74.75,
   75,
   75.25,
   75.5,
   75.75,
   76,
   76.25,
   76.5,
   76.75,
   77,
   77.25,
   77.5,
   77.75,
   78,
   78.25,
   78.5,
   78.75,
   79,
   79.25,
   79.5,
   79.75,
   80,
   80.25,
   80.5,
   80.75,
   81,
   81.25,
   81.5,
   81.75,
   82,
   82.25,
   82.5,
   82.75,
   83,
   83.25,
   83.5,
   83.75,
   84,
   84.25,
   84.5,
   84.75,
   85,
   85.25,
   85.5,
   85.75,
   86,
   86.25,
   86.5,
   86.75,
   87,
   87.25,
   87.5,
   87.75,
   88,
   88.25,
   88.5,
   88.75,
   89,
   89.25,
   89.5,
   89.75,
   90,
   90.25,
   90.5,
   90.75,
   91,
   91.25,
   91.5,
   91.75,
   92,
   92.25,
   92.5,
   92.75,
   93,
   93.25,
   93.5,
   93.75,
   94,
   94.25,
   94.5,
   94.75,
   95,
   95.25,
   95.5,
   95.75,
   96,
   96.25,
   96.5,
   96.75};
   Double_t _fy2[388] = {
   0,
   1.793194,
   0,
   0,
   0,
   1.745013,
   0,
   0,
   0,
   1.723215,
   0,
   0,
   0,
   1.725077,
   0,
   0,
   0,
   1.770891,
   0,
   0,
   0,
   1.801229,
   0,
   0,
   0,
   1.663974,
   0,
   0,
   0,
   1.687181,
   0,
   0,
   0,
   1.813427,
   0,
   0,
   0,
   1.755357,
   0,
   0,
   0,
   1.672404,
   0,
   0,
   0,
   1.647295,
   0,
   0,
   0,
   1.788041,
   0,
   0,
   0,
   1.706955,
   0,
   0,
   0,
   1.744,
   0,
   0,
   0,
   1.574732,
   0,
   0,
   0,
   1.89641,
   0,
   0,
   0,
   1.744706,
   0,
   0,
   0,
   1.765845,
   0,
   0,
   0,
   1.512909,
   0,
   0,
   0,
   1.988209,
   0,
   0,
   0,
   1.616111,
   0,
   0,
   0,
   1.908541,
   0,
   0,
   0,
   1.847673,
   0,
   0,
   0,
   1.636006,
   0,
   0,
   0,
   1.364016,
   0,
   0,
   0,
   2.36183,
   0,
   0,
   0,
   1.846521,
   0,
   0,
   0,
   1.747626,
   0,
   0,
   0,
   1.863647,
   0,
   0,
   0,
   1.607228,
   0,
   0,
   0,
   1.471379,
   0,
   0,
   0,
   1.994681,
   0,
   0,
   0,
   1.965854,
   0,
   0,
   0,
   1.689656,
   0,
   0,
   0,
   1.776403,
   0,
   0,
   0,
   1.727344,
   0,
   0,
   0,
   2.06311,
   0,
   0,
   0,
   1.590442,
   0,
   0,
   0,
   1.883738,
   0,
   0,
   0,
   1.686677,
   0,
   0,
   0,
   1.899619,
   0,
   0,
   0,
   1.525807,
   0,
   0,
   0,
   1.458276,
   0,
   0,
   0,
   2.058002,
   0,
   0,
   0,
   1.608436,
   0,
   0,
   0,
   1.703153,
   0,
   0,
   0,
   1.706898,
   0,
   0,
   0,
   1.813727,
   0,
   0,
   0,
   1.799925,
   0,
   0,
   0,
   1.696767,
   0,
   0,
   0,
   1.70735,
   0,
   0,
   0,
   1.848835,
   0,
   0,
   0,
   1.786599,
   0,
   0,
   0,
   1.731383,
   0,
   0,
   0,
   1.808656,
   0,
   0,
   0,
   1.75869,
   0,
   0,
   0,
   1.793429,
   0,
   0,
   0,
   1.733756,
   0,
   0,
   0,
   1.905591,
   0,
   0,
   0,
   1.598979,
   0,
   0,
   0,
   1.670103,
   0,
   0,
   0,
   1.80871,
   0,
   0,
   0,
   1.720612,
   0,
   0,
   0,
   1.799289,
   0,
   0,
   0,
   1.986651,
   0,
   0,
   0,
   1.763627,
   0,
   0,
   0,
   2.023563,
   0,
   0,
   0,
   1.652608,
   0,
   0,
   0,
   1.758201,
   0,
   0,
   0,
   1.725373,
   0,
   0,
   0,
   1.848052,
   0,
   0,
   0,
   1.722873,
   0,
   0,
   0,
   1.703275,
   0,
   0,
   0,
   1.910067,
   0,
   0,
   0,
   1.644867,
   0,
   0,
   0,
   1.878179,
   0,
   0,
   0,
   1.678643,
   0,
   0,
   0,
   1.78014,
   0,
   0,
   0,
   1.62599,
   0,
   0,
   0,
   1.801973,
   0,
   0,
   0,
   1.624681,
   0,
   0,
   0,
   2.011649,
   0,
   0,
   0,
   1.766912,
   0,
   0,
   0,
   1.64628,
   0,
   0,
   0,
   1.725798,
   0,
   0,
   0,
   1.699908,
   0,
   0,
   0,
   1.908938,
   0,
   0,
   0,
   1.819354,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0,
   0,
   -1,
   0,
   0};
   graph = new TGraph(388,_fx2,_fy2);
   graph->SetName("");
   graph->SetTitle("");

   ci = TColor::GetColor("#0000ff");
   graph->SetFillColor(ci);
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","",388,0,106.425);
   Graph_Graph2->SetMinimum(-1.336183);
   Graph_Graph2->SetMaximum(2.698013);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->SetLineStyle(0);
   Graph_Graph2->SetMarkerStyle(20);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   graph->Draw("b");
   
   Double_t _fx3[388] = {
   0,
   0.25,
   0.5,
   0.75,
   1,
   1.25,
   1.5,
   1.75,
   2,
   2.25,
   2.5,
   2.75,
   3,
   3.25,
   3.5,
   3.75,
   4,
   4.25,
   4.5,
   4.75,
   5,
   5.25,
   5.5,
   5.75,
   6,
   6.25,
   6.5,
   6.75,
   7,
   7.25,
   7.5,
   7.75,
   8,
   8.25,
   8.5,
   8.75,
   9,
   9.25,
   9.5,
   9.75,
   10,
   10.25,
   10.5,
   10.75,
   11,
   11.25,
   11.5,
   11.75,
   12,
   12.25,
   12.5,
   12.75,
   13,
   13.25,
   13.5,
   13.75,
   14,
   14.25,
   14.5,
   14.75,
   15,
   15.25,
   15.5,
   15.75,
   16,
   16.25,
   16.5,
   16.75,
   17,
   17.25,
   17.5,
   17.75,
   18,
   18.25,
   18.5,
   18.75,
   19,
   19.25,
   19.5,
   19.75,
   20,
   20.25,
   20.5,
   20.75,
   21,
   21.25,
   21.5,
   21.75,
   22,
   22.25,
   22.5,
   22.75,
   23,
   23.25,
   23.5,
   23.75,
   24,
   24.25,
   24.5,
   24.75,
   25,
   25.25,
   25.5,
   25.75,
   26,
   26.25,
   26.5,
   26.75,
   27,
   27.25,
   27.5,
   27.75,
   28,
   28.25,
   28.5,
   28.75,
   29,
   29.25,
   29.5,
   29.75,
   30,
   30.25,
   30.5,
   30.75,
   31,
   31.25,
   31.5,
   31.75,
   32,
   32.25,
   32.5,
   32.75,
   33,
   33.25,
   33.5,
   33.75,
   34,
   34.25,
   34.5,
   34.75,
   35,
   35.25,
   35.5,
   35.75,
   36,
   36.25,
   36.5,
   36.75,
   37,
   37.25,
   37.5,
   37.75,
   38,
   38.25,
   38.5,
   38.75,
   39,
   39.25,
   39.5,
   39.75,
   40,
   40.25,
   40.5,
   40.75,
   41,
   41.25,
   41.5,
   41.75,
   42,
   42.25,
   42.5,
   42.75,
   43,
   43.25,
   43.5,
   43.75,
   44,
   44.25,
   44.5,
   44.75,
   45,
   45.25,
   45.5,
   45.75,
   46,
   46.25,
   46.5,
   46.75,
   47,
   47.25,
   47.5,
   47.75,
   48,
   48.25,
   48.5,
   48.75,
   49,
   49.25,
   49.5,
   49.75,
   50,
   50.25,
   50.5,
   50.75,
   51,
   51.25,
   51.5,
   51.75,
   52,
   52.25,
   52.5,
   52.75,
   53,
   53.25,
   53.5,
   53.75,
   54,
   54.25,
   54.5,
   54.75,
   55,
   55.25,
   55.5,
   55.75,
   56,
   56.25,
   56.5,
   56.75,
   57,
   57.25,
   57.5,
   57.75,
   58,
   58.25,
   58.5,
   58.75,
   59,
   59.25,
   59.5,
   59.75,
   60,
   60.25,
   60.5,
   60.75,
   61,
   61.25,
   61.5,
   61.75,
   62,
   62.25,
   62.5,
   62.75,
   63,
   63.25,
   63.5,
   63.75,
   64,
   64.25,
   64.5,
   64.75,
   65,
   65.25,
   65.5,
   65.75,
   66,
   66.25,
   66.5,
   66.75,
   67,
   67.25,
   67.5,
   67.75,
   68,
   68.25,
   68.5,
   68.75,
   69,
   69.25,
   69.5,
   69.75,
   70,
   70.25,
   70.5,
   70.75,
   71,
   71.25,
   71.5,
   71.75,
   72,
   72.25,
   72.5,
   72.75,
   73,
   73.25,
   73.5,
   73.75,
   74,
   74.25,
   74.5,
   74.75,
   75,
   75.25,
   75.5,
   75.75,
   76,
   76.25,
   76.5,
   76.75,
   77,
   77.25,
   77.5,
   77.75,
   78,
   78.25,
   78.5,
   78.75,
   79,
   79.25,
   79.5,
   79.75,
   80,
   80.25,
   80.5,
   80.75,
   81,
   81.25,
   81.5,
   81.75,
   82,
   82.25,
   82.5,
   82.75,
   83,
   83.25,
   83.5,
   83.75,
   84,
   84.25,
   84.5,
   84.75,
   85,
   85.25,
   85.5,
   85.75,
   86,
   86.25,
   86.5,
   86.75,
   87,
   87.25,
   87.5,
   87.75,
   88,
   88.25,
   88.5,
   88.75,
   89,
   89.25,
   89.5,
   89.75,
   90,
   90.25,
   90.5,
   90.75,
   91,
   91.25,
   91.5,
   91.75,
   92,
   92.25,
   92.5,
   92.75,
   93,
   93.25,
   93.5,
   93.75,
   94,
   94.25,
   94.5,
   94.75,
   95,
   95.25,
   95.5,
   95.75,
   96,
   96.25,
   96.5,
   96.75};
   Double_t _fy3[388] = {
   0,
   0,
   0.8984559,
   0,
   0,
   0,
   0.8810389,
   0,
   0,
   0,
   1.405748,
   0,
   0,
   0,
   0.8914118,
   0,
   0,
   0,
   1.23001,
   0,
   0,
   0,
   0.9279767,
   0,
   0,
   0,
   0.9048479,
   0,
   0,
   0,
   0.8184274,
   0,
   0,
   0,
   1.258698,
   0,
   0,
   0,
   1.229273,
   0,
   0,
   0,
   0.8743006,
   0,
   0,
   0,
   0.8064661,
   0,
   0,
   0,
   1.02134,
   0,
   0,
   0,
   0.9380833,
   0,
   0,
   0,
   0.8770471,
   0,
   0,
   0,
   0.8248462,
   0,
   0,
   0,
   1.057573,
   0,
   0,
   0,
   1.031629,
   0,
   0,
   0,
   0.8258432,
   0,
   0,
   0,
   1.076291,
   0,
   0,
   0,
   0.8371065,
   0,
   0,
   0,
   0.8472071,
   0,
   0,
   0,
   1.021135,
   0,
   0,
   0,
   0.8974882,
   0,
   0,
   0,
   0.9572017,
   0,
   0,
   0,
   0.9718724,
   0,
   0,
   0,
   0.8834793,
   0,
   0,
   0,
   0.9696515,
   0,
   0,
   0,
   0.9502458,
   0,
   0,
   0,
   1.034244,
   0,
   0,
   0,
   0.7826826,
   0,
   0,
   0,
   0.9672693,
   0,
   0,
   0,
   0.9321529,
   0,
   0,
   0,
   1.026007,
   0,
   0,
   0,
   1.001507,
   0,
   0,
   0,
   0.8540208,
   0,
   0,
   0,
   0.9798483,
   0,
   0,
   0,
   0.9492418,
   0,
   0,
   0,
   1.027735,
   0,
   0,
   0,
   1.11495,
   0,
   0,
   0,
   0.8693197,
   0,
   0,
   0,
   0.8624279,
   0,
   0,
   0,
   1.275924,
   0,
   0,
   0,
   0.8403504,
   0,
   0,
   0,
   0.8895638,
   0,
   0,
   0,
   0.91713,
   0,
   0,
   0,
   0.9802552,
   0,
   0,
   0,
   0.945511,
   0,
   0,
   0,
   0.9776251,
   0,
   0,
   0,
   0.8600559,
   0,
   0,
   0,
   1.20792,
   0,
   0,
   0,
   1.114756,
   0,
   0,
   0,
   0.9552594,
   0,
   0,
   0,
   1.066912,
   0,
   0,
   0,
   0.7906763,
   0,
   0,
   0,
   0.9463086,
   0,
   0,
   0,
   0.9023462,
   0,
   0,
   0,
   0.8667895,
   0,
   0,
   0,
   0.958402,
   0,
   0,
   0,
   0.9018064,
   0,
   0,
   0,
   0.9139152,
   0,
   0,
   0,
   0.914645,
   0,
   0,
   0,
   0.8775891,
   0,
   0,
   0,
   0.8233625,
   0,
   0,
   0,
   1.027757,
   0,
   0,
   0,
   0.8763443,
   0,
   0,
   0,
   0.9338654,
   0,
   0,
   0,
   0.8917236,
   0,
   0,
   0,
   1.011749,
   0,
   0,
   0,
   0.9912957,
   0,
   0,
   0,
   0.9056764,
   0,
   0,
   0,
   0.9308974,
   0,
   0,
   0,
   0.8716561,
   0,
   0,
   0,
   1.004349,
   0,
   0,
   0,
   0.802622,
   0,
   0,
   0,
   0.8365791,
   0,
   0,
   0,
   1.182211,
   0,
   0,
   0,
   0.9153915,
   0,
   0,
   0,
   0.861983,
   0,
   0,
   0,
   0.8554623,
   0,
   0,
   0,
   0.9712143,
   0,
   0,
   0,
   1.047847,
   0,
   0,
   0,
   0.8988491,
   0,
   0,
   0,
   0.9906385,
   0,
   0,
   0,
   0.8631291,
   0,
   0,
   0,
   0.8744086,
   0,
   0,
   0,
   0.8702947,
   0,
   0,
   0,
   0.9196787,
   0,
   0,
   0,
   0.8101196,
   0,
   0,
   0,
   0.9426185,
   0,
   0,
   0,
   0.7679548,
   0,
   0,
   0,
   0.813182,
   0,
   0,
   0,
   0.9151202,
   0,
   0,
   0,
   0.9297261,
   0,
   0,
   0,
   0.9535207,
   0,
   0,
   0,
   0.9113499,
   0,
   0,
   0,
   0.9589558,
   0};
   graph = new TGraph(388,_fx3,_fy3);
   graph->SetName("");
   graph->SetTitle("");

   ci = TColor::GetColor("#009900");
   graph->SetFillColor(ci);
   graph->SetFillStyle(1000);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","",388,0,106.425);
   Graph_Graph3->SetMinimum(0);
   Graph_Graph3->SetMaximum(1.546323);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);
   Graph_Graph3->SetLineStyle(0);
   Graph_Graph3->SetMarkerStyle(20);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph3->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph3->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("b");
   
   TLegend *leg = new TLegend(0.22,0.64,0.42,0.84,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("","CT14","lpf");

   ci = TColor::GetColor("#ff0000");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("","CT14+EPPS16","lpf");

   ci = TColor::GetColor("#009900");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("","CT14+nCTEQ15","lpf");

   ci = TColor::GetColor("#0000ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   TLine *line = new TLine(0,1,100,1);
   line->SetLineWidth(3);
   line->Draw();
   TLatex *   tex = new TLatex(0.96,0.9424,"pPb 173.4 nb^{-1}             #sqrt{s_{NN}} = 8.16 TeV");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9424,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.914,0.874,"");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.914,0.802,"Preliminary");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
