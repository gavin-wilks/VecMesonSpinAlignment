#include <TString.h>
#include <TFile.h>
#include <TGraph.h>

#include <string>
#include <ostream>

void plotTheory()
{
  string par[3] = {"120","128","136"};
  TGraph *g_Theory[3];

  for(int i_par= 0; i_par< 3; ++i_par)
  {
    string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Theory%s.txt",par[i_par].c_str());
    cout << "inputfilename = " << inputfile.c_str() << endl;
    g_Theory[i_par] = new TGraph();

    FILE *fp = fopen(inputfile.c_str(), "r");
    if(fp == NULL)
    {
      perror("Error opening file!");
    }
    else
    {
      float energy, rho00;
      char line[100];
      Int_t line_counter = 0;

      while(fgets(line,100,fp))
      {
	sscanf(&line[0], "%f %f", &energy, &rho00);
	cout << "energy = " << energy << ", rho00 = " << rho00 << endl;
	g_Theory[i_par]->SetPoint(line_counter,energy,rho00);
	line_counter++;
      }
    }
  }

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Theory.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_par = 0; i_par < 3; ++i_par)
  {
    string GraphName = Form("g_Theory_%s",par[i_par].c_str());
    g_Theory[i_par]->SetName(GraphName.c_str());
    g_Theory[i_par]->Write();
  }
  File_OutPut->Close();
}
