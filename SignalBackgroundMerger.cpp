#include <iostream>

#include "Merger.h"
#include "Merger_parser.h"


// =============================================================
int main(int argc, char* argv[]) {

  auto t0 = std::chrono::high_resolution_clock::now();

  // Parse the command line arguments
  argparse::ArgumentParser* args = digestArgs(argc, argv);

  // Create an instance of SignalBackgroundMerger
  std::string outputFile   = args->get<std::string>("--outputFile");
  bool        rootFormat   = args->get<bool>("--rootFormat");
  double      intWindow    = args->get<double>("--intWindow");
  int         rngSeed      = args->get<int>("--rngSeed");
  bool        verbose      = args->get<bool>("--verbose");
  bool        squashTime   = args->get<bool>("--squashTime");
  double      bunchSpacing = args->get<double>("--bunchSpacing");
  Merger sbm(outputFile, rootFormat, intWindow, rngSeed, verbose, squashTime, bunchSpacing);

  // Add sources
  std::string signalFile = args->get<std::string>("--signalFile");
  double signalFreq = args->get<double>("--signalFreq");
  int signalCorrelation = args->get<int>("--signalCorrelation");
  sbm.addSource(signalFile, signalFreq,0,signalCorrelation);

  std::string bg1File = args->get<std::string>("--bg1File");
  double bg1Freq = args->get<double>("--bg1Freq");
  int bg1Correlation = args->get<int>("--bg1Correlation");
  sbm.addSource(bg1File, bg1Freq,1,bg1Correlation);
  
  std::string bg2File = args->get<std::string>("--bg2File");
  double bg2Freq = args->get<double>("--bg2Freq");
  int bg2Correlation = args->get<int>("--bg2Correlation");
  sbm.addSource(bg2File, bg2Freq,2,bg2Correlation);
  
  std::string bg3File = args->get<std::string>("--bg3File");
  double bg3Freq = args->get<double>("--bg3Freq");
  int bg3Correlation = args->get<int>("--bg3Correlation");
  sbm.addSource(bg3File, bg3Freq,3,bg3Correlation);

  // Merge the sources
  int  nSlices    = args->get<int>("--nSlices");
  sbm.merge(nSlices);

  // // DEBUG
  // auto f = new TFile("f.root","RECREATE");
  // auto p = new TH1I("p","poisson",100,-1,99);
  // auto e = new TH1I("e","exponential",100,-1,99);
  // int N = 100000;
  // float rate = 20;
  // float length = 2;
  // std::poisson_distribution<> d( rate * length );
  
  // for (int i=0; i< N; ++i){
  //   auto np = d(sbm.rng);
  //   p->Fill(np);
  //   auto pTimes = sbm.poissonTimes(rate, length);
  //   auto ne = pTimes.size();
  //   e->Fill(ne);
  //   // cout << np << "  " << ne << endl;
  // }
  // sbm.p->Write();
  // sbm.e->Write();  
  // sbm.f->Close();

  std::cout << "\n==================================================================\n";
  std::cout << "Overall running time: " << std::round(std::chrono::duration<double, std::chrono::minutes::period>(std::chrono::high_resolution_clock::now() - t0).count()) << " min" << std::endl;
  
  
  return 0;
}
