
using std::cout;
using std::cerr;
using std::endl;
using std::string;

#include "HEPMC_Source.h"

// =============================================================
/**
    Combine signal and up to three background HEPMC files.
    
    Typical usage:
    ./SignalBackgroundMerger --signalFile dis.hepmc --signalFreq 0 \
            --bg1File hgas.hepmc --bg1Freq 1852 \
            --bg2File egas.hepmc --bg2Freq 1852 \
            --bg3File sr.hepmc
**/    
class Merger {

private:
  // more private data at the end; pulling these more complicated objects up for readability
  
  // std::map<std::string, std::pair< std::shared_ptr<HepMC3::Reader>,double> > sigDict;
  // std::map<std::string, std::pair< std::shared_ptr<HepMC3::Reader>,double> > freqDict;
  

public:
    // Constructor
  Merger(std::string outputFile, bool rootFormat, double intWindow, int rngSeed, bool verbose, bool squashTime, double bunchSpacing);

  void merge(int nSlices); 
  void printBanner(int nSlices);
  
  // ---------------------------------------------------------------------------
  // Add data source
  // ---------------------------------------------------------------------------
  void addSource(const std::string fileName, double freq, int sourceNo);

  // ---------------------------------------------------------------------------
  void squawk(int i);

  std::unique_ptr<HepMC3::GenEvent> mergeSlice(int i);

  // ---------------------------------------------------------------------------
  void insertHepmcEvent( const HepMC3::GenEvent& inevt,
			 std::unique_ptr<HepMC3::GenEvent>& hepSlice, double time=0);
  
  private:
    std::shared_ptr<HepMC3::Writer> outFile;
    std::mt19937 rng;
    
    bool   rootFormat;
    double m_intWindow;
    bool   m_squashTime;
    bool   m_verbose;
    double m_bunchSpacing;
    std::string outputFileName;
    
    const double c_light = 299.792458; // speed of light = 299.792458 mm/ns to get mm

    std::vector<HEPMC_Source> sources;
    

  // // DEBUG  
  // TFile *f;
  // TH1I* p;
  // TH1I* e;

  
};
