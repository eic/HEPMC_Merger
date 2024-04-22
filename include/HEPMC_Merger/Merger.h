
#pragma once
#include <iostream>
#include <chrono>
#include <string>
#include <random>
#include "HEPMC_Source.h"
#include <HepMC3/WriterAscii.h>
#include <HepMC3/WriterRootTree.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

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

public:
    // Constructor
  Merger(std::string outputFile, bool rootFormat, double intWindow, int rngSeed, bool verbose, bool squashTime, double bunchSpacing);

  void merge(int nSlices); 
  void printBanner(int nSlices);
  
  // ---------------------------------------------------------------------------
  // Add data source
  // ---------------------------------------------------------------------------
  void addSource(const std::string fileName, double freq, int sourceNo, int beamCorrelation);

  void addEvents(HEPMC_Source&, std::unique_ptr<HepMC3::GenEvent>&);

  // ---------------------------------------------------------------------------
  void squawk(int i);

  std::unique_ptr<HepMC3::GenEvent> mergeSlice(int i);

  // ---------------------------------------------------------------------------
  void insertHepmcEvent( const HepMC3::GenEvent& inevt,
			 std::unique_ptr<HepMC3::GenEvent>& hepSlice, double time=0, int sourceNo=0);
  
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
      
};
