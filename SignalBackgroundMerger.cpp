#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include <HepMC3/ReaderFactory.h>

#include "argparse/argparse.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;


// =============================================================
/**
    Combine signal and up to three background HEPMC files.
    
    Typical usage:
    ./SignalBackgroundMerger --signalFile dis.hepmc --signalFreq 0 \
            --bg1File hgas.hepmc --bg1Freq 1852 \
            --bg2File egas.hepmc --bg2Freq 1852 \
            --bg3File sr.hepmc
**/    
class SignalBackgroundMerger {
public:

  SignalBackgroundMerger(int argc, char* argv[]) {
    // Parse arguments, print banner, open files, initialize rng
    digestArgs(argc, argv);
    rng.seed( rngSeed );
    banner();
  }
  // ---------------------------------------------------------------------------
  void digestArgs(int argc, char* argv[]) {
    // Handle the command line tedium
    // ArgumentParser is meant to be used in a single function.
    // ArgumentParser internally uses std::string_views,
    // references, iterators, etc.
    // Many of these elements become invalidated after a copy or move.
    argparse::ArgumentParser args ("Merge signal events with up to three background sources.");
    
    args.add_argument("-i", "--signalFile")
      .default_value(std::string("small_ep_noradcor.10x100_q2_10_100_run001.hepmc"))
      .help("Name of the HEPMC file with the signal events");
    
    args.add_argument("-sf", "--signalFreq")
      .default_value(0.0)
      .help("Poisson-mu of the signal frequency in ns. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS rate to randomize.");
    
    args.add_argument("-bg1", "--bg1File")
      .default_value(std::string("small_hgas_100GeV_HiAc_25mrad.Asciiv3.hepmc"))
      .help("Name of the first HEPMC file with background events");
    
    args.add_argument("-bf1", "--bg1Freq")
      .default_value(31347.0)
      .help("Poisson-mu of the first background frequency in ns. Default is the estimated hadron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.");
    
    args.add_argument("-bg2", "--bg2File")
      .default_value(std::string("small_beam_gas_ep_10GeV_foam_emin10keV_vtx.hepmc"))
      .help("Name of the second HEPMC file with background events");
    
    args.add_argument("-bf2", "--bg2Freq")
      .default_value(333.0)
      .help("Poisson-mu of the second background frequency in ns. Default is the estimated electron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.");
    
    args.add_argument("-bg3", "--bg3File")
      .default_value(std::string("small_SR_single_2.5A_10GeV.hepmc"))
      .help("Name of the third HEPMC file with background events");
    
    args.add_argument("-bf3", "--bg3Freq")
      .default_value(0.0)
      .help("Poisson-mu of the third background frequency in ns. Default is 0 to use the weights in the corresponding input file. Set to a value >0 to specify a poisson mu instead.");
    
    args.add_argument("-o", "--outputFile")
      .default_value(std::string(""))
      .help("Specify the output file name. By default it will be auto-generated.");
    
    args.add_argument("-w", "--intWindow")
      .default_value(2000.0)
      .help("Length of the integration window in nanoseconds. Default is 2000.");
    
    args.add_argument("-N", "--nSlices")
      .default_value(-1)
      .help("Number of sampled time slices ('events'). Default is 10. If set to -1, all events in the signal file will be used and background files cycled as needed.");
    
    args.add_argument("--squashTime")
      .default_value(false)
      .implicit_value(true)
      .help("Integration is performed but no time information is associated to vertices.");
    
    args.add_argument("--rngSeed")
      .default_value(0)
      .action([](const std::string& value) { return std::stoi(value); })
      .help("Random seed, default is None");
    
    args.add_argument("-v", "--verbose")
      .default_value(false)
      .implicit_value(true)
      .help("Display details for every slice.");
    
    try {
      args.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
      std::cout << err.what() << std::endl;
      std::cout << args;
      exit(0);
    }
    // Access arguments using args.get method
    signalFile = args.get<std::string>("--signalFile");
    signalFreq = args.get<double>("--signalFreq");
    bg1File = args.get<std::string>("--bg1File");
    bg1Freq = args.get<double>("--bg1Freq");
    bg2File = args.get<std::string>("--bg2File");
    bg2Freq = args.get<double>("--bg2Freq");
    bg3File = args.get<std::string>("--bg3File");
    bg3Freq = args.get<double>("--bg3Freq");
    
    outputFile = args.get<std::string>("--outputFile");
    intWindow  = args.get<double>("--intWindow");
    nSlices    = args.get<int>("--nSlices");
    squashTime = args.get<bool>("--squashTime");
    rngSeed    = args.get<int>("--rngSeed"); 
  }
  
  // ---------------------------------------------------------------------------
  void banner() {
    // Print banner
    std::cout << "==================================================================" << std::endl;
    std::cout << "=== EPIC HEPMC MERGER ===" << std::endl;
    std::cout << "authors: Benjamen Sterwerf* (bsterwerf@berkeley.edu), Kolja Kauder** (kkauder@bnl.gov), Reynier Cruz-Torres***" << std::endl;
    std::cout << "* University of California, Berkeley" << std::endl;
    std::cout << "** Brookhaven National Laboratory" << std::endl;
    std::cout << "*** formerly Lawrence Berkeley National Laboratory" << std::endl;
    std::cout << "\nFor more information, run \n./signal_background_merger --help" << std::endl;
    
    std::string freqTerm = signalFreq > 0 ? std::to_string(signalFreq) + " ns" : "(one event per time slice)";
    std::cout << "Signal events file and frequency:\n";
    std::cout << "\t- " << signalFile << "\t" << freqTerm << "\n";
    
    std::cout << "\nBackground files and their respective frequencies:\n";
    if (!bg1File.empty()) {
      freqTerm = bg1Freq > 0 ? std::to_string(bg1Freq) + " ns" : "(from weights)";
      std::cout << "\t- " << bg1File << "\t" << freqTerm << "\n";
    }
    if (!bg2File.empty()) {
      freqTerm = bg2Freq > 0 ? std::to_string(bg2Freq) + " ns" : "(from weights)";
      std::cout << "\t- " << bg2File << "\t" << freqTerm << "\n";
    }
    if (!bg3File.empty()) {
      freqTerm = bg3Freq > 0 ? std::to_string(bg3Freq) + " ns" : "(from weights)";
      std::cout << "\t- " << bg3File << "\t" << freqTerm << "\n";
    }	
  }
  
  // ---------------------------------------------------------------------------  
  void makeDicts(const std::string& fileName, double freq, bool signal=false) {
    if (fileName.empty()) return;
    
    std::shared_ptr<HepMC3::Reader> adapter;
    try {
      adapter = HepMC3::deduce_reader(fileName);
      if (!adapter) {
	throw std::runtime_error("Failed to open file");
      }
    } catch (const std::runtime_error& e) {
      std::cerr << "Opening " << fileName << " failed: " << e.what() << std::endl;
      exit(1);
    }
    
    if (signal) {
      sigDict[fileName] = {adapter, freq};
      return;
    }
    
    if (freq <= 0) {
      std::cout << "Reading in all events from " << fileName << std::endl;
      std::vector<HepMC3::GenEvent> events;
      std::vector<double> weights;
      double weight;
      HepMC3::GenEvent event;

      HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
      adapter->read_event(evt);

	
      // while (File >> event >> weight) {
      // 	if (weight > 0) {
      // 	  events.push_back(event);
      // 	  weights.push_back(weight);
      // 	}
      // }
      // File.close();

      std::sort(weights.begin(), weights.end());
      double avgRate = std::accumulate(weights.begin(), weights.end(), 0.0) / weights.size();
      avgRate *= 1e-9;
      std::cout << "Average rate is " << avgRate << " GHz" << std::endl;
      for (auto& weight : weights) {
	weight /= avgRate;
      }
      weightDict[fileName] = {events, weights};
      return;
    }

    // Not signal and not weighted --> update freqDict
    freqDict[fileName] = {adapter, freq};
  }
  // ---------------------------------------------------------------------------

private:
  std::default_random_engine rng;
  string signalFile, bg1File, bg2File, bg3File;
  double signalFreq, bg1Freq, bg2Freq, bg3Freq;
  string outputFile;
  double intWindow;
  int nSlices; // should be long, but argparse cannot read that
  bool squashTime;
  int rngSeed;  // should be unsigned, but argparse cannot read that

  std::map<std::string, std::pair< std::shared_ptr<HepMC3::Reader>,double> > sigDict;
  std::map<std::string, std::pair<std::vector<HepMC3::GenEvent>, std::vector<double>>> weightDict;
  std::map<std::string, std::pair<std::shared_ptr<HepMC3::Reader>,double> > freqDict;

};

// =============================================================
int main(int argc, char* argv[]) {
  cout << argc << "  " << argv[0] << endl;

  // Create an instance of SignalBackgroundMerger
  SignalBackgroundMerger merger (argc, argv);


  
  
  return 0;
}
