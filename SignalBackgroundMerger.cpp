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
#include <chrono>
#include <cmath>
#include <random>
#include <sys/resource.h>

#include <HepMC3/ReaderFactory.h>
#include <HepMC3/WriterAscii.h>
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"

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
    auto t0 = std::chrono::high_resolution_clock::now();

    // Parse arguments, print banner, open files, initialize rng
    digestArgs(argc, argv);
    rng.seed( rngSeed );
    banner();
    if (outputFile != "" ) {
      outputFileName = outputFile;
    } else {
      outputFileName = nameGen();
    }
    cout << "Writing to " << outputFileName << endl;

    makeDicts ( signalFile, signalFreq, true );
    makeDicts ( bg1File, bg1Freq );
    makeDicts ( bg2File, bg2Freq );
    makeDicts ( bg3File, bg3Freq );

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initiation time: " << std::round(std::chrono::duration<double, std::chrono::seconds::period>(t1 - t0).count()) << " sec" << std::endl;
    std::cout << "\n==================================================================\n" << std::endl;

  }

  void merge(){
    auto t1 = std::chrono::high_resolution_clock::now();

    // Open output file    
    std::shared_ptr<HepMC3::Writer> f;
    if (rootFormat)
      f = std::make_shared<HepMC3::WriterRootTree>(outputFileName);
    else
      f = std::make_shared<HepMC3::WriterAscii>(outputFileName);
    
    // Slice loop
    int i = 0;
    for (i = 0; i<nSlices; ++i ) {
      if (i % 100 == 0 || verbose ) squawk(i);
      auto hepSlice = mergeSlice(i);
      if (!hepSlice) {
	std::cout << "Exhausted signal source." << std::endl;
	break;
      }
      hepSlice->set_event_number(i);
      f->write_event(*hepSlice);
    }
    std::cout << "Finished all requested slices." << std::endl;

    int slicesDone = i;
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Slice loop time: " << std::round(std::chrono::duration<double, std::chrono::minutes::period>(t2 - t1).count()) << " min" << std::endl;
    std::cout << " -- " << std::round(std::chrono::duration<double>(t2 - t1).count() / i) << " sec / slice" << std::endl;

    // clean up, close all files
    auto it = sigDict.find(signalFile);
    if (it != sigDict.end()) {
      it->second.first->close();
    }
    for (auto& it : freqDict) {
      it.second.first->close();
    }
    f->close();
	
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
      .scan<'g', double>()
      .help("Poisson-mu of the signal frequency in ns. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS rate to randomize.");
    
    args.add_argument("-bg1", "--bg1File")
      .default_value(std::string("small_hgas_100GeV_HiAc_25mrad.Asciiv3.hepmc"))
      .help("Name of the first HEPMC file with background events");
    
    args.add_argument("-bf1", "--bg1Freq")
      .default_value(31347.0)
      .scan<'g', double>()
      .help("Poisson-mu of the first background frequency in ns. Default is the estimated hadron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.");
    
    args.add_argument("-bg2", "--bg2File")
      .default_value(std::string("small_beam_gas_ep_10GeV_foam_emin10keV_vtx.hepmc"))
      .help("Name of the second HEPMC file with background events");
    
    args.add_argument("-bf2", "--bg2Freq")
      .default_value(333.0)
      .scan<'g', double>()
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

    args.add_argument("-r", "--rootFormat")
      .default_value(false)
      .implicit_value(true)
      .help("Use hepmc.root output format, default is false.");
	
    args.add_argument("-w", "--intWindow")
      .default_value(2000.0)
      .help("Length of the integration window in nanoseconds. Default is 2000.");
    
    args.add_argument("-N", "--nSlices")
      .default_value(-1)
      .scan<'i', int>()
      .help("Number of sampled time slices ('events'). Default is -1. If set to -1, all events in the signal file will be used and background files cycled as needed.");
    
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
    rootFormat = args.get<bool>("--rootFormat");
    intWindow  = args.get<double>("--intWindow");
    nSlices    = args.get<int>("--nSlices");
    squashTime = args.get<bool>("--squashTime");
    rngSeed    = args.get<int>("--rngSeed");
    verbose    = args.get<bool>("--verbose");

      
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

    std::cout << "Number of Slices:" << nSlices << endl;
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

      std::vector<weightedEvent> weightedEvents;
	// HepMC3::GenEvent event;


      while(!adapter->failed()) {
	HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
	adapter->read_event(evt);
	// HepMC3::Print::listing(evt);
	// HepMC3::Print::content(evt);

	// remove events with 0 weight - note that this does change avgRate = <weight> (by a little)
	if (evt.weight() > 0){
	  weightedEvents.push_back({evt,evt.weight()});
	}
      }
      adapter->close();
      
      double avgRate = 0.0;
      for ( auto w : weightedEvents ){ avgRate += w.second;}
      avgRate /= weightedEvents.size();
      avgRate *= 1e-9; // convert to 1/ns == GHz
      std::cout << "Average rate is " << avgRate << " GHz" << std::endl;
      for ( auto w : weightedEvents ) {
	w.second /= avgRate;
      }

      // Optional to sort, may speed up random drawing
      std::sort(weightedEvents.begin(), weightedEvents.end(), [](auto &left, auto &right) {
	  return left.second < right.second;
	});
      
      weightDict[fileName] = {weightedEvents, avgRate};

      return;
    }

    // Not signal and not weighted --> update freqDict
    freqDict[fileName] = {adapter, freq};
  }


   // ---------------------------------------------------------------------------
  bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
      return false;
    }
  }
  
  // ---------------------------------------------------------------------------
  std::string nameGen() {
    std::string name = signalFile;
    if (nSlices > 0) {
        size_t pos = name.find(".hepmc");
        if (pos != std::string::npos) {
	  name.replace(pos, 6, "_n_" + std::to_string(nSlices) + ".hepmc");
        }
    }
    if ( rootFormat && !hasEnding(name,".root")){
      name.append(".root");
    }
    name = "bgmerged_" + name;
    return name;
}

  // ---------------------------------------------------------------------------
  void squawk(int i) {
    struct rusage r_usage;
    // NOTE: Reported in kB on Linux, bytes in Mac/Darwin
    // Could try to explicitly catch __linux__ as well
    // Unclear in BSD, I've seen conflicting reports
    float mbsize = 1024;
#ifdef __MACH__
    mbsize = 1024 * 1024;
#endif
  
    getrusage(RUSAGE_SELF, &r_usage);
    
    std::cout << "Working on slice " << i + 1 << std::endl;
    std::cout << "Resident Memory " << r_usage.ru_maxrss / mbsize << " MB" << std::endl;
  }
  // ---------------------------------------------------------------------------

  std::unique_ptr<HepMC3::GenEvent> mergeSlice(int i) {
    auto hepSlice = std::make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
    
    addFreqEvents(signalFile, hepSlice, true);
    
    for (const auto& fileName : freqDict) {
      addFreqEvents(fileName.first, hepSlice);
    }
    
    for (const auto& fileName : weightDict) {
      addWeightedEvents(fileName.first, hepSlice);
    }

    return hepSlice;
  };

  // ---------------------------------------------------------------------------

  void addFreqEvents(std::string fileName, std::unique_ptr<HepMC3::GenEvent>& hepSlice, bool signal = false) {
    double freq;
    std::shared_ptr<HepMC3::Reader> adapter;
    
    if (signal) {
      std::tie(adapter, freq) = sigDict[fileName];
    } else {
      std::tie(adapter, freq) = freqDict[fileName];
    }
    
    // First, create a timeline
    // Signals can be different
    std::vector<double> slice;
    std::uniform_real_distribution<> uni(0, intWindow);
    if (freq == 0) {
      if (!signal) {
	std::cerr << "frequency can't be 0 for background files" << std::endl;
	exit(1);
      }
      // exactly one signal event, at an arbitrary point
      slice.push_back(uni(rng));
    } else {
      // Generate poisson-distributed times to place events
      slice = poissonTimes(freq, intWindow);
    }
    
    if ( verbose) std::cout << "Placing " << slice.size() << " events from " << fileName << std::endl;
    if (slice.empty()) return;
    
    // Insert events at all specified locations
    for (double time : slice) {
      if(adapter->failed()) {
	try{
	  if (signal) {
	    // Exhausted signal events
	    throw std::ifstream::failure("EOF");
	  } else {
	    // background file reached its end, reset to the start
	    std::cout << "Cycling back to the start of " << fileName << std::endl;
	    adapter->close();
	    adapter = HepMC3::deduce_reader(fileName);	    
	  }
	} catch (std::ifstream::failure& e) {
	  continue; // just need to suppress the error
        }
      }

      HepMC3::GenEvent inevt;
      adapter->read_event(inevt);
      
      // Unit conversion
      double timeHepmc = c_light * time;

      std::vector<HepMC3::GenParticlePtr> particles;
      std::vector<HepMC3::GenVertexPtr> vertices;

      // Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
      for (auto& vertex : inevt.vertices()) {
	HepMC3::FourVector position = vertex->position();
	if (!squashTime) {
	  position.set_t(position.t() + timeHepmc);
	}
	auto v1 = std::make_shared<HepMC3::GenVertex>(position);
	vertices.push_back(v1);
      }

      // copies the particles and attaches them to their corresponding vertices
      for (auto& particle : inevt.particles()) {
	HepMC3::FourVector momentum = particle->momentum();
	int status = particle->status();
	int pid = particle->pid();
	auto p1 = std::make_shared<HepMC3::GenParticle> (momentum, pid, status);
	p1->set_generated_mass(particle->generated_mass());
	particles.push_back(p1);
	
	// since the beam particles do not have a production vertex they cannot be attached to a production vertex
	if (particle->production_vertex()->id() < 0) {
	  int production_vertex = particle->production_vertex()->id();
	  vertices[abs(production_vertex) - 1]->add_particle_out(p1);
	  hepSlice->add_particle(p1);
	}
	
	// Adds particles with an end vertex to their end vertices
	if (particle->end_vertex()) {
	  int end_vertex = particle->end_vertex()->id();
	  vertices[abs(end_vertex) - 1]->add_particle_in(p1);
	}
      }

      // Adds the vertices with the attached particles to the event
      for (auto& vertex : vertices) {
	hepSlice->add_vertex(vertex);
      }
    }

    return;
  }

  // ---------------------------------------------------------------------------

  void addWeightedEvents(std::string fileName, std::unique_ptr<HepMC3::GenEvent>& hepSlice, bool signal = false) {
    std::vector<weightedEvent> events;
    double avgRate;

    std::tie(events, avgRate) = weightDict[fileName];

    // How many events? Assume Poisson distribution
    int nEvents;
    std::exponential_distribution<> d( intWindow * avgRate );
    while (true) {
      nEvents = static_cast<int>(d(rng));
      if (nEvents > events.size()) {
	std::cout << "WARNING: Trying to place " << nEvents << " events from " << fileName
		  << " but the file doesn't have enough. Rerolling, but this is not physical." << std::endl;
	continue;
      }
      break;
    }
    if (verbose) std::cout << "Placing " << nEvents << " events from " << fileName << std::endl;

    // // Get events
    // std::vector<weightedEvent> toPlace = rng.choice(events, nEvents, false);

    // // Place at random times
    // std::vector<double> times;
    // std::uniform_real_distribution<> uni(0, intWindow);
    // if (!squashTime) {
    //   times = rng.uniform(0, intWindow, nEvents);
    // }

    // for (Event& inevt : toPlace) {
    //     double Time = 0;
    //     if (!squash) {
    //         Time = times.back();
    //         times.pop_back();
    //     }

    //     std::vector<HepMC3::GenParticle> particles;
    //     std::vector<HepMC3::GenVertex> vertices;

    //     // Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
    //     for (auto& vertex : inevt.vertices) {
    //         HepMC3::FourVector position = vertex->position();
    //         if (!squash) {
    //             // Unit conversion
    //             double TimeHepmc = c_light * Time;
    //             position.set_t(position.t() + TimeHepmc);
    //         }
    //         HepMC3::GenVertex v1(position);
    //         vertices.push_back(v1);
    //     }

    //     // copies the particles and attaches them to their corresponding vertices
    //     for (auto& particle : inevt.particles) {
    //         HepMC3::FourVector momentum = particle->momentum();
    //         int status = particle->status();
    //         int pid = particle->pid();
    //         HepMC3::GenParticle p1(momentum, pid, status);
    //         p1.set_generated_mass(particle->generated_mass());
    //         particles.push_back(p1);

    //         // since the beam particles do not have a production vertex they cannot be attached to a production vertex
    //         if (particle->production_vertex()->id() < 0) {
    //             int production_vertex = particle->production_vertex()->id();
    //             vertices[abs(production_vertex) - 1].add_particle_out(p1);
    //             hepSlice.add_particle(p1);
    //         }

    //         // Adds particles with an end vertex to their end vertices
    //         if (particle->end_vertex()) {
    //             int end_vertex = particle->end_vertex()->id();
    //             vertices[abs(end_vertex) - 1].add_particle_in(p1);
    //         }
    //     }

    //     // Adds the vertices with the attached particles to the event
    //     for (auto& vertex : vertices) {
    //         hepSlice.add_vertex(vertex);
    //     }
    // }

    return;
}

  // ---------------------------------------------------------------------------

  std::vector<double> poissonTimes(double mu, double endTime) {
    std::exponential_distribution<> exp(1.0 / mu);
    
    double t = 0;
    std::vector<double> ret;
    while (true) {
      double delt = exp(rng);
      t += delt;
      if (t >= endTime) {
	break;
      }
      ret.push_back(t);
    }
    return ret;
}
  // ---------------------------------------------------------------------------

  
private:
  std::mt19937 rng;
  string signalFile, bg1File, bg2File, bg3File;
  double signalFreq, bg1Freq, bg2Freq, bg3Freq;
  string outputFile;
  string outputFileName;
  bool rootFormat;
  double intWindow;
  int nSlices; // should be long, but argparse cannot read that
  bool squashTime;
  int rngSeed;  // should be unsigned, but argparse cannot read that
  bool verbose;
  
  std::map<std::string, std::pair< std::shared_ptr<HepMC3::Reader>,double> > sigDict;
  std::map<std::string, std::pair< std::shared_ptr<HepMC3::Reader>,double> > freqDict;

  typedef std::pair<HepMC3::GenEvent,double> weightedEvent;
  std::map<std::string, std::pair<std::vector<weightedEvent>,double>> weightDict;

  const double c_light = 299.792458; // speed of light = 299.792458 mm/ns to get mm

};

// =============================================================
int main(int argc, char* argv[]) {

  auto t0 = std::chrono::high_resolution_clock::now();
  // Create an instance of SignalBackgroundMerger
  SignalBackgroundMerger sbm (argc, argv);
  sbm.merge();

  std::cout << "\n==================================================================\n";
  std::cout << "Overall running time: " << std::round(std::chrono::duration<double, std::chrono::minutes::period>(std::chrono::high_resolution_clock::now() - t0).count()) << " min" << std::endl;
  
  
  return 0;
}
