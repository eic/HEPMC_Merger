#ifdef __MACH__
#include <mach/mach.h>
#endif

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
#include <tuple>
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
    Combine signal and up to four background HEPMC files.
    
    Typical usage:
    ./SignalBackgroundMerger --signalFile dis.hepmc3.tree.root --signalFreq 0 \
            --bgFile hgas.hepmc3.tree.root 2000 0 2000 \
	    --bgFile egastouschk.hepmc3.tree.root 20 0 3000 \
            --bgFile egascouloumb.hepmc3.tree.root 20 0 4000 \
	    --bgFile egasbrems.hepmc3.tree.root 20 0 5000 \
            --bgFile synrad.hepmc3.tree.root 25 0 6000
**/    

struct BackgroundConfig {
    std::string file;
    double frequency=0;
    int skip=0;
    int status=0;
} ;

class SignalBackgroundMerger {

private:
  // more private data at the end; pulling these more complicated objects up for readability
  std::shared_ptr<HepMC3::Reader> sigAdapter;
  double sigFreq = 0;
  int sigStatus = 0;
  std::map<std::string, std::shared_ptr<HepMC3::Reader>> freqAdapters;
  std::map<std::string, double> freqs;
  std::map<std::string, int> baseStatuses;

  std::map<std::string,
	  std::tuple<std::vector<HepMC3::GenEvent>,
		      std::piecewise_constant_distribution<>,
		      double>
	   > weightDict;

  // just keep count of some numbers, could be more sophisticated
  typedef struct{
    long eventCount;
    long particleCount;
  } stats;
  std::map<std::string, stats > infoDict;

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
    std::cout << "\n==================================================================\n";
    cout << "Writing to " << outputFileName << endl;

    PrepData ( signalFile, signalFreq, signalSkip, signalStatus, true );
    for (const auto& bg : backgroundFiles) {
	PrepData ( bg.file, bg.frequency, bg.skip, bg.status, false );
    }
    
    
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initiation time: " << std::round(std::chrono::duration<double, std::chrono::seconds::period>(t1 - t0).count()) << " sec" << std::endl;
    std::cout << "\n==================================================================\n" << std::endl;

  }

  // Helper to parse raw strings into BackgroundConfig structs
  std::vector<BackgroundConfig>
  parse_backgrounds(const std::vector<std::string> &raw_args_list) {
    std::vector<BackgroundConfig> backgrounds;
    auto is_pure_integer = [](const std::string &str) {
      if (str.empty()) return false;
      for (char c : str) {
        if (!std::isdigit(c)) return false;
      }
      return true;
    };
	  
    // Group strings into sets of 2-4 arguments per background
    for (size_t i = 0; i < raw_args_list.size();) {
      // Determine how many arguments this background has
      size_t args_count = 2; // minimum

      // Look ahead to see if next strings can be parsed as numbers (skip/status)
      if (i + 2 < raw_args_list.size() && is_pure_integer(raw_args_list[i + 2])) {
          args_count = 3;

          if (i + 3 < raw_args_list.size() && is_pure_integer(raw_args_list[i + 3])) {
              args_count = 4;
          }
      }

      // Ensure we don't go beyond the vector bounds
      if (i + args_count > raw_args_list.size()) {
        args_count = raw_args_list.size() - i;
      }

      if (args_count < 2) {
        throw std::runtime_error("Background file " +
                                 std::to_string(backgrounds.size()) +
                                 " must have at least 2 arguments");
      }

      try {
        BackgroundConfig bg;
        bg.file = raw_args_list[i];
        bg.frequency = std::stod(raw_args_list[i + 1]);
        bg.skip = (args_count > 2) ? std::stoi(raw_args_list[i + 2]) : 0;
        bg.status = (args_count > 3) ? std::stoi(raw_args_list[i + 3]) : 0;
        backgrounds.push_back(bg);
      } catch (const std::exception &e) {
        throw std::runtime_error("Error parsing background file " +
                                 std::to_string(backgrounds.size()) + ": " +
                                 e.what());
      }

      i += args_count;
    }

    return backgrounds;
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
      if (i % 1000 == 0 || verbose ) squawk(i);
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
    std::cout << " -- " << std::round(std::chrono::duration<double, std::chrono::microseconds::period>(t2 - t1).count() / i) << " us / slice" << std::endl;

    for (auto info : infoDict) {
      std::cout << "From " << info.first << std::endl;
      std::cout << "  placed " << info.second.eventCount << " events" << std::endl; 
      std::cout << "  --> on average " << std::setprecision(3) << info.second.eventCount / float(nSlices) << std::endl;
      std::cout << "  placed " << info.second.particleCount << " final state particles" << std::endl;
      std::cout << "  --> on average " << std::setprecision(3) << info.second.particleCount / float(nSlices) << std::endl;
      
    }

    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);

    // NOTE: Reported in kB on Linux, bytes in Mac/Darwin
    // Could try to explicitly catch __linux__ as well
    // Unclear in BSD, I've seen conflicting reports
#ifdef __MACH__
    float mbsize = 1024 * 1024;
#else // Linux
    float mbsize = 1024;
#endif
  

    std::cout << endl << "Maximum Resident Memory " << r_usage.ru_maxrss / mbsize << " MB" << std::endl;
    // clean up, close all files
    sigAdapter->close();
    for (auto& it : freqAdapters) {
      it.second->close();
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
    argparse::ArgumentParser args ("Merge signal events with up to four background sources.");
    
    args.add_argument("-i", "--signalFile")
      .default_value(std::string("root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1/pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root"))
      .help("Name of the HEPMC file with the signal events");
    
    args.add_argument("-sf", "--signalFreq")
      .default_value(0.0)
      .scan<'g', double>()
      .help("Signal frequency in kHz. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS rate to randomize.");
    
    args.add_argument("-S", "--signalSkip")
      .default_value(0)
    .scan<'i', int>()
    .help("Number of signals events to skip. Default is 0.");

    args.add_argument("-St", "--signalStatus")
      .default_value(0)
    .scan<'i', int>()
    .help("Apply shift on particle generatorStatus code for signal. Default is 0. ");

    args.add_argument("-b","--bgFile")
      .nargs(2,4)
      .append()
      .help("Tuple with name of the HEPMC file with background events, background frequency in kHz, number of background events to skip (default 0), shift on particle generatorStatus code (default 0).");
    
    args.add_argument("-o", "--outputFile")
      .default_value(std::string("bgmerged.hepmc3.tree.root"))
      .help("Specify the output file name. By default bgmerged.hepmc3.tree.root is used");

    args.add_argument("-r", "--rootFormat")
      .default_value(true)
      .implicit_value(true)
      .help("Use hepmc.root output format, default is true.");
	
    args.add_argument("-w", "--intWindow")
      .default_value(2000.0)
      .scan<'g', double>()
      .help("Length of the integration window in nanoseconds. Default is 2000.");
    
    args.add_argument("-N", "--nSlices")
      .default_value(10000)
      .scan<'i', int>()
      .help("Number of sampled time slices ('events'). Default is 10000. If set to -1, all events in the signal file will be used and background files cycled as needed.");
    
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
    signalSkip = args.get<int>("--signalSkip");
    signalStatus = args.get<int>("--signalStatus");
    backgroundFiles = parse_backgrounds(args.get<std::vector<std::string>>("--bgFile"));
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

    std::string statusMessage = "Shifting all particle status codes from this source by ";
    std::vector<int> statusList_stable, statusList_decay;

    std::cout << "Number of Slices:" << nSlices << endl;
    std::string freqTerm = signalFreq > 0 ? std::to_string(signalFreq) + " kHz" : "(one event per time slice)";
    std::string statusTerm = signalStatus > 0 ? statusMessage + std::to_string(signalStatus): "";
    if (signalStatus>0){
      statusList_stable.push_back(signalStatus+1);
      statusList_decay.push_back(signalStatus+2);
    }
    std::cout << "Signal events file and frequency:\n";
    std::cout << "\t- " << signalFile << "\t" << freqTerm << "\n" << statusTerm << "\n";
    
    std::cout << "\nBackground files and their respective frequencies:\n";

    for (const auto& bg : backgroundFiles) {
      if (!bg.file.empty()) {
        freqTerm = bg.frequency > 0 ? std::to_string(bg.frequency) + " kHz" : "(from weights)";
        statusTerm = bg.status > 0 ? statusMessage + std::to_string(bg.status) : "";
        std::cout << "\t- " << bg.file << "\t" << freqTerm << "\n" << statusTerm << "\n";
        if (bg.status>0){
          statusList_stable.push_back(bg.status+1);
          statusList_decay.push_back(bg.status+2);
        }
      }
    }
   	
    auto join = [](const std::vector<int>& vec) {
        return std::accumulate(vec.begin(), vec.end(), std::string(),
            [](const std::string& a, int b) {
                return a.empty() ? std::to_string(b) : a + " " + std::to_string(b);
            });
    };
    std::string stableStatuses = join(statusList_stable);
    std::string decayStatuses = join(statusList_decay);
    std::string message = "\n!!!Attention!!!\n To proceed the shifted particles statuses in DD4hep, please add the following options to ddsim:\n"
                          "--physics.alternativeStableStatuses=\"" + stableStatuses + 
                          "\"  --physics.alternativeDecayStatuses=\"" + decayStatuses + "\"\n";
    std::cout << message<<std::endl;           
  }
  
  // ---------------------------------------------------------------------------  
  void PrepData(const std::string& fileName, double freq, int skip=0, int baseStatus=0, bool signal=false) {
    if (fileName.empty()) return;

    cout << "Prepping " << fileName << endl;
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
    
    infoDict[fileName] = {0,0};

    if (signal) {
      sigAdapter = adapter;
      sigFreq = freq;
      sigStatus = baseStatus;
      sigAdapter->skip(skip);
      return;
    }

    // Now catch the weighted case
    if (freq <= 0) {
      std::cout << "Reading in all events from " << fileName << std::endl;
      std::vector<HepMC3::GenEvent> events;
      std::vector<double> weights;

      while(!adapter->failed()) {
	    HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
	    adapter->read_event(evt);

    	// remove events with 0 weight - note that this does change avgRate = <weight> (by a little)
	    if (double w=evt.weight() > 0){
	      events.push_back(evt);
	      weights.push_back(evt.weight());
	      }
      }
      adapter->close();
      
      double avgRate = 0.0;
      for ( auto w : weights ){ avgRate += w;}
      avgRate /= weights.size();
      avgRate *= 1e-9; // convert to 1/ns == GHz
      std::cout << "Average rate is " << avgRate << " GHz" << std::endl;

      std::vector<int> indices (weights.size());
      std::iota (std::begin(indices), std::end(indices), 0); // [ 0 , ... , N ] <- draw randomly from this
      
      // Replacing python's compact toPlace = self.rng.choice( a=events, size=nEvents, p=probs, replace=False )
      // is tricky. Possibly more elegant or faster versions exist,
      // https://stackoverflow.com/questions/42926209/equivalent-function-to-numpy-random-choice-in-c
      // we'll do it rather bluntly, since the need for this code should go away soon with new SR infrastructure
      // https://stackoverflow.com/questions/1761626/weighted-random-numbers
      // Normalizing is not necessary for this method 
      // for ( auto& w : weights ) {
      // 	w /= avgRate;
      // }
      std::piecewise_constant_distribution<> weightedDist(std::begin(indices),std::end(indices),
							  std::begin(weights));
      weightDict[fileName] = { std::make_tuple(events, weightedDist, avgRate) };

      return;
    }

    // Not signal and not weighted --> prepare frequency backgrounds
    adapter->skip(skip);
    freqAdapters[fileName] = adapter;
    freqs[fileName] = freq;
    baseStatuses[fileName] = baseStatus;
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
    // Generate a name for the output file
    // It's too simplistic for input with directories
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

    // More fine-grained info about current usage
#ifdef __MACH__
    task_basic_info_data_t info;
    mach_msg_type_number_t size = sizeof(info);
    kern_return_t kerr = task_info(mach_task_self(),
                                   TASK_BASIC_INFO,
                                   (task_info_t)&info,
                                   &size);

    long memory_usage = -1;
    if (kerr == KERN_SUCCESS) {
      memory_usage = info.resident_size  / 1024 / 1024;
    }
#else // Linux
    std::ifstream statm("/proc/self/statm");
    long size, resident, share, text, lib, data, dt;
    statm >> size >> resident >> share >> text >> lib >> data >> dt;
    statm.close();

    long page_size = sysconf(_SC_PAGESIZE);  // in case x86-64 is configured to use 2MB pages
    long memory_usage = resident * page_size  / 1024 / 1024 ;    
#endif
  
    
    std::cout << "Working on slice " << i + 1 << std::endl;
    std::cout << "Current memory usage: " << memory_usage << " MB" << std::endl;

  }
  // ---------------------------------------------------------------------------

  std::unique_ptr<HepMC3::GenEvent> mergeSlice(int i) {
    auto hepSlice = std::make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
    
    addFreqEvents(signalFile, sigAdapter, sigFreq, hepSlice, signalStatus, true);
    
    for (const auto& freqBgs : freqAdapters) {
      auto fileName=freqBgs.first;
      addFreqEvents(fileName, freqAdapters[fileName], freqs[fileName], hepSlice, baseStatuses[fileName], false);
    }
    
    for (const auto& fileName : weightDict) {
      addWeightedEvents(fileName.first, hepSlice, baseStatuses[fileName.first]);
    }

    return hepSlice;
  };

  // ---------------------------------------------------------------------------

  void addFreqEvents(std::string fileName, std::shared_ptr<HepMC3::Reader>& adapter, const double freq,
		     std::unique_ptr<HepMC3::GenEvent>& hepSlice, int baseStatus = 0, bool signal = false) {

    // First, create a timeline
    // Signals can be different
    std::vector<double> timeline;

    std::uniform_real_distribution<> uni(0, intWindow);
    if (freq == 0){
      if (!signal) {
        std::cerr << "frequency can't be 0 for background files" << std::endl;
        exit(1);
      }
      // exactly one signal event, at an arbitrary point
      timeline.push_back(uni(rng));
    } else {
      // Generate poisson-distributed times to place events
      timeline = poissonTimes(freq, intWindow);
    }
    
    if ( verbose) std::cout << "Placing " << timeline.size() << " events from " << fileName << std::endl;

    if (timeline.empty()) return;
    long particleCount = 0;

    // Insert events at all specified locations
    for (double time : timeline) {
      if(adapter->failed()) {
      try{
        if (signal) { // Exhausted signal events
	        throw std::ifstream::failure("EOF");
	      } else { // background file reached its end, reset to the start
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

      if (squashTime) time = 0;
      particleCount += insertHepmcEvent( inevt, hepSlice, time, baseStatus, signal);
    }

    infoDict[fileName].eventCount += timeline.size();
    infoDict[fileName].particleCount += particleCount;


    return;
  }

  // ---------------------------------------------------------------------------

  void addWeightedEvents(std::string fileName, std::unique_ptr<HepMC3::GenEvent>& hepSlice, int baseStatus=0, bool signal = false) {
    auto& [events, weightedDist, avgRate ] = weightDict[fileName];

    // How many events? Assume Poisson distribution
    int nEvents;
    std::poisson_distribution<> d( intWindow * avgRate );

    // Small SR files may not have enough photons (example or test files). Could use them all or reroll
    // Choosing the latter, neither is physical
    while (true) {
      nEvents = d(rng);
      if (nEvents > events.size()) {
	      std::cout << "WARNING: Trying to place " << nEvents << " events from " << fileName
       	      	  << " but the file doesn't have enough. Rerolling, but this is not physical." << std::endl;
      	continue;
      }
      break;
    }

    if (verbose) std::cout << "Placing " << nEvents << " events from " << fileName << std::endl;
    
    // Get randomized event indices
    // Note: Could change to drawing without replacing ( if ( not in toPLace) ...) , not worth the effort
    std::vector<HepMC3::GenEvent> toPlace(nEvents);
    for ( auto& e : toPlace ){
      auto i = static_cast<int> (weightedDist(rng));
      e = events.at(i);
    }
    
    // Place at random times
    std::vector<double> timeline;
    std::uniform_real_distribution<> uni(0, intWindow);
    long particleCount = 0;
    if (!squashTime) {
      for ( auto& e : toPlace ){
	      double time = squashTime ? 0 : uni(rng);
        particleCount += insertHepmcEvent( e, hepSlice, time, baseStatus, signal);
      }
    }

    infoDict[fileName].eventCount += nEvents;
    infoDict[fileName].particleCount += particleCount;

    return;
}

  // ---------------------------------------------------------------------------
  long insertHepmcEvent( const HepMC3::GenEvent& inevt,
			 std::unique_ptr<HepMC3::GenEvent>& hepSlice, double time=0, int baseStatus=0, bool signal = false) {
    // Unit conversion
    double timeHepmc = c_light * time;
    
    std::vector<HepMC3::GenParticlePtr> particles;
    std::vector<HepMC3::GenVertexPtr> vertices;

    // Stores the vertices of the event inside a vertex container. These vertices are in increasing order
    // so we can index them with [abs(vertex_id)-1]
    for (auto& vertex : inevt.vertices()) {
      HepMC3::FourVector position = vertex->position();
      position.set_t(position.t() + timeHepmc);
      auto v1 = std::make_shared<HepMC3::GenVertex>(position);
      vertices.push_back(v1);
    }
      
    // copies the particles and attaches them to their corresponding vertices
    long finalParticleCount = 0;
    for (auto& particle : inevt.particles()) {
      HepMC3::FourVector momentum = particle->momentum();
      int status = particle->status();
      if (status == 1 ) finalParticleCount++;
      int pid = particle->pid();
      status += baseStatus;
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
	      vertices.at(abs(end_vertex) - 1)->add_particle_in(p1);	
      }
    }

    // Adds the vertices with the attached particles to the event
    for (auto& vertex : vertices) {
      hepSlice->add_vertex(vertex);
    }
    
    return finalParticleCount;
  }

  // ---------------------------------------------------------------------------

  std::vector<double> poissonTimes(double mu, double endTime) {
    std::exponential_distribution<> exp(mu);
    
    double t = 0;
    std::vector<double> ret;
    while (true) {
      double delt = exp(rng)*1e6;
      // cout << delt <<endl;
      t += delt;
      if (t >= endTime) {
	break;
      }
      ret.push_back(t);
    }
    return ret;
}
  // ---------------------------------------------------------------------------

  
  // private:
  std::mt19937 rng;
  string signalFile;
  double signalFreq;
  int signalSkip;
  int signalStatus;
  std::vector<BackgroundConfig> backgroundFiles;  
  string outputFile;
  string outputFileName;
  bool rootFormat;
  double intWindow;
  int nSlices; // should be long, but argparse cannot read that
  bool squashTime;
  int rngSeed;  // should be unsigned, but argparse cannot read that
  bool verbose;
  
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
