
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
  Merger(std::string outputFile, bool rootFormat, double intWindow, int rngSeed, bool verbose, bool squashTime):
        m_intWindow(intWindow), m_verbose(verbose), m_squashTime(squashTime);

  void merge(int nSlices); 
  void printBanner(int nSlices);
  
  // ---------------------------------------------------------------------------
  // Add data source
  // ---------------------------------------------------------------------------
  void addSource(const std::string fileName, double freq);

  // ---------------------------------------------------------------------------  
  void PrepData(const std::string& fileName, double freq, bool signal=false) {
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
    
    if (signal) {
      sigAdapter = adapter;
      sigFreq = freq;
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
	// HepMC3::Print::listing(evt);
	// HepMC3::Print::content(evt);

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
    freqAdapters[fileName] = adapter;
    freqs[fileName] = freq;
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

  std::unique_ptr<HepMC3::GenEvent> mergeSlice(int i);

  // ---------------------------------------------------------------------------

  void addFreqEvents(std::string fileName, std::shared_ptr<HepMC3::Reader>& adapter, const double freq,
		     std::unique_ptr<HepMC3::GenEvent>& hepSlice, bool signal = false) {

    // First, create a timeline
    // Signals can be different
    std::vector<double> timeline;
    std::uniform_real_distribution<> uni(0, intWindow);
    if (freq == 0) {
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

    // // DEBUG
    // e->Fill(timeline.size());
    // std::poisson_distribution<> d( freq*1e-6 * intWindow );
    // auto np = d(rng);
    // p->Fill(np);
    // // /DEBUG
    
    
    if (timeline.empty()) return;

    // Insert events at all specified locations
    for (double time : timeline) {
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

      if (squashTime) time = 0;
      insertHepmcEvent( inevt, hepSlice, time);
      
    }
    return;
  }

  // ---------------------------------------------------------------------------

  void addWeightedEvents(std::string fileName, std::unique_ptr<HepMC3::GenEvent>& hepSlice, bool signal = false) {
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
    if (!squashTime) {
      for ( auto& e : toPlace ){
	double time = squashTime ? 0 : uni(rng);
	insertHepmcEvent( e, hepSlice, time);
      }
    }

    return;
}

  // ---------------------------------------------------------------------------
  void insertHepmcEvent( const HepMC3::GenEvent& inevt,
			 std::unique_ptr<HepMC3::GenEvent>& hepSlice, double time=0) {
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
	vertices.at(abs(end_vertex) - 1)->add_particle_in(p1);	
      }
    }

    // Adds the vertices with the attached particles to the event
    for (auto& vertex : vertices) {
      hepSlice->add_vertex(vertex);
    }
    
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

  
  private:
    std::shared_ptr<HepMC3::Writer> outFile;
    std::mt19937 rng;
    
    bool   rootFormat;
    double m_intWindow;
    bool   m_squashTime;
    bool   m_verbose;
    
    const double c_light = 299.792458; // speed of light = 299.792458 mm/ns to get mm

    std::vector<HEPMC_Source> sources;

    std::map<std::string,
        std::tuple<std::vector<HepMC3::GenEvent>,
                std::piecewise_constant_distribution<>,
                double>
        > weightDict;
  // // DEBUG  
  // TFile *f;
  // TH1I* p;
  // TH1I* e;

  
};
