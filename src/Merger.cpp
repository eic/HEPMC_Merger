#include "Merger.h"
#include <sys/resource.h>

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
Merger::Merger(std::string outputFile, bool rootFormat, double intWindow, int rngSeed, bool verbose, bool squashTime, double bunchSpacing) :
    outputFileName(outputFile), m_intWindow(intWindow), m_verbose(verbose), m_squashTime(squashTime), m_bunchSpacing(bunchSpacing) {

    auto t0 = std::chrono::high_resolution_clock::now();

    // initialize rng
    rng.seed( rngSeed );
    
    // Setup output file 
    if (outputFileName.empty()) {
        if(rootFormat) outputFileName = "bgmerged.root";
        else outputFileName = "bgmerged.hepmc";
    }

    // Open output file
    if (rootFormat) {
        outFile = std::make_shared<HepMC3::WriterRootTree>(outputFileName);
    } else {
        outFile = std::make_shared<HepMC3::WriterAscii>(outputFileName);
    }

    std::cout << "Writing to " << outputFileName << std::endl;

    // // DEBUG
    // f = new TFile("f.root","RECREATE");
    // p = new TH1I("p","poisson",20,-1,19);
    // e = new TH1I("e","exponential",20,-1,19);
    // // /DEBUG

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initiation time: " << std::round(std::chrono::duration<double, std::chrono::seconds::period>(t1 - t0).count()) << " sec" << std::endl;
    std::cout << "\n==================================================================\n" << std::endl;

}


// --------------------------------------------------------------------------- 
// Perform merge
// ---------------------------------------------------------------------------
void Merger::merge(int nSlices) {

    auto t1 = std::chrono::high_resolution_clock::now();

    // Print banner
    printBanner(nSlices);

    // Slice loop
    int i = 0;
    for (i = 0; i<nSlices; ++i ) {
        if (i % 100 == 0 || m_verbose ) squawk(i);
        auto hepSlice = mergeSlice(i);
        if (!hepSlice) {
        std::cout << "Exhausted signal source." << std::endl;
        break;
        }
        hepSlice->set_event_number(i);
        outFile->write_event(*hepSlice);
    }
    std::cout << "Finished all requested slices." << std::endl;

    int slicesDone = i;
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Slice loop time: " << std::round(std::chrono::duration<double, std::chrono::minutes::period>(t2 - t1).count()) << " min" << std::endl;
    std::cout << " -- " << std::round(std::chrono::duration<double>(t2 - t1).count() / i) << " sec / slice" << std::endl;

    // Clean up, close all files
    for(auto source : sources) {
        source.close();
    }

    outFile->close();

}

// ---------------------------------------------------------------------------
// Merge single slice
// ---------------------------------------------------------------------------
std::unique_ptr<HepMC3::GenEvent> Merger::mergeSlice(int i) {

    auto hepSlice = std::make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);

    for (auto source : sources) {
        addEvents(source, hepSlice);
    }
    

    // addFreqEvents(signalFile, sigAdapter, sigFreq, hepSlice, true);

    // for (const auto& freqBgs : freqAdapters) {
    //     auto fileName=freqBgs.first;
    //     addFreqEvents(fileName, freqAdapters[fileName], freqs[fileName], hepSlice);
    // }

    // for (const auto& fileName : weightDict) {
    //     addWeightedEvents(fileName.first, hepSlice);
    // }

    return hepSlice;
}

// ---------------------------------------------------------------------------
// Print banner  
// ---------------------------------------------------------------------------
void Merger::printBanner(int nSlices){
    std::cout << "==================================================================" << std::endl;
    std::cout << "=== EPIC HEPMC MERGER ===" << std::endl;
    std::cout << "authors: Benjamen Sterwerf* (bsterwerf@berkeley.edu), Kolja Kauder** (kkauder@bnl.gov), Reynier Cruz-Torres***" << std::endl;
    std::cout << "* University of California, Berkeley" << std::endl;
    std::cout << "** Brookhaven National Laboratory" << std::endl;
    std::cout << "*** formerly Lawrence Berkeley National Laboratory" << std::endl;
    std::cout << "\nFor more information, run \n./signal_background_merger --help" << std::endl;

    std::cout << "Number of Slices:" << nSlices << std::endl;
    std::string freqTerm = sources[0].getFreq() > 0 ? std::to_string(sources[0].getFreq()) + " GHz" : "(one event per time slice)";
    std::cout << "Signal events file and frequency:\n";
    std::cout << "\t- " << sources[0].getFileName() << "\t" << freqTerm << "\n";

    std::cout << "\nBackground files and their respective frequencies:\n";
    for(auto it = sources.begin() + 1; it != sources.end(); ++it) {
        auto& source = *it;
        freqTerm = source.getFreq() > 0 ? std::to_string(source.getFreq()) + " GHz" : "(one event per time slice)";
        std::cout << "\t- " << source.getFileName() << "\t" << freqTerm << "\n";
    }
}

// ---------------------------------------------------------------------------
// Add data source
// ---------------------------------------------------------------------------
void Merger::addSource(const std::string fileName, double freq, int sourceNo) {
    if (fileName.empty()) return;
    sources.push_back(HEPMC_Source(fileName, freq, sourceNo));
}

// ---------------------------------------------------------------------------
// squark
// ---------------------------------------------------------------------------
void Merger::squawk(int i) {
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
// Add events
// --------------------------------------------------------------------------- 
void Merger::addEvents(HEPMC_Source& source, std::unique_ptr<HepMC3::GenEvent>& hepSlice) {
    
    // First, create a timeline
    std::vector<double> timeline = source.GenerateSampleTimes(m_intWindow,  m_bunchSpacing, rng);
    
    if (m_verbose) std::cout << "Placing " << timeline.size() << " events from " << source.getFileName() << std::endl;
    
    // Loop over the timeline
    for (auto time : timeline) {
        // Read the event
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);

        evt = source.getNextEvent(rng);

        insertHepmcEvent( evt, hepSlice, time);
    }
}

    
// ---------------------------------------------------------------------------
void Merger::insertHepmcEvent( const HepMC3::GenEvent& inevt,
            std::unique_ptr<HepMC3::GenEvent>& hepSlice, double time) {
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