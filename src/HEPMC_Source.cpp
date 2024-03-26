#include "HEPMC_Source.h"

// ---------------------------------------------------------------------------
//Constructor
// ---------------------------------------------------------------------------
HEPMC_Source::HEPMC_Source(std::string fileName, double freq, int sourceNo): m_fileName(fileName), m_freq(freq), m_sourceNo(sourceNo){
    
    try {
        adapter = HepMC3::deduce_reader(m_fileName);
        if (!adapter) {
            throw std::runtime_error("Failed to open file");
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Opening " << fileName << " failed: " << e.what() << std::endl;
        exit(1);
    }

    // Check if the file is weighted
    if (m_freq < 0) {
        isWeighted = true;
        SetupWeights();
    } else {
        isWeighted = false;
    }        

}

// ---------------------------------------------------------------------------
// Setup the a weighted data source
// ---------------------------------------------------------------------------
void HEPMC_Source::SetupWeights() {
    // Set up the weights
    // Read the file and store the events
    std::vector<double> weights;
    double weightSum = 0;
    while (!adapter->failed()) {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
        if(evt.weight() <= 0) continue;
        
        eventList.push_back(evt);
        weights.push_back(evt.weight());
        weightSum += evt.weight();
    }
    
    m_freq = (weightSum/eventList.size())*1e-9; // Calculate average frequency and convert to GHz

    // Create the distribution
    weightedDist(weights);
}

// ---------------------------------------------------------------------------
// Generate a timeline of event times
// ---------------------------------------------------------------------------
std::vector<double> HEPMC_Source::GenerateSampleTimes(double intWindow, double bunchSpacing) {
    
    int nEvents = 0;
    if (m_freq == 0){
        nEvents = 1;
    }
    else {
        std::poisson_distribution<> d(intWindow * m_freq);
        nEvents = d(rng);
    }

    std::vector<double> timeline;
    
    // Set up the distribution
    if(bunchCorrelated) {
        int nBunches = int(intWindow/bunchSpacing);
        std::uniform_int_distribution<>  dist(0, nBunches);
    }
    else {
        std::uniform_real_distribution<> dist(0, intWindow);
    }

    // Generate the timeline
    for(int i = 0; i < nEvents; i++) {
        double time = dist(rng);
        if(bunchCorrelated) {
            time *= bunchSpacing;
        }
        timeline.push_back(time);
    }

    return timeline;
}