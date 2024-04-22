#pragma once

#include <string>
#include <vector>
#include <random>
#include <memory>
#include <HepMC3/Reader.h>
#include <HepMC3/Units.h>
#include <HepMC3/GenEvent.h>

class HEPMC_Source {

public:
    HEPMC_Source(std::string fileName, double freq, int sourceNo, int beamCorrelation);

    void close() { adapter->close(); }

    inline std::string getFileName() const { return m_fileName; }
    inline double      getFreq()     const { return m_freq; }
    inline int         getSourceNo() const { return m_sourceNo; }
    inline bool        getIsWeighted()  const { return isWeighted; }
    inline bool        getIsBunchCorrelated() const { return m_bunchCorrelated; }
    inline int         getBeamCorrelation() const { return m_beamCorrelation; }

    void SetupWeights();
    std::vector<double> GenerateSampleTimes(double intWindow, double bunchSpacing, std::mt19937& rng);
    HepMC3::GenEvent getNextEvent(std::mt19937& rng);

private:
    std::shared_ptr<HepMC3::Reader> adapter;

    std::string m_fileName;
    double m_freq; // frequency of the source kHz
    int    m_sourceNo; // source reference number

    bool m_bunchCorrelated = false;

    // Some sort of matrix to do a transformation based on vertex 4-vector?
    // Maybe just a time correction to correlate with z position? Needs energy and mass of particle
    int m_beamCorrelation = 1; // 0 = none, 1 = electron correlation, 2 = ion correlation

    // Some sort of gaussian distribution to smear the time out by
    double gausWidth; // ns

    // Variables for calculating the weight
    bool isWeighted = false;
    std::vector<HepMC3::GenEvent> eventList;
    std::discrete_distribution<> weightedDist; // Random number distributor to select event from eventList based on its weight 
    
    
};