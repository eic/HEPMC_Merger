#include Merger.h


// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
Merger(std::string outputFile, bool rootFormat, double intWindow, int rngSeed, bool verbose, bool squashTime) {

    auto t0 = std::chrono::high_resolution_clock::now();

    // initialize rng
    rng.seed( rngSeed );

    // Setup output file 
    if (outputFile != "" ) {
        outputFileName = outputFile;
    } else {
        outputFileName = nameGen();
    }

    // Open output file
    if (rootFormat) {
        outFile = std::make_shared<HepMC3::WriterRootTree>(outputFileName);
    } else {
        outFile = std::make_shared<HepMC3::WriterAscii>(outputFileName);
    }

    cout << "Writing to " << outputFileName << endl;

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
    banner(nSlices);

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
        outFile->write_event(*hepSlice);
    }
    std::cout << "Finished all requested slices." << std::endl;

    int slicesDone = i;
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Slice loop time: " << std::round(std::chrono::duration<double, std::chrono::minutes::period>(t2 - t1).count()) << " min" << std::endl;
    std::cout << " -- " << std::round(std::chrono::duration<double>(t2 - t1).count() / i) << " sec / slice" << std::endl;

    // clean up, close all files
    sigAdapter->close();
    for (auto& it : freqAdapters) {
        it.second->close();
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
    

    addFreqEvents(signalFile, sigAdapter, sigFreq, hepSlice, true);

    for (const auto& freqBgs : freqAdapters) {
        auto fileName=freqBgs.first;
        addFreqEvents(fileName, freqAdapters[fileName], freqs[fileName], hepSlice);
    }

    for (const auto& fileName : weightDict) {
        addWeightedEvents(fileName.first, hepSlice);
    }

    return hepSlice;
};


// ---------------------------------------------------------------------------
// Print banner  
// ---------------------------------------------------------------------------
void Merger::printBanner(){
    std::cout << "==================================================================" << std::endl;
    std::cout << "=== EPIC HEPMC MERGER ===" << std::endl;
    std::cout << "authors: Benjamen Sterwerf* (bsterwerf@berkeley.edu), Kolja Kauder** (kkauder@bnl.gov), Reynier Cruz-Torres***" << std::endl;
    std::cout << "* University of California, Berkeley" << std::endl;
    std::cout << "** Brookhaven National Laboratory" << std::endl;
    std::cout << "*** formerly Lawrence Berkeley National Laboratory" << std::endl;
    std::cout << "\nFor more information, run \n./signal_background_merger --help" << std::endl;

    std::cout << "Number of Slices:" << nSlices << endl;
    std::string freqTerm = sources[0].getFreq() > 0 ? std::to_string(sources[0].getFreq()) + " kHz" : "(one event per time slice)";
    std::cout << "Signal events file and frequency:\n";
    std::cout << "\t- " << sources[0].getFileName() << "\t" << freqTerm << "\n";

    std::cout << "\nBackground files and their respective frequencies:\n";
    for(auto it = sources.begin() + 1; it != sources.end(); ++it) {
        auto& source = *it;
        freqTerm = source.getFreq() > 0 ? std::to_string(source.getFreq()) + " kHz" : "(one event per time slice)";
        std::cout << "\t- " << source.getFileName() << "\t" << freqTerm << " kHz" << "\n";
    }
}


// ---------------------------------------------------------------------------
// Add data source
// ---------------------------------------------------------------------------
void Merger::addSource(const std::string fileName, double freq) {
    if (fileName.empty()) return;
    sources.push_back(HEPMC_Source(fileName, freq));
}