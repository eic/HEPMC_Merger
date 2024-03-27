#include <HepMC3/ReaderFactory.h>

// Horrible hack until hepmc3 version is updated to be less nasty.
std::shared_ptr<HepMC3::Reader> deduce_reader_wrapper(const std::string& filename) {
    return HepMC3::deduce_reader(filename);
}