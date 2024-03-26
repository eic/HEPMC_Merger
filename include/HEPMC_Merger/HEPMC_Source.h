
class HEPMC_Source {

public:
    HEPMC_Source(std::string fileName, double freq, int sourceNo): m_fileName(fileName), m_freq(freq), m_sourceNo(sourceNo){
        
        try {
            adapter = HepMC3::deduce_reader(m_fileName);
            if (!adapter) {
                throw std::runtime_error("Failed to open file");
            }
        } catch (const std::runtime_error& e) {
            std::cerr << "Opening " << fileName << " failed: " << e.what() << std::endl;
            exit(1);
        }
        
    }

    std::string getFileName() const {
        return m_fileName;
    }

    double getFreq() const {
        return m_freq;
    }

    int getSourceNo() const {
        return m_sourceNo;
    }

private:
    std::shared_ptr<HepMC3::Reader> adapter;
    string m_fileName;
    double m_freq; // frequency of the source kHz
    int    m_sourceNo; // source reference number

};