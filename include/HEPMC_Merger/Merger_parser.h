#pragma once
#include "argparse.hpp"
// =============================================================
// digest command line arguments
// =============================================================
argparse::ArgumentParser* digestArgs(int argc, char* argv[]) {
    // Handle the command line tedium
    // ArgumentParser is meant to be used in a single function.
    // ArgumentParser internally uses std::string_views,
    // references, iterators, etc.
    // Many of these elements become invalidated after a copy or move.
    
  auto args = new argparse::ArgumentParser("Merge signal events with up to three background sources.");
    
  args->add_argument("-i", "--signalFile")
    .default_value(std::string("small_ep_noradcor.10x100_q2_10_100_run001.hepmc"))
    .help("Name of the HEPMC file with the signal events");

  args->add_argument("-sf", "--signalFreq")
    .default_value(0.0)
    .scan<'g', double>()
    .help("Signal frequency in kHz. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS rate to randomize.");
  
  args->add_argument("-bg1", "--bg1File")
    .default_value(std::string("small_hgas_100GeV_HiAc_25mrad.Asciiv3.hepmc"))
    .help("Name of the first HEPMC file with background events");
  
  args->add_argument("-bf1", "--bg1Freq")
    .default_value(342.8)
    .scan<'g', double>()
    .help("First background frequency in kHz. Default is the estimated hadron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.");
  
  args->add_argument("-bg2", "--bg2File")
    .default_value(std::string("small_beam_gas_ep_10GeV_foam_emin10keV_vtx.hepmc"))
    .help("Name of the second HEPMC file with background events");
  
  args->add_argument("-bf2", "--bg2Freq")
    .default_value(3177.25)
    .scan<'g', double>()
    .help("Second background frequency in kHz. Default is the estimated electron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.");
  
  args->add_argument("-bg3", "--bg3File")
    .default_value(std::string("small_SR_single_2.5A_10GeV.hepmc"))
    .help("Name of the third HEPMC file with background events");
  
  args->add_argument("-bf3", "--bg3Freq")
    .default_value(0.0)
    .scan<'g', double>()
    .help("Third background frequency in kHz. Default is 0 to use the weights in the corresponding input file. Set to a value >0 to specify a frequency instead.");
  
  args->add_argument("-o", "--outputFile")
    .default_value(std::string(""))
    .help("Specify the output file name. By default it will be auto-generated.");

  args->add_argument("-r", "--rootFormat")
    .default_value(false)
    .implicit_value(true)
    .help("Use hepmc.root output format, default is false.");

  args->add_argument("-w", "--intWindow")
    .default_value(2000.0)
    .scan<'g', double>()
    .help("Length of the integration window in nanoseconds. Default is 2000.");

  args->add_argument("-bs", "--bunchSpacing")
    .default_value(10.0)
    .scan<'g', double>()
    .help("Time between bunch crossings. Default is 10.");
  
  args->add_argument("-N", "--nSlices")
    .default_value(-1)
    .scan<'i', int>()
    .help("Number of sampled time slices ('events'). Default is -1. If set to -1, all events in the signal file will be used and background files cycled as needed.");
  
  args->add_argument("--squashTime")
    .default_value(false)
    .implicit_value(true)
    .help("Integration is performed but no time information is associated to vertices.");
  
  args->add_argument("--rngSeed")
    .default_value(0)
    .action([](const std::string& value) { return std::stoi(value); })
    .help("Random seed, default is None");
  
  args->add_argument("-v", "--verbose")
    .default_value(false)
    .implicit_value(true)
    .help("Display details for every slice.");
  
  try {
    args->parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cout << err.what() << std::endl;
    std::cout << args;
    exit(0);
  }
     
  return args;

}