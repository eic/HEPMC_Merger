# EPIC HEPMC MERGER
## Original author: Benjamen Sterwerf, UC Berkeley (bsterwerf@berkeley.edu)
## Maintainer of the C++ implementation: Kolja Kauder (kkauder@bn=l.gov)
This code is used to combine HEPMC files. You can run this code as:

```bash
    ./SignalBackgroundMerger --signalFile dis.hepmc --signalFreq 0 \
            --bg1File hgas.hepmc --bg1Freq 31.9 \
            --bg2File egas.hepmc --bg2Freq 3177.25 \
            --bg3File sr.hepmc
```

Without input flags, `./SignalBackgroundMerger` will create 10k time frames of 2 nanosecond legths streaming in input files from XrootD and uses frequencies from https://wiki.bnl.gov/EPIC/index.php?title=Background

For more information run:

```bash
./SignalBackgroundMerger --help
```

 ## Software:
- Standalone framework independent C++ code
- Supports all HepMC input and output formats such as `hepmc3.tree.root`
- Supports streaming I/O via XrootD
- Supports weighted inputs such as coming from legacy SR data (though this is currenbtly not recommended. Instead use newer SR input once available)
- A legacy python version exists in the python subdirectory. It is known to treat weighted input slightly wrong.

## Installation
Prerequisites: ROOT, HepMC
Includes the header Argument Parser library argparse by Pranav Srinivas Kumar under the MIT license

```bash
git clone https://github.com/eic/HEPMC_Merger.git
mkdir -p HEPMC_Merger/build
cd HEPMC_Merger/build
cmake ..
make
```


<!-- Preinstalled versions are available in [ATHENA containers](https://doc.athena-eic.org/en/latest/overview/containers.html).  -->


## Logic:
For a given time frame, the number of signal and background events to place is drawn from a Poisson distribution corresponding to its frequency. Events are then distributed uniformly within the frame. 

If the signal frequency is set to 0, exactly one signal event will be placed in the slice.

If a background frequency is set to 0, events are instead drawn from a weighted distribution that must be provided by the weights in the input file. This functionality should not be needed much longer since Synchrotron radiation files will in the future be provided in the same format as all other input files.

Note: When the end of a background file is reached, it cycles back to the beginning. 


