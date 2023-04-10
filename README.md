# EPIC HEPMC MERGER
## author: Benjamen Sterwerf, UC Berkeley (bsterwerf@berkeley.edu)

This code is used to combine HEPMC files. You can run this code as:

```bash
python signal_background_merger.py --Signal_File signal.hepmc --Background_Files background_1.hepmc background_2.hepmc
```

For more information run:

```bash
python signal_background_merger.py --help
```

## Dependencies

We need to install the ```pyhepmc``` module (e.g. ```pip install pyhepmc```).

## Variables

1. one_background_particle
	- default=True
	- Set to True if all the background vertices in the file are to be connected temporally and spatially and False for backgrounds with many particles in the event like SR (synchrotron radiation)
2. Int_Window
	- default=0
	- length of the integration window in nanoseconds. If set to a value > 0, the signal event will either be set to a random time within the integration window 
3. Signal_File
	- default = 'dummy_signal.hepmc'
	- Name of the HEPMC file with signal events
4. Background_Files
	- default = ['dummy_background_1.hepmc','dummy_background_2.hepmc']
	- Names of the HEPMC files with background events

