# Signal and Background Merger Code

this code is used to combine two HEPMC files.

## Dependencies

We need to install the ```pyhepmc``` module (e.g. ```pip install pyhepmc```).

## Variables

1. one_background_particle
	- default=True
	- Set to True if all the background vertices in the file are to be connected temporally and spatially and False for backgrounds with many particles in the event like SR (synchrotron radiation)
2. shifter
	- default=True
	- depending on how integration frames are decided, the signal event will either be positioned at time 0 (False) or set to a random time within the integration window
3. Int_Window
	- default=2000
	- length of the integration window in nanoseconds
4. Signal_File
	- default='Test_DIS_event.hepmc'
	- Name of the HEPMC file with the signal events
5. Background_File
	- default= 'Test_Back_event.hepmc'
	- Name of the HEPMC file with the background events

