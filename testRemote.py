import subprocess

subprocess.run(["./build/SignalBackgroundMerger",
                "--nSlices", "10",
                "--intWindow", "2000.0",
                "--signalFile", "root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1/pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run5.ab.hepmc3.tree.root", 
                "--signalFreq", "0",
                "--bg1File", "root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/100GeV/pythia8.306-1.0_ProtonBeamGas_100GeV_run091.hepmc3.tree.root",
                "--bg1Freq", "18520",
                "--bg2File", "root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/electron/beam_gas_ep_10GeV_foam_emin10keV_10Mevt_vtx_cs_info.hepmc3.tree.root", 
                "--bg2Freq", "18520",
                "--bg3File", "root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/BACKGROUNDS/SYNRAD/SR_single_1.8M_2.5A_10GeV.hepmc3.tree.root",
                "--outputFile", "mergedC.root",
                "--rootFormat",
                "-v"])
