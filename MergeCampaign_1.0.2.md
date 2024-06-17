Used: https://github.com/eic/HEPMC_Merger/releases/tag/v1.0.2

Input:
```
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1/pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root .
```

Electron gas:
```
xrdcp root://dtn-eic.jlab.org//work/eic2/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/electron/GETaLM1.0.0-1.0/10GeV/GETaLM1.0.0-1.0_ElectronBeamGas_10GeV_foam_emin10keV_run001.hepmc3.tree.root .
```

Proton gas:
```
for f in `xrdfs root://dtn-eic.jlab.org ls /work/eic2/EPIC/EVGEN/BACKGROUNDS/BEAMGAS/proton/pythia8.306-1.0/100GeV | grep 'pythia8.306-1.0_ProtonBeamGas_100GeV_run00.*.hepmc3.tree.root'`
do
xrdcp root://dtn-eic.jlab.org/$f .
done
hadd pythia8.306-1.0_ProtonBeamGas_100GeV_runs00.hepmc3.tree.root pythia8.306-1.0_ProtonBeamGas_100GeV_run00*
```

Merge:
```
./build/SignalBackgroundMerger -N 25000 \
-i pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root -sf 184 \
-bg1 pythia8.306-1.0_ProtonBeamGas_100GeV_runs00.hepmc3.tree.root -bf1 31.9 \
-bg2 GETaLM1.0.0-1.0_ElectronBeamGas_10GeV_foam_emin10keV_run001.hepmc3.tree.root -bf2 3177.25 \
-bg3 "" -r \
-w 2000.0 --rngSeed 42 \
-o HEPMC_merger-1.0.2_bgmerged_25k_RealisticSignalPerFrame_MinBias_pythia6_10x100_egas_bgas.hepmc3.tree.root
```

```
./build/SignalBackgroundMerger -N 10000 \
-i pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root -sf 0 \
-bg1 pythia8.306-1.0_ProtonBeamGas_100GeV_runs00.hepmc3.tree.root -bf1 31.9 \
-bg2 GETaLM1.0.0-1.0_ElectronBeamGas_10GeV_foam_emin10keV_run001.hepmc3.tree.root -bf2 3177.25 \
-bg3 "" -r \
-w 2000.0 --rngSeed 42 \
-o HEPMC_merger-1.0.2_bgmerged_10k_1SignalPerFrame_MinBias_pythia6_10x100_egas_bgas.hepmc3.tree.root
```

