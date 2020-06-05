# phi-meson Spin Alignment Code

### TreeProduction:
> - reconstrut event plane => ZDC-SMD (1st) and TPC (2nd)
> - reconstrut phi-meson => Same Event and Mixed Event
> - save event plane and reconstruted phi-meson into TTree

### FillSpinAlignment:
> - read in phi-meson TTree
> - boost K+ back into phi-meson rest frame
> - correlate phi-meson and event plane with K+ (cos(theta\*))
> - save histogram

### CalSpinAlignment:
> - subtract background from Mixed Event
> - extract raw spin alignment signal
> - apply TPC efficiency correction
> - apply MC event plane resolution correction
> - calculate systematic error

### Embedding
> - read in STAR embedding TTree
> - find MC tracks and corresponding RC tracks
> - Save all relavent information to TNtuple

### RcPhiEffCorr
> - read in TNtuple from Embedding process
> - apply same cut as picoDst data analysis
> - extract TPC efficiency for K+/K-
> - generate phi-meson through MC and let it decay with PHYTIA
> - apply K+/K- TPC effciency to decay daugthers and reconstruct phi-meson
> - extract TPC efficiency for phi-meson

### McPhiResCorr
> - generate phi-meson through MC with STAR published spectra and v2
> - let phi-meson decay with PHYTIA and boost K+ back to phi-meson rest frame
> - correlate phi-meson with fixed event plane (Psi = 0)
> - smear event plane with measured event plane resolution
> - correlate rho_00^phy with rho_00^obs to extract event plane resolution correction factor

### Utility
> - constant used in `VecMesonSpinAlignment`
> - functions used in `VecMesonSpinAlignment`
> - custom defined fype used in `VecMesonSpinAlignment`

### PlotMacro
> - macros to plot QA and final figures for rho_00 vs energy with sysErrors
> - macros to plot figures in Analysis Note and Paper Draft

### figures
> - place to save QA and final plots

### Procedure to run the code
- This package is based on LBNL PicoDst which is not available anymore.
- `TreeProduction` was used to produce the event plane and phi-meson TTree.
- The phi-meson TTree of Same Event and Mixed Events have been saved for the analysis.
- All results can be re-produced from phi-meson TTree.
- The saved TTree can be found at: 
  - /star/data01/pwg/sunxuhit/AuAu200GeV/SpinAlignment/Phi/Forest
  - /star/data01/pwg/sunxuhit/AuAu62GeV/SpinAlignment/Phi/Forest
  - /star/data01/pwg/sunxuhit/AuAu39GeV/SpinAlignment/Phi/Forest
  - /star/data01/pwg/sunxuhit/AuAu27GeV/SpinAlignment/Phi/Forest
  - /star/data01/pwg/sunxuhit/AuAu19GeV/SpinAlignment/Phi/Forest
  - /star/data01/pwg/sunxuhit/AuAu11GeV/SpinAlignment/Phi/Forest
- The corresponding run list can be found at:
  - /star/data01/pwg/sunxuhit/AuAu200GeV/SpinAlignment/Phi/List
  - /star/data01/pwg/sunxuhit/AuAu62GeV/SpinAlignment/Phi/List
  - /star/data01/pwg/sunxuhit/AuAu39GeV/SpinAlignment/Phi/List
  - /star/data01/pwg/sunxuhit/AuAu27GeV/SpinAlignment/Phi/List
  - /star/data01/pwg/sunxuhit/AuAu19GeV/SpinAlignment/Phi/List
  - /star/data01/pwg/sunxuhit/AuAu11GeV/SpinAlignment/Phi/List
- `FillSpinAlignment` is to produce histograms of phi-meson invariant mass
  - please change outputfile in `VecMesonSpinAlignment/FillSpinAlignment/StRoot/StVecMesonAna/StVecMesonAna.cxx` first
  - use cons to under `FillSpinAlignment` to compile the code
  - To run the code: root4star -b -q FillSpinAlignment.C\(energy,flag_ME,List,StartEvent,StopEvent,pid\)
    - mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
    - flag_ME: 0 for Same Event, 1 for Mixed Event
    - List: different number for different TTree list
    - StartEvent: first event number
    - StopEvent: last event number
    - pid: 0 for phi meson
  - The output file will have historamgs of invariant mass of K+ and K- for different pT, centrality and (cos(theta\*)) with different systematic cuts
- `CalSpinAlignment` is to calculate raw rho\_00 as a function of pT, centrality and (cos(theta\*)) with different systematic cuts
  - please change outputfile in each macro below
  - `subBackGround.C` subtract background from Mixed Event
    - root4star -b -q subBackGround.C++\(energy,pid\)
  - `calSpinAlignmentSys.C` extract raw spin alignment signal
    - root4star -b -q calSpinAlignmentSys.C++\(energy,pid\)
  - raw rho00 results can be found at:
    - /star/data01/pwg/sunxuhit/AuAu200GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
    - /star/data01/pwg/sunxuhit/AuAu62GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
    - /star/data01/pwg/sunxuhit/AuAu39GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
    - /star/data01/pwg/sunxuhit/AuAu27GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
    - /star/data01/pwg/sunxuhit/AuAu19GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
    - /star/data01/pwg/sunxuhit/AuAu11GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root
- all the correction (efficiency, event plane resolution and acceptance effect) can be found in Chensheng Zhou's code
