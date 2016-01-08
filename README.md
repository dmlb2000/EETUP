# EMSL Experiment/Theory Unification Project

Originally cloned from [EMSLHub](https://emslhub.emsl.pnl.gov/tools/nmrtools),
this project is trying to support workflows that bridge scientific theory with
experimental results.

## Current Workflows

### NMR and NWChem Workflow

This workflow utilizes an NWChem output to generate an NMR signal that can be
compared to NMR singals. NMRPipe is the NMR analysis tool used on the
experimental results. NMRPipe supports various instruments see their [home
page](http://spin.niddk.nih.gov/NMRPipe/) for more details.

The software needed to make this workflow run is found in the following
repositories:

 1. [SWADL](https://github.com/dmlb2000/swadl-library)
 1. [Kepler SWADL Components](https://github.com/dmlb2000/kepler-swadl)
 1. [NWChem CML](https://github.com/dmlb2000/nwchem-cml)


