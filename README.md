# Emod3d_FWI
0. Make sure the main path folder has the following structure:
  AdjSims  FwdSims  Kernels  Model  PART1  PART2  StatInfo
  - The running folder: Kernels/Iters/iter1/
  - The parallel kernel calculation folder: PARTi/Kernels/Iters/iter1/
  - Observed data uploaded to: Kernels/Vel_ob_240s_HH_13s
  - SRF files for 13 events uploaded to: Kernels/Iters/iter1/SRF_13s
  - Source/ station coordinate files uploaded to: StatInfo
  - Model (Vs, Vp, rho) uploaded to: Model/Models/

1. Run Master_iterations_OPT_1_slurm_moduled.py to check for forward simulation and the parameter/ source files generated in advance in: Kernels/Iters/iter1/FWD/All_par and  Kernels/Iters/iter1/FWD/All_srf. Check log file in Rlog for more information of each event's forward simulation. Check waveform output at Kernels/Vel_opt/

2. Run Master_INV_pool_moduled.py for iteratively inversion. All itertion data saved in Kernels/Iters/iter1/Dump
