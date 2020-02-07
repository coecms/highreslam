Reconfiguration Modifications
=============================

Two modifications are made to the normal reconfiguration

1. Rather than reading the land mask on one processor, scattering, then
   gathering again back to rank 0 in rcf_setup_lsm_out.f90, the land mask is
   just read on rank 0. This avoids a strange alignment error in the gather.

2. The 'spiral_circle' extrapolation method used to setup coastal grid points
   that don't have a corresponding point on the source grid is performed
   offline. This allows us to use the high-memory low cpu count 'hugemem' queue
   for the majority of the procssing, and the normal-memory high cpu count
   'normal' queue for the coastal extrapolation which is easily parallelized.

# Setup

In the Rose suite, run task 'fcm_make' to download the source code onto Gadi

# Run `0_patch_write_spiral_circle.sh $RUNID`

This will apply a patch to the reconfiguration source code of your run

# Run Rose tasks 'fcm_make2' and 'Aus_d0036_RA2M_um_recon'

This will write the information from the spiral_circle function to the job run
dir (e.g. `~/cylc-run/u-bq574/work/20170326T1800Z/Aus_d0036_RA2M_um_recon/`)

# Run `1_offline_spiral_circle.sh $RUNID`

This will build the spiral_circle offline processor, and submit a queue job to
process the outputs of write_spiral_circle

# Run `2_patch_read_spiral_circle.sh $RUNID`

This will apply a patch to the reconfiguration source code of your run

# Run Rose tasks 'fcm_make2' and 'Aus_d0036_RA2M_um_recon'

This will read the data that has been processed offline and continue with the
rest of the reconfiguration
