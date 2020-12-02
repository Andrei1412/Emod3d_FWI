#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH --error=test.e   # stderr file
#SBATCH --output=test.o   # stdout file
#SBATCH --exclusive

###***#SBATCH --job-name=nr02   # job name
###***#SBATCH --partition=   # queue to run in

echo "Starting" $SLURM_JOB_ID `date`
echo "Initiated on `hostname`"
echo ""
cd "$SLURM_SUBMIT_DIR"           # connect to working directory of sbatch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

set VERSION = 3.0.7-dev
set EMOD_DIR = /home/rgraves/Mpi/Emod3d/V3.0.7-dev

set OUTBIN = ../../../FwdSims/V3.0.7-dc_stf_file/OutBin
set FWD = ../../../FwdSims/V3.0.7-dc_stf_file

set MY_RUN_ID = `echo $user $SLURM_JOB_ID | gawk '{split($2,a,".");printf "%s-%s\n",$1,a[1];}'`
echo "MY_RUN_ID=" $MY_RUN_ID

#Run fwd emod3d for all sources for pyflex pick-up/ misfit calculation and optimization
set SOURCES = ( 1   2   3   4   5   6   7   8   9   10  11  12  13)

set f = 0
foreach source ( $SOURCES )
@ f ++
#e3d.par file for each source has been prepared in advanced in Master_iterations_OPT_1_slurm_Moduled.py
#python -c "import os;
#     os.system(cp FWD/All_par/e3d_mysource_'+str($source)+'_opt.par FWD/e3d_mysource_i_opt.par);
#     os.system(cp FWD/All_srf/srf_source_'+str($source)+'_opt.srf srf_file.srf);"
echo "source=" ${source}
cp FWD/All_par/e3d_mysource_${source}_opt.par FWD/e3d_mysource_i_opt.par
cp FWD/All_srf/srf_source_${source}.srf srf_file.srf
rm ${OUTBIN}/*.*

srun $EMOD_DIR/emod3d-mpi par=FWD/e3d_mysource_i_opt.par < /dev/null
echo "emod3d 1 finished"

#Store waveform in ../../Vel_opt/Vel_ob_i for i-source
set Vel_ob_i = ../../../Kernels/Vel_opt/Vel_ob_${source}
rm ${Vel_ob_i}/*.*
python -c "from qcore import timeseries as ts; ts.LFSeis('$OUTBIN').all2txt(prefix='$Vel_ob_i/')"
echo "winbin aio finished"

end

echo "Done" `date`
exit
