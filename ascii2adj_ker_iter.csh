#! /bin/csh

# Script to reformat a series of adjoint ascii files into a single binary file for emod3d-mpi

set BPATH = /home/rgraves/Bin

# specify output binary file and directory
set ADJOINT_FILE = my_adjoint_file.e3d
set OUTDIR = ../../../AdjSims/V3.0.7-a2a_xyz/Adj-SourceBin

\mkdir -p $OUTDIR

# required simulation parameters need to be the same as forward runs
set NT = 1500
set DT = 0.16
#set NT = 3000
#set DT = 0.08
set HH = 4.0
set MODEL_ROT = 0.0

# specification of input files and associated meta data follows

set INDIR = ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii
set STATCORDS = ../../../../StatInfo/STATION_dh_4km.txt
#set STATCORDS = /scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/INV_MarlVM_DH_4km_test_checker_LOWER/StatInfo/STATION_checker.txt
# note: station lon,lat info not used currently, but would be best to provide
#       correct locations for possible double-checking
set STATS = ( QRZ NNZ WEL DSZ THZ KHZ INZ LTZ GVZ OXZ )
set SLONS = ( 0 0 0 0 0 0 0 0 0 0 )
set SLATS = ( 0 0 0 0 0 0 0 0 0 0 )

# specification of components follows
#
# valid choices for adjoint components are: vx, vy and vz, which are then inserted at
# vx, vy, and vz nodes, respectively in emod3d-mpi
#
# COMPS are the components for the ascii input files and should match the .<comp> values given 
# in the filenames

# example for only 2 components
set ADJOINT_COMPS = ( vx vz )
set COMPS = ( x z )

# example for only 1 components
set ADJOINT_COMPS = ( vy )
set COMPS = ( y )

# example for all 3 components
set ADJOINT_COMPS = ( vx vy vz )
set COMPS = ( x y z )

# Done with input specifications, processing stuff follows

set ADJ_COMPS_IN = $ADJOINT_COMPS[1]
foreach comp ( $ADJOINT_COMPS[2-] )

set ADJ_COMPS_IN = ${ADJ_COMPS_IN},$comp

end

# temporary file for listing inputs, will be created by script
set FILELIST = ../../../AdjSims/V3.0.7-a2a_xyz/a2a_filelist.txt

\rm $FILELIST
set s = 0
foreach stat ( $STATS )
@ s ++

gawk -v ss=$stat -v nt=$NT -v slon=$SLONS[$s] -v slat=$SLATS[$s] '{ \
if(NF>=4 && ss==$4){ix=$1;iy=$2;iz=$3;}} \
END{printf "%5d %5d %5d %6d %.5f %.5f %s",ix,iy,iz,nt,slon,slat,ss;}' $STATCORDS >> $FILELIST

set a = 0
foreach comp ( $ADJOINT_COMPS )
@ a ++

echo -n " ${INDIR}/${stat}.${COMPS[$a]}" >> $FILELIST

end

echo "" >> $FILELIST

end

${BPATH}/ascii2adj filelist=$FILELIST outfile=$OUTDIR/$ADJOINT_FILE \
                   adjoint_comps=${ADJ_COMPS_IN} dt=$DT h=$HH modelrot=$MODEL_ROT >& ../../../AdjSims/V3.0.7-a2a_xyz/a2a.out
