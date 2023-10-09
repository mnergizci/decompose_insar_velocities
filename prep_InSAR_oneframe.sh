#!/bin/bash

location=$1 # a1 or a285
batch_id=$2  # licsdir=[a1/a285]/FINALS/${batch_id}
frame=$3
geocdir=$4 
div_id=$5  # divdir=/nfs/a285/homes/eejdm/DIV/SouthIsland/${batch_id}${div_id}
seismic=$6

if [ ! -z ${seismic} ]; then
	div_id=${batch_id}${div_id}_${seismic}
else
        div_id=${batch_id}${div_id}
fi

if [ $location == 'a1' ]; then
	licsdir=/nfs/a1/eejdm/FINALS/${batch_id}
else
	licsdir=/nfs/a285/homes/eejdm/FINALS/${batch_id}
fi
divdir=/nfs/a1/eejdm/DIV/${div_id}

mkdir -p ${divdir}
if [ ! -d  ${divdir}/GNSS ]; then
  echo Adding GNSS to ${divdir}
  cp -rf /nfs/a1/eejdm/DIV/GNSS ${divdir}
  echo GNSS
fi

if [ ! -f ${divdir}/DIV.conf ]; then
  cp -rf /nfs/a1/eejdm/DIV/DIV.conf ${divdir}
fi

echo ${frame}
mkdir -p ${divdir}/${frame}

LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/E.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/E.geo.tif;
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/N.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/N.geo.tif;
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/U.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/U.geo.tif;

if [ -z ${seismic} ]; then

  tsdir=TS_${geocdir}
  LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/vel.geo.tif;
  cp ${licsdir}/${frame}/$tsdir/results/vel.png ${divdir}/${frame}/vel.geo.png
  LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vstd -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/vstd.geo.tif;
  cp ${licsdir}/${frame}/$tsdir/results/vstd.png ${divdir}/${frame}/vstd.geo.png
  LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/mask -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/mask.geo.tif;
  cp ${licsdir}/${frame}/$tsdir/results/mask.png ${divdir}/${frame}/mask.geo.png

else

  tsdir=TS_${geocdir}/results/seismic_vels
  if [ -d ${tsdir} ]; then
    echo ''
    echo 'Preparing data from LiCSBAS_seismic_inv.py'
    echo ''
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/pre_vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/pre_vel.png ${divdir}/${frame}/pre_vel.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/post_vel${seismic} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/post_vel${seismic}.png ${divdir}/${frame}/post_vel.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/pre_vel_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel_err.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/pre_vel_err.png ${divdir}/${frame}/pre_vel_err.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/post_vel${seismic}_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel${seismic}_err.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/post_vel${seismic}_err.png ${divdir}/${frame}/post_vel${seismic}_err.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/mask -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/mask.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/mask.png ${divdir}/${frame}/mask.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/coseismic${seismic} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/coseismic${seismic}.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/coseismic${seismic}.png ${divdir}/${frame}/coseismic${seismic}.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/coseismic${seismic}_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/coseismic${seismic}_err.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/coseismic${seismic}_err.png ${divdir}/${frame}/coseismic${seismic}_err.geo.png
  else
    tsdir=TS_${geocdir}
    echo ''
    echo 'No outputs from LiCSBAS_seismic_inv.py found. Simulating...'
    echo ''
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/vel.png ${divdir}/${frame}/pre_vel.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/vel.png ${divdir}/${frame}/post_vel.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vstd -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel_err.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/vstd.png ${divdir}/${frame}/pre_vel_err.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vstd -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel${seismic}_err.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/vstd.png ${divdir}/${frame}/post_vel${seismic}_err.geo.png
    LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/mask -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/mask.geo.tif;
    cp ${licsdir}/${frame}/$tsdir/results/mask.png ${divdir}/${frame}/mask.geo.png
  fi
  
cp ${licsdir}/${frame}/TS_${geocdir}/network/network13_nobad.png ${divdir}/${frame}/network.png;
  
fi


mask_param='coh_avg'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='n_unw'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='maxTlen'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='n_gap'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='stc'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='n_ifg_noloop'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='n_loop_err'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='resid_rms'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

mask_param='vstd'
LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/${mask_param} -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/${mask_param}.geo.tif;
cp ${licsdir}/${frame}/$tsdir/results/${mask_param}.png ${divdir}/${frame}/${mask_param}.geo.png

