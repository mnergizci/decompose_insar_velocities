#!/bin/bash

batch_id=default_lpc3rad
seismic='y'

geocdir=$1

if [ ${seismic} == 'y' ]; then
	div_id=${batch_id}_seismic
	frames=TOPS.txt
else
	frames=DIV_frames.txt
fi

licsdir=/nfs/a285/homes/eejdm/FINALS/${batch_id}
divdir=/nfs/a285/homes/eejdm/DIV/SouthIsland/${div_id}

mkdir -p ${divdir}
cp -rf /nfs/a285/homes/eejdm/DIV/SouthIsland/GNSS ${divdir}

if [ ! -f ${divdir}/DIV.conf ]; then
cp -rf /nfs/a285/homes/eejdm/DIV/SouthIsland/DIV.conf ${divdir}
fi

echo GNSS added to ${divdir}


for frame in $(cat ${frames}); do
	echo $frame
	mkdir -p ${divdir}/${frame}
	
        LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/E.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/E.geo.tif;
        LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/N.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/N.geo.tif;
        LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/${geocdir}/U.geo -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/U.geo.tif;
        cp ${licsdir}/${frame}/$tsdir/network/network13_nobad.png ${divdir}/${frame}/network.png;

	if [ ! ${seismic} == 'y' ]; then

	        tsdir=TS_${geocdir}
		LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/vel.geo.tif;
		cp ${licsdir}/${frame}/$tsdir/results/vel.png ${divdir}/${frame}/vel.geo.png
		LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/vstd -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/vstd.geo.tif;
        	cp ${licsdir}/${frame}/$tsdir/results/vstd.png ${divdir}/${frame}/vstd.geo.png
		LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/mask -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/mask.geo.tif;
        	cp ${licsdir}/${frame}/$tsdir/results/mask.png ${divdir}/${frame}/mask.geo.png

	else

		tsdir=TS_${geocdir}/results/seismic_vels	
		LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/pre_vel -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/pre_vel.png ${divdir}/${frame}/pre_vel.geo.png
                LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/post_vel20161113 -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/post_vel20161113.png ${divdir}/${frame}/post_vel.geo.png
		LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/pre_vel_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/pre_vel_err.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/pre_vel_err.png ${divdir}/${frame}/pre_vel_err.geo.png
                LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/post_vel20161113_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/post_vel20161113_err.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/post_vel20161113_err.png ${divdir}/${frame}/post_vel20161113_err.geo.png
                LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/results/mask -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/mask.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/results/mask.png ${divdir}/${frame}/mask.geo.png
                LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/coseismic20161113 -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/coseismic20161113.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/coseismic20161113.png ${divdir}/${frame}/coseismic20161113.geo.png
                LiCSBAS_flt2geotiff.py -i ${licsdir}/${frame}/$tsdir/coseismic20161113_err -p ${licsdir}/${frame}/${geocdir}/EQA.dem_par  -o ${divdir}/${frame}/coseismic20161113_err.geo.tif;
                cp ${licsdir}/${frame}/$tsdir/coseismic20161113_err.png ${divdir}/${frame}/coseismic20161113_err.geo.png
	fi
done
