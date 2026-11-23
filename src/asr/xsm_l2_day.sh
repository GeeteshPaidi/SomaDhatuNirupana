#!/bin/bash
# usage: run_day.sh /mnt/c/.../ch2_xsm_YYYYMMDD_v1.zip
set -e
export xsmdas=$HOME/xsm/xsmdas
export PATH="$xsmdas/bin:$xsmdas/scripts:$PATH"
export LD_LIBRARY_PATH="$xsmdas/lib/:${LD_LIBRARY_PATH:-}"
z=$1; b=$(basename $z .zip)
OUT=/mnt/c/Users/geete/ch2_data/xsm_l2
[ -s $OUT/${b}_level2.pha ] && { echo "skip $b"; exit 0; }
W=$(mktemp -d); trap "rm -rf $W" EXIT; cd $W
cp $xsmdas/pfiles/*.par $W/; export PFILES=$W
python3 -c "import zipfile,sys;zipfile.ZipFile(sys.argv[1]).extractall('.')" $z
R=$(dirname $(find . -name '*_level1.fits'))
L1=$R/${b}_level1.fits; HK=$R/${b}_level1.hk; SA=$R/${b}_level1.sa
xsmgtigen hkfile=$HK safile=$SA gtifile=$b.gti filterfile=$xsmdas/pfiles/std_filterfile.dat usergtifile=- usesunang=yes clobber=yes > log 2>&1
xsmgenspec l1file=$L1 specfile=${b}_level2.pha spectype=time-resolved tstart=0 tstop=0 tbinsize=60 hkfile=$HK safile=$SA gtifile=$b.gti genarf=no areascal=yes chantype=PI clobber=yes >> log 2>&1
cp ${b}_level2.pha $OUT/
tail -3 $W/log; cd /; rm -rf $W; echo "done $b"
