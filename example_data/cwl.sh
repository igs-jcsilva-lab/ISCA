#!/bin/sh

export WDIR=/tmp/example
export TMPDIR=$WDIR/tmp
export OUTDIR=$WDIR/out
export TMPOUTDIR=$WDIR/tmp/out

mkdir $TMPDIR $OUTDIR $TMPOUTDIR

/usr/local/packages/python-3.7.7/bin/cwltool \
  --tmpdir-prefix $TMPDIR \
  --tmp-outdir-prefix $TMPOUTDIR \
  --outdir $OUTDIR \
  --disable-color \
  /usr/local/scratch/TA-test/TA1/targeted_assembly.cwl \
  $WDIR/targeted_assembly.yml

