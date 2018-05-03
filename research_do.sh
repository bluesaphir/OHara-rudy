#!/bin/sh

INITIAL="-60"
FINAL="60"
STRIDE="10"

for ITER in `seq $INITIAL $STRIDE $FINAL`
do
	mkdir -p "d""$ITER"
	cd "d""$ITER"
	cp ../state_2000_IKR_from_ORd.dat ./
	mkdir -p tests_ctrl
	../integration "$ITER"
	cd ..
	echo "Calculation with ""$ITER"" is done *******"
done

