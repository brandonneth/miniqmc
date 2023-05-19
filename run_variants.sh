#!/bin/bash

#Runs the variants of the miniqmc application as indicated by the command line arguments.
# the -v flag indicates a variant to run, with possible values 0 through 6.
# the -r flag indicates that the raja only variant should be run as well.
# the --append flag indicates that the result should be appended to the output file.
# the -n flag indicates the number of executions to time
# the -g flag indicates the shape to use for the run. this should be in quotations.
# the -d flag indicates the delay size. Default is 32.
VARIANTS=""
APPEND=0
NUMRUNS=1
DELAY=32
SHAPE="1 2 2"
OUTFILE=results.csv
while test $# -gt 0
do
        case "$1" in
                --append) echo "Appending to results file";
                        APPEND=1
                        ;;
                -v) echo "Adding variant $2";
                        VARIANTS="$VARIANTS $2";
                        shift
                        ;;
                -r) echo "Running";
                        VARIANTS="$VARIANTS ";
                        ;;
		-n) echo "Running $2 executions";
			NUMRUNS=$2
			shift
			;;
		-g) echo "Running with group shape $2";
			SHAPE=$2
			shift
			;;
		-d) echo "Using delay count of $2";
			DELAY=$2
			shift
			;;
		-o) echo "output file $2";
			OUTFILE=$2
			shift
			;;
                *) echo "unknown argument: $1"; exit;
        esac
        shift
done


if [[ $APPEND -ne 1 ]] ; then
rm $OUTFILE
touch $OUTFILE
echo "Shape, Delay, Variant, Time" > $OUTFILE
fi


for variant in $VARIANTS; do
echo "Running variant $variant";
pushd ~/miniqmc/build
	rm -rf ./*
	cmake -DENABLE_OPENMP=On -DBLT_CXX_STD=c++17 -DCMAKE_INSTALL_PREFIX=$HOME/libs/quartz .. -DVARIANT=$variant
	make -j10

	for run in $(seq $NUMRUNS); do
		echo ./bin/miniqmc -g "$SHAPE" -k $DELAY
		./bin/miniqmc -g "$SHAPE" -k $DELAY 
		xmlname="info"
		for dim in $SHAPE; do
			xmlname=${xmlname}_$dim
		done
		xmlname=${xmlname}.xml
		
		time=$(python3 ../extract_update.py $xmlname)
		echo "$SHAPE, $DELAY, $variant, $time"
		echo "$SHAPE, $DELAY, $variant, $time" >> ../$OUTFILE
	done
popd
done





