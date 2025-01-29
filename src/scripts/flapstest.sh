# flaps test script
# This script runs a sample of simple flaps jobs to test the functionality
# of flaps.  The script compares the results to expected results and outputs
# summary results for the test.
#
export TEST_STATUS=0  # 0 - Test ran successfully, 1 - Test not run, 2 - Test Failed
export FILE=0
export CMD_RESULT=0
export TEST_BASE=0
export TEST_COUNT=0
export FAILED_COUNT=0
export SUCCESS_COUNT=0
export NOTTESTED_COUNT=0
export FLAPS_CMD=0
export MAIL_CMD=0
export DEBUG_PRINT=1

debug () {
	if [ ${DEBUG_PRINT}z != z ] ; then
		echo $@
	fi
}

usage () {
   echo "Usage: `basename $0` [-a][-c config][-d][-r froot]" >&2
   echo " options:" >&2
   echo "    -d directory Directory where the tests are found" >&2 
   exit 1
}

testit () {
	debug ">Entering testit for file: " $FILE  
	TEST_BASE=`basename $FILE .ax`
	if [ ${FILE##*.} = ax ] ; then
		# the file ended with .ax
		TEST_COUNT=$(($TEST_COUNT+1))
		debug ">  Will attempt to run test on: " $FILE
		echo "  " >> $WORK_DIR/flapstest.sum
		echo "Testing job: ${FILE}" >> $WORK_DIR/flapstest.sum
      jobtest
		if [ $TEST_STATUS -eq 0 ] ; then
			echo "   Test for job: ${FILE}  was successful" >> $WORK_DIR/flapstest.sum
			debug ">>>  Test for job: ${FILE}  was successful"
			SUCCESS_COUNT=$(($SUCCESS_COUNT+1))
		elif [ $TEST_STATUS -eq 1 ] ; then
			echo "   Test for job: ${FILE} was not run: no validation files were present" >> $WORK_DIR/flapstest.sum 
			debug ">>>  Test for job: ${FILE} was not run: no validation files were present" 
			NOTTESTED_COUNT=$(($NOTTESTED_COUNT+1))
		else 
			echo "   Tests for job: ${FILE} failed" >> $WORK_DIR/flapstest.sum
			debug ">>>  Tests for job: ${FILE} failed"
			FAILED_COUNT=$(($FAILED_COUNT+1))
		fi
	fi

}
jobtest () {
	debug ">  Entering jobtest for: " $FILE
#  initially set TEST_STATUS to 1 (test not run)
	TEST_STATUS=1

	debug ">    in directory: $PWD"
	
	debug ">    Checking for presence of validation files"

	for VALID_FILE in `ls *.valid` ; do
		debug ">  ${VALID_FILE}"
		COMPARE_FILE=`basename ${VALID_FILE} .valid`
		debug ">  $COMPARE_FILE"
		debug ">    removing $COMPARE_FILE"
	   rm $COMPARE_FILE
# validation files found we can run the test
		TEST_STATUS=0
	done

	if [ ${TEST_STATUS} -ne 1 ];then
	   debug ">    Validation files were present"
		debug ">    Running test: ${FLAPS_CMD} ${TEST_BASE}.ax > ${TEST_BASE}.out 2> ${TEST_BASE}.err"
		${FLAPS_CMD} ${TEST_BASE}.ax > ${TEST_BASE}.out 2> ${TEST_BASE}.err
		CMD_RESULT=$?
		if [ ${CMD_RESULT} -ne 0 ] ; then
		   TEST_STATUS=2
		fi
		debug ">    Returned: $CMD_RESULT"

		for VALID_FILE in *.valid; do
	   	COMPARE_FILE=`basename ${VALID_FILE} .valid`
			debug ">    Validating files with diff: ${COMPARE_FILE} ${VALID_FILE}"
			FILE_EXTENSION=${VALID_FILE%%.valid}
			FILE_EXTENSION=${FILE_EXTENSION##*.}
			if [ -f ${COMPARE_FILE} ];then
				if [ ${FILE_EXTENSION} = "esa" ]; then
					../bin/compare_esa.py 0.002 ${COMPARE_FILE} ${VALID_FILE} 2>> $WORK_DIR/flapstest.sum
					CMD_RESULT=$(($?))
				else
				   diff ${COMPARE_FILE} ${VALID_FILE} 
					CMD_RESULT=$(($?))
				fi
				echo "cmd_result: ${CMD_RESULT}"
			else
				echo "   ${COMPARE_FILE} is missing but required to validate ${VALID_FILE}" >> $WORK_DIR/flapstest.sum
				TEST_STATUS=2
			fi
	   	if [ ${CMD_RESULT} -ne 0 ] ; then
				echo "inside test for nonzero cmd_result"
	      	TEST_STATUS=2
				echo "   diff of ${COMPARE_FILE} and ${VALID_FILE} showed a difference" >> $WORK_DIR/flapstest.sum
				echo "   ${CMD_RESULT}" >> $WORK_DIR/flapstest.sum
			fi
	   	debug ">    Returned: $CMD_RESULT"
		done
	fi
}										
# iterate through all of the files in the path provided as an argument
# if no path is provided the current directory is used.

while getopts ":d:r:" opt
do
   case $opt in
      d) DIR_NAME=$OPTARG 
			export DIR_NAME
         ;;
		r) FROOT=$OPTARG
		   export FROOT
			;;
      ?) usage;;
   esac
done

export WORK_DIR=$PWD 
if [[ -n $DIR_NAME ]] ; then
   if [ `echo $DIR_NAME | cut -c 1` = "/" ] ; then
		export BASE_DIR=${DIR_NAME}"/"
	else
		export BASE_DIR=${PWD}/${DIR_NAME}"/"
	fi
else
	export BASE_DIR=${PWD}"/" 
fi
debug ">Working directory " $WORK_DIR
debug ">Base directory for testing " $BASE_DIR
cd $BASE_DIR
#set -x
echo "flaps test shell output"  > $WORK_DIR/flapstest.sum
date >> $WORK_DIR/flapstest.sum
echo "------------------------------------------" >> $WORK_DIR/flapstest.sum

# Set platform specific environment variables that specify
# the flaps and mail commands
#---------------------------------------------------------
os="`uname -s`"
if [ -f /unicos ];then
   os=cray;
fi
case $os in
   AIX)
      MAIL_CMD=/usr/bin/mail
		FLAPS_CMD="/site/bin/flaps"
      ;;
   cray)
      MAIL_CMD=/bin/mail
		FLAPS_CMD="/u/share/flaps"
      ;;
   IRIX64)
      MAIL_CMD=/usr/sbin/Mail
		FLAPS_CMD="/u/share/flaps"
      ;;
   Linux)
      MAIL_CMD=/usr/bin/Mail
		FLAPS_CMD="/boeing/bin/flaps -v alpha"
      ;;
esac

for DIR in `ls`  
do
	debug ">dirname" $DIR
	if [[ -d $DIR ]] ; then

		# the file is a directory
		cd $DIR

		for FILE in `ls` 
		do
			if [ ${FILE##*.} = ax ] ; then
				testit 
			fi
      done
      cd ..
	elif [ ${DIR##*.} = ax ] ; then
		FILE=$DIR 
		testit 
	fi

done
echo " " >> $WORK_DIR/flapstest.sum
echo "Successful tests: ${SUCCESS_COUNT} " >> $WORK_DIR/flapstest.sum
echo "Skipped tests:    ${NOTTESTED_COUNT} " >> $WORK_DIR/flapstest.sum
echo "Failed tests:     ${FAILED_COUNT} " >> $WORK_DIR/flapstest.sum
echo "----------------------- " >> $WORK_DIR/flapstest.sum
echo "Total tests:      ${TEST_COUNT} " >> $WORK_DIR/flapstest.sum
# Send the error messages via e-mail
if [ $TEST_STATUS != 0 ] ; then
	   ssh tallis "$MAIL_CMD $LOGNAME < $WORK_DIR/flapstest.sum"
fi
 
