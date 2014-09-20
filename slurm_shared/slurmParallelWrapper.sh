#!/bin/sh
# This wrapper script is intended to be submitted to SLURM to support parallel
# execution.
#
# This script uses the following environment variables set by the submit MATLAB code:
# MDCE_CMR            - the value of ClusterMatlabRoot (may be empty)
# MDCE_MATLAB_EXE     - the MATLAB executable to use
# MDCE_MATLAB_ARGS    - the MATLAB args to use
#
# The following environment variables are forwarded through mpiexec:
# MDCE_DECODE_FUNCTION     - the decode function to use
# MDCE_STORAGE_LOCATION    - used by decode function 
# MDCE_STORAGE_CONSTRUCTOR - used by decode function 
# MDCE_JOB_LOCATION        - used by decode function 

# Copyright 2010 MathWorks, Inc.

# Create full paths to mw_smpd/mw_mpiexec if needed
FULL_SMPD=${MDCE_CMR:+${MDCE_CMR}/bin/}mw_smpd
FULL_MPIEXEC=${MDCE_CMR:+${MDCE_CMR}/bin/}mw_mpiexec
SMPD_LAUNCHED_HOSTS=""
MPIEXEC_CODE=0

# Work out where we need to launch SMPDs given our hosts file - defines
# SMPD_HOSTS
chooseSmpdHosts() {
echo -e "\nchooseSmpdHosts()"

    # We need the SLURM_NODELIST value - the following line either echoes the value,
    # or aborts.
    echo Node file: ${SLURM_NODELIST:?"Node file undefined"}
    # We must launch SMPD on each unique host that this job is to run on. We need
    # this information as a single line of text, and so we pipe the output of "uniq"
    # through "tr" to convert newlines to spaces
    SMPD_HOSTS=`echo ${SLURM_NODELIST} | uniq | tr ',' ' '`
    SMPD_HOSTS=`scontrol show hostname ${SLURM_NODELIST} | uniq | tr '\n' ' '`
}

# Work out which port to use for SMPD
chooseSmpdPort() {
echo -e "\nchooseSmpdPort()"

    # Choose unique port for SMPD to run on. SLURM_JOBID is something like
    # 15.slurm-server-host.domain.com, so we extract the numeric part of that
    # using sed.
    JOB_NUM=`echo ${SLURM_JOBID:?"SLURM_JOBID undefined"} | sed 's#^\([0-9][0-9]*\).*$#\1#'`
    # Base smpd_port on the numeric part of the above
    SMPD_PORT=`expr $JOB_NUM % 10000 + 20000`
}

# Use ssh to launch the SMPD daemons on each processor
launchSmpds() {
echo -e "\nlaunchSmpds()"

    # Launch the SMPD processes on all hosts using SSH
    echo "Starting SMPD on ${SMPD_HOSTS} ..."


#    echo "srun --ntasks=${SLURM_NNODES} \"${FULL_SMPD}\" -s -phrase MATLAB -port ${SMPD_PORT}&"
#    srun --ntasks=${SLURM_NNODES} \"${FULL_SMPD}\" -s -phrase MATLAB -port ${SMPD_PORT} &


    for host in ${SMPD_HOSTS}
    do
        # This script assumes that SSH is set up to work without passwords between
        # all nodes on the cluster
	echo ssh $host \"${FULL_SMPD}\" -s -phrase MATLAB -port ${SMPD_PORT}
	ssh $host \"${FULL_SMPD}\" -s -phrase MATLAB -port ${SMPD_PORT}
	ssh_return=${?}
	if [ ${ssh_return} -ne 0 ]; then
            echo "Launching smpd failed for node: ${host}"
            exit 1
	else
            SMPD_LAUNCHED_HOSTS="${SMPD_LAUNCHED_HOSTS} ${host}"
	fi
    done
    echo "All SMPDs launched"
}

# Work out how many processes to launch - set MACHINE_ARG
chooseMachineArg() {
echo -e "\nchooseMachineArg()"

    MACHINE_ARG="-n ${SLURM_NTASKS}"
    MACHINE_ARG="-n ${SLURM_NPROCS}"
    echo "Machine args: $MACHINE_ARG"
}

runMpiexec() {
echo -e "\nrunMpiexec()"

    # As a debug stage: echo the command line...
    echo -e \"${FULL_MPIEXEC}\" -phrase MATLAB -port ${SMPD_PORT} \
        -l ${MACHINE_ARG} -genvlist \
        MDCE_DECODE_FUNCTION,MDCE_STORAGE_LOCATION,MDCE_STORAGE_CONSTRUCTOR,MDCE_JOB_LOCATION \
        \"${MDCE_MATLAB_EXE}\" ${MDCE_MATLAB_ARGS} "\n"
    # ...and then execute it
    eval \"${FULL_MPIEXEC}\" -phrase MATLAB -port ${SMPD_PORT} \
        -l ${MACHINE_ARG} -genvlist \
        MDCE_DECODE_FUNCTION,MDCE_STORAGE_LOCATION,MDCE_STORAGE_CONSTRUCTOR,MDCE_JOB_LOCATION \
        \"${MDCE_MATLAB_EXE}\" ${MDCE_MATLAB_ARGS}
    MPIEXEC_CODE=${?}
    echo {$MPIEXEC_CODE}
}

# Now that we have launched the SMPDs, we must install a trap to ensure that
# they are closed either in the case of normal exit, or job cancellation:
# Default value of the return code
cleanupAndExit() {
echo -e "\ncleanupAndExit()"

    echo ""
    echo "Stopping SMPD on ${SMPD_LAUNCHED_HOSTS} ..."


#    if [ -z ${SLURM_JOBID} ]; then
#	echo "Job did not properly start.  No need to bring down smpd."
#    else
#	scancel  ${SLURM_JOBID}.0
## srun  --ntasks_per_node=1 \"${FULL_SMPD}\" -shutdown -phrase MATLAB -port ${SMPD_PORT}
#    fi


    for host in ${SMPD_LAUNCHED_HOSTS}
    do
        echo ssh $host \"${FULL_SMPD}\" -shutdown -phrase MATLAB -port ${SMPD_PORT}
        ssh $host \"${FULL_SMPD}\" -shutdown -phrase MATLAB -port ${SMPD_PORT}
    done
    echo "Exiting with code: ${MPIEXEC_CODE}"
    exit ${MPIEXEC_CODE}
}

# Define the order in which we execute the stages defined above
MAIN() {
    trap "cleanupAndExit" 0 1 2 15
    chooseSmpdHosts
    chooseSmpdPort
    launchSmpds
    chooseMachineArg
    runMpiexec
    exit ${MPIEXEC_CODE}
}

# Call the MAIN loop
MAIN
