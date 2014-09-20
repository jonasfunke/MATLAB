function slurmParallelSubmitFcn( scheduler, job, props, varargin )
%SLURMPARALLELSUBMITFCN Submit parallel job to SLURM.

% Copyright 2010 MathWorks, Inc.

setenv('MDCE_DECODE_FUNCTION', 'slurmParallelDecode')
setenv('MDCE_STORAGE_LOCATION', props.StorageLocation)
setenv('MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor)
setenv('MDCE_JOB_LOCATION', props.JobLocation)
% Ask the workers to print debug messages by default
setenv('MDCE_DEBUG', 'true')

% Set this so that the script knows where to find MATLAB, MW_SMPD and MW_MPIEXEC
% on the cluster. This might be empty - the wrapper script will deal with that
% eventuality.
setenv('MDCE_CMR', scheduler.ClusterMatlabRoot)

% Set this so that the script knows where to find MATLAB, SMPD and MPIEXEC on
% the cluster. This might be empty - the wrapper script must deal with that.
setenv('MDCE_MATLAB_EXE', props.MatlabExecutable)
setenv('MDCE_MATLAB_ARGS', props.MatlabArguments)

% The wrapper script is in the same directory as this MATLAB file
dirpart = fileparts(mfilename('fullpath'));
scriptName = fullfile(dirpart, 'slurmParallelWrapper.sh');

% Forward the total number of tasks we're expecting to launch
setenv('MDCE_TOTAL_TASKS', num2str( props.NumberOfTasks ))

% Choose a file for the output.
logFile = fullfile( scheduler.DataLocation, sprintf('Job%d', job.ID), ...
    sprintf( 'Job%d.out', job.ID ) );
errFile = fullfile( scheduler.DataLocation, sprintf('Job%d', job.ID), ...
    sprintf( 'Job%d.err', job.ID ) );
cmdLine = sprintf( ...
    'sbatch --ntasks=%d --job-name Job%d --output="%s" --error="%s" "%s"', ...
    props.NumberOfTasks, job.ID, logFile, errFile, scriptName );

[s, w] = system(cmdLine);
if s==0
    % The output of successful submissions shows the SLURM job identifier
    fprintf(1, 'Job output will be written to: %s\nSRUN output: %s\n', ...
        logFile, w);
    w = strtrim(w);
    while true
        % In order to be able to cancel jobs in the furture, we
        % need to keep track of the job ids.
        [jobid, w] = strtok(w); %#ok<STTOK>
        if isempty(w), break, end
    end
    scheduler.setJobSchedulerData(job, jobid)
else
    % Report an error if the script did not execute correctly.
    warning( 'distcompexamples:generic:SLURM', ...
        'Submit failed with the following message:\n%s', w)
end
