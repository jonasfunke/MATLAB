function slurmSimpleSubmitFcn(scheduler, job, props, varargin)
%SLURMSIMPLESUBMITFCN Submit job to SLURM.

% Copyright 2010 MathWorks, Inc.

setenv('MDCE_DECODE_FUNCTION', 'slurmSimpleDecode')
setenv('MDCE_STORAGE_LOCATION', props.StorageLocation)
setenv('MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor)
setenv('MDCE_JOB_LOCATION', props.JobLocation)
% Ask the workers to print debug messages by default
setenv('MDCE_DEBUG', 'true')

% Tell the script what it needs to run. These two properties will
% incorporate ClusterMatlabRoot if it is set.
setenv('MDCE_MATLAB_EXE', props.MatlabExecutable)
setenv('MDCE_MATLAB_ARGS', props.MatlabArguments)

% The wrapper script is in the same directory as this MATLAB file
dirpart = fileparts(mfilename('fullpath'));
scriptName = fullfile(dirpart, 'slurmSimpleWrapper.sh');

% Submit the wrapper script to SLURM once for each task, supplying a different
% environment each time.
ntasks = props.NumberOfTasks;
jobids = nan(ntasks,1);
for i = 1:ntasks
    fprintf('Submitting task %i\n', i);
    setenv('MDCE_TASK_LOCATION', props.TaskLocations{i});
    % Choose a file for the output.
    logFile = fullfile( scheduler.DataLocation, sprintf('Job%d', job.ID), ...
        sprintf( 'Job%d_Task%d.out', job.ID, job.Tasks(i).ID ) );
    errFile = fullfile( scheduler.DataLocation, sprintf('Job%d', job.ID), ...
        sprintf( 'Job%d_Task%d.err', job.ID, job.Tasks(i).ID ) );
    cmdLine = sprintf( 'sbatch --job-name=Job%d.%d --output="%s" --error="%s" "%s"', ...
        job.ID, job.Tasks(i).ID, logFile, errFile, scriptName );
    
    [s, w] = system(cmdLine);
    if s==0
        % The output of successful submissions shows the SLURM job identifier
        fprintf(1, 'Job output will be written to: %s\nSRUN output: %s\n', ...
            logFile, w);
        w = strtrim(w);
        while true
            % In order to be able to cancel jobs in the furture, we
            % need to keep track of the job ids.
            [jobid, w] = strtok(w); %#ok
            if isempty(w), break, end
        end
        jobids(i,1) = str2double(jobid);
    else
        warning( 'distcompexamples:generic:SLURM', ...
            'Submit failed with the following message:\n%s', w)
    end
end

scheduler.setJobSchedulerData(job, jobids)
