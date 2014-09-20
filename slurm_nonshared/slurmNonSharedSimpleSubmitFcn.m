function slurmNonSharedSimpleSubmitFcn(scheduler, job, props, ...
    clusterHost, remoteDataLocation)
%SLURMNONSHAREDSIMPLESUBMITFCN Submits a MATLAB job to a SLURM scheduler in the
% absence of a shared file system between the MATLAB client and the
% SLURM scheduler.
%
% See also slurmNonSharedSimpleDecodeFcn.

% Copyright 2010 MathWorks, Inc.
disp('slurmNonSharedSimpleSubmitFcn')

if ~ischar(clusterHost)
    error('distcomp:genericscheduler:SubmitFcnError', ...
        'Hostname must be a string');
end
if ~ischar(remoteDataLocation)
    error('distcomp:genericscheduler:SubmitFcnError', ...
        'Remote Data Location must be a string');
end
scheduler.UserData = { clusterHost ; remoteDataLocation };
localDataLocation = scheduler.JobStorageLocation; %scheduler.DataLocation;

 

% Set the name of the decode function which will be executed by
% the worker. The decode function must be on the path of the MATLAB
% worker when it starts up. This is typically done by placing the decode
% function in MATLABROOT/toolbox/local on the cluster nodes, or by
% prefixing commandToRun (created below) with a command to cd to the
% directory where the decode function's MATLAB file exists.
decodeFcn = 'slurmSimpleDecode';

% Read the number of tasks which are to be created. This property
% cannot be changed.
numberOfTasks = props.NumberOfTasks;

% A unique file and directory name for the job. This is used to create
% files and a directory under scheduler.DataLocation
jobLocation = props.JobLocation;

% A cell array of unique file names for tasks. These are used to create
% files under jobLocation
taskLocations = props.TaskLocations;

% Since SLURM jobs will be submitted from a UNIX host on the cluster,
% a single quoted string will protect the MATLAB command.
quote = '''';

% The MATLAB command to be run on a cluster node to execute a task.
commandToRun = [quote props.MatlabExecutable quote ' ' ...
    props.MatlabArguments];


% Copy the matlab_metadata.mat file to the remote host.
localMetaDataFile = [localDataLocation '/matlab_metadata.mat'];
copyDataToCluster(localMetaDataFile, remoteDataLocation, clusterHost);

% Copy the local job directory to the remote host.
localJobDirectory = [localDataLocation '/' jobLocation];
copyDataToCluster(localJobDirectory, remoteDataLocation, clusterHost);

% Copy the local job files to the remote host.
localJobFiles = [localDataLocation '/' jobLocation '.*'];
copyDataToCluster(localJobFiles, remoteDataLocation, clusterHost);



% Submit tasks which the scheduler will execute by starting MATLAB workers.
for i = 1:numberOfTasks
    taskLocation = taskLocations{i};
    remoteJobDir = [remoteDataLocation '/' jobLocation];
    
    % Create a script to submit a SLURM job.
    localScript = createSLURMSubmitScript(commandToRun, decodeFcn, ...
        props.StorageConstructor, remoteDataLocation, ...
        jobLocation, taskLocation, quote, job.ID, job.Tasks(i).ID);
    [~, scriptName] = fileparts(localScript);
    
    % Copy the script to the remote host.
    copyDataToCluster(localScript, remoteJobDir, clusterHost);
    
    % Create the command to run on the remote host.
    remoteScriptLocation = [ remoteJobDir '/' scriptName ];
    remoteCommand = sprintf('sbatch --ntasks=%d %s', props.NumberOfTasks, ...
        remoteScriptLocation);
    
    % Execute the submit command on the remote host.
    runCmdOnCluster(remoteCommand, clusterHost)
    
    % Delete the local copy of the script
    delete(localScript)
end


function filename = createSLURMSubmitScript(commandToRun, decodeFunction, ...
    storageConstructor, remoteDataLocation, jobLocation, ...
    taskLocation, quote, jobID, TaskID)
%Create a SLURM submit script that forwards the required environment
%variables and runs a MATLAB worker.
logFileLocation = [quote remoteDataLocation '/' taskLocation '.log' quote];
errFileLocation = [quote remoteDataLocation '/' taskLocation '.err' quote];

setEnv = sprintf([ ...
    '#!/bin/sh', '\n', ...
    '#SBATCH --job-name=Job', num2str(jobID), '.', num2str(TaskID), '\n', ...
    '#SBATCH --output=', logFileLocation, '\n', ...
    '#SBATCH --error=', errFileLocation, '\n', ...
    '#SBATCH --cpus-per-task=1\n', ... %added
    '#SBATCH --nodelist=compute-0-15,compute-0-16\n', ... % added
    ['export MDCE_DECODE_FUNCTION=' decodeFunction, '\n' ...
    'export MDCE_STORAGE_LOCATION=' remoteDataLocation, '\n' ...
    'export MDCE_STORAGE_CONSTRUCTOR=' storageConstructor, '\n' ...
    'export MDCE_JOB_LOCATION=' jobLocation, '\n' ...
    'export MDCE_TASK_LOCATION=' taskLocation, '\n' ...
    'export MDCE_DEBUG=true\n']
    ]);

% Content of script.
scriptContent = sprintf('%s%s \n', setEnv, commandToRun);

% Create script.
filename = tempname;
% Open file in binary mode to make it cross-platform.
fid = fopen(filename, 'w');
if fid<0
    error(['Failed to open file: ' filename])
end
fprintf(fid, scriptContent);
fclose(fid);
