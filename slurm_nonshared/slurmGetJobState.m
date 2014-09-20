function state = slurmGetJobState(scheduler, job, state)
%slurmGetJobState Gets the state of a job on a cluster.
%
% Set your schedulers's GetJobStateFcn to this function using the following
% command (see README):
%     set(sched, 'GetJobStateFcn', @slurmGetJobState)

%  Copyright 2010 MathWorks, Inc.
disp(['slurmGetJobState: ' state])
%scheduler
%job


mlock;
persistent jobsToMonitorNames;

if isempty(jobsToMonitorNames)
    jobsToMonitorNames = {};
end

if isempty(scheduler.UserData)
    if ~iscell(scheduler.CommunicatingSubmitFcn) || length(scheduler.CommunicatingSubmitFcn) < 3
        error('distcomp:genericscheduler:SubmitFcnError',...
            'SubmitFcn must include clusterHost and remoteDataLocation as extra arguments.');
    end
    scheduler.UserData = { scheduler.CommunicatingSubmitFcn{2} ; scheduler.CommunicatingSubmitFcn{3} };
end


jobName = job.Name; %job.pGetEntityLocation;

if strcmp(state, 'finished')
    jobsToMonitorNames = setxor(jobsToMonitorNames, {jobName});
    disp('Job finished... deleted from jobsMonitorNames.')
else
    if strcmp(state, 'queued') || strcmp(state, 'running')
        if ~ismember(jobName, jobsToMonitorNames)
            disp('Job is not in Monitor')
            jobsToMonitorNames = union(jobsToMonitorNames, {jobName});
            jobStateTimer = timer('Period', 10.0, 'ExecutionMode', 'fixedRate');
            set(jobStateTimer, 'TimerFcn', { @copyJobFilesIfFinished, scheduler, job });
            start(jobStateTimer);
        end
    end
end
