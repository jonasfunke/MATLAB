function slurmDestroyJob(scheduler, job)
%slurmDestroyJob Destroys job files on a remote host.
%
% Set your schedulers's DestroyJobFcn to this function using the following
% command (see README):
%     set(sched, 'DestroyJobFcn', @slurmDestroyJob);

%  Copyright 2010 MathWorks, Inc.

jobids = scheduler.getJobSchedulerData(job);
for jidx = 1:length(jobids)
    cmd = sprintf('scancel %d', jobids(jidx));
    s = system(cmd);
    if s~=0
        warning(['Failed to cancel job: ' num2str(jobids(jidx))]) %#ok<WNTAG>
    end
end
