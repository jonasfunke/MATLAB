%% 
clear all, close all, clc

mypool = parpool('InstallTest',1)


%%
clear all, close all, clc

cd('/Users/jonasfunke/Documents/MATLAB/parallel_test')
mycluster = parcluster('InstallTest');

%parpool(mycluster, 2);

n_worker=[ 1 ];
N = length(n_worker);
r = zeros(N,2);
%for n=1:N
n = 1;
    disp(['------------------------------- Running with ' num2str(n_worker(n)) ' workers -------------------------------'])
    %j = batch(mycluster, @test_function,1,'CaptureDiary',true, 'CurrentDirectory', '.','AdditionalPaths', {'/nfs/matlabuser'}, 'Pool', n_worker(n));
    j = batch(mycluster, @test_function,1,'CaptureDiary',true, 'CurrentDirectory', '.','AdditionalPaths', {'/nfs/matlabuser'}, 'Pool', 2);

    wait(j)
    r(n,1) = n_worker(n);
    r(n,2) = j.fetchOutputs{1};% Get results into a cell array
   % delete(j)
%end


%%
clear all, close all, clc

cd('/Users/jonasfunke/Documents/MATLAB/parallel_test')
mycluster = parcluster('InstallTest');

%parpool(mycluster, 2);

n_worker=[63 49 ];
N = length(n_worker);
r = zeros(N,2);
for n=1:N
    disp(['------------------------------- Running with ' num2str(n_worker(n)) ' workers -------------------------------'])
    j = batch(mycluster, @test_function,1,'CaptureDiary',true, 'CurrentDirectory', '.','AdditionalPaths', {'/nfs/matlabuser'}, 'Pool', n_worker(n));
    wait(j)
    r(n,1) = n_worker(n);
    r(n,2) = j.fetchOutputs{1};% Get results into a cell array
    delete(j)
end

%%
N = 100000
a = zeros(N,1); 
parfor I = 1:N 
    a(I) = max(eig(rand(M)));
end
%%