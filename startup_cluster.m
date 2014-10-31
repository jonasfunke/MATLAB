%% cluster startup file
disp('Welcome to the High Performance Computing Cluster')

base_dir = '/Users/jonasfunke/Documents/Temporary/testfolder'; % base directory for matlabuser
excluded_directories = {'otherscripts', 'root'}; % diresctories that do not belong to a user


folders = dir(base_dir);
folders = folders([folders.isdir]); % get only folders
folders(strncmp({folders.name}, '.', 1)) = []; % remove folders that start with a '.'

% exclude direcotries that do not belong to a user
for i=1:length(excluded_dirs)
    folders(strncmp({folders.name}, excluded_dirs{i}, 1)) = []; % remove folders that start with a '.'
end

%% select a user
[selection, ok] = listdlg('PromptString','Select a user:', 'SelectionMode','single', 'ListString',{folders.name});
            
if ok
    selected_dir = [base_dir filesep folders(selection).name filesep 'MATLAB'];
    disp(['Selected MATLAB directory: ' selected_dir])
    
    % change to selected directory
    try
        cd(selected_dir)
        userpath(selected_dir) % set new userpath
        disp(['Userpath set to: ' selected_dir ])
    catch err
        rethrow(err)
        if strcmp(err.identifier, 'MATLAB:cd:NonExistentDirectory')
            disp('No MATLAB directory present. Create a directory called MATLAB!')
        end
    end
else
    disp('No directory selected: you need to start manually')
end
    