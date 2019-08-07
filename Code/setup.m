%% SETUP

% run this once

listOfFolders = dir; 
for i=1:length(listOfFolders)
    if listOfFolders(i).isdir == 1 && ~strcmp(listOfFolders(i).name, '.') && ~strcmp(listOfFolders(i).name, '..') && ~strcmp(listOfFolders(i).name, '.vscode')
        addpath(genpath(listOfFolders(i).name))

    end
end
