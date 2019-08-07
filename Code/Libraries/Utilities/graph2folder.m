function graph2folder(graphpath, folderpath, newname)

% GRAPH2FOLDER moves a graph to a folder. If newname is present, the graph changes name to newname

if nargin<3
    status = copyfile(graphpath, folderpath)
    if status
        disp(['Graph ' graphpath ' copied to ' folderpath]);
    end
else
    destination = [folderpath '/' newname];
    status = copyfile(graphpath, destination)
    if status
        disp(['Graph ' graphpath ' copied to ' folderpath ' with new name ' newname]);
    end

end

    

end
