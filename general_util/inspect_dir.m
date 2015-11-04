%
% Returns cell array of file names in 
% a directory
%
function names = inspect_dir(dir_name)
    listing = dir(dir_name);
    if size(listing,1)==0
       error('No such directory'); 
    end
    len = length(listing);
    except = {'.','..','.DS_Store'}; 
    names = [];
    for i=1:len
        this_name = listing(i).name;
        if ~any(strcmp(this_name,except))
            names = [names; {this_name}];
        end        
    end
end