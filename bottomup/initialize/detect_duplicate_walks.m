%
% Takes a cell array of random walks, and
% returns a boolean vector of the indices "to remove"
% that are duplicates
%
function tormv = detect_duplicate_walks(walks)

    nwalks = length(walks);
    EQ = false(nwalks);
    for i=1:nwalks
        for j=i+1:nwalks
            EQ(i,j) = isequal_walk(walks{i},walks{j});
        end
    end
    tormv = any(EQ,1);
    
end