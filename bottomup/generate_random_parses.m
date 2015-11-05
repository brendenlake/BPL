% Mostly bottom-up method for generating a set
% of good candidate parses to optimize further.
%
% Input
%  I: [n x n bool] binary image
%  lib: library
%  ninit: number of parses to choose
%  verbose: true/false
%
% Output
%  bestMP: [ninit x 1 cell] list of best candidate motor programs
%  score_sorted: [n x 1] score of all parses considered in sorted order
%
function [bestMP,score_sorted] = generate_random_parses(I,lib,ninit,verbose)
    if ~exist('verbose','var')
       verbose = false; 
    end    

    % Check that image is in the right format    
    assert(UtilImage.check_black_is_true(I));
    
    % If we have a blank image
    if sum(sum(I))==0
       bestMP = [];
       return
    end
    
    % Get character skeleton from the fast bottom-up method
    G = extract_skeleton(I,verbose);
    
    % Create a set of random parses through random walks
    ps = defaultps_bottomup;
    if verbose, fprintf(1,'\ngenerating random walks...\n'); end
    RW = RandomWalker(G);
    PP = ProcessParses(I,lib,verbose);
    
    % Add deterministic minimum angle walks
    for i=1:ps.nwalk_det
        PP.add(RW.det_walk);
    end
    
    % Sample random walks until we reach capacity.
    walk_count = PP.nwalks;
    while (PP.nl < ps.max_nstroke) && (walk_count < ps.max_nwalk)
        list_walks = RW.sample(1);
        PP.add(list_walks{1});
        walk_count = walk_count + 1;
        if verbose && mod(walk_count,5)==0
            fprintf(1,'%d,',walk_count);
            if mod(walk_count,20)==0
               fprintf(1,'\n'); 
            end
        end
    end
    if verbose, fprintf(1,'done.\n'); end
        
    % Turn parses into motor programs with additional search 
    % and processing
    [bestMP,score_sorted] = parses_to_MPs(I,PP,ninit,lib,verbose);
    
end