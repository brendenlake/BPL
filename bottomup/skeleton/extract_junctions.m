% Extract critical features (junctions and endpoints).
%
% Rougly follows procedure of Liu et al. "Identification of Fork Points..." IEEE
% TAPMI
%
% For a circular stroke, there are no "features" by this definition.
% Thus, for any partition of the pixels into regions, there should
% be at least one feature for the tracing algorithm to find.
%
% Input
%  T: [n x n boolean] thinned image.
%    images are binary, where true means "black"
%
% Output
%  SN: [n x n boolean] extracted features.
function SN = extract_junctions(T)

    SE = bwmorph(T,'endpoints');
    SB = T; % black pixels
    sz = size(T,1);
    
    lutS3 = makelut( @(P)fS3(P) , 3);
    S3 = applylut(T,lutS3);

    % final criteria
    SN = SE | (SB & S3);
    
    % Check to see that each connected component has a feature.
    % This is necessary to process circles in the image.
    CC = bwconncomp(T,8);
    nCC = CC.NumObjects;
    for c=1:nCC
       
       pid = CC.PixelIdxList{c};       
       
       % We have a circle. Circles are generally drawn from the
       % top, we choose the top pixel here
       if sum(SN(pid))==0
          [irow,icol] = ind2sub(sz,pid);
          sel = argmin(irow);
          SN(pid(sel)) = true;           
       end
       
    end
       
end

% See Liu et al.
function Y=fS3(P)
    sz = size(P);
    assert(isequal(sz,[3 3]));
    
    % Get cross number
    NC = fNC(P);
    
    % Count black pixels
    PM = P;
    PM(2,2) = false;
    NB = sum(PM(:));
    
    % Criteria
    Y = (NC >= 3-eps) || (NB >= 4-eps);
end

% See Liu et al.
function Y=fNC(P)       
    sum = 0;
    for i=0:7
        sum = sum + abs( P(fIP(i+1)) - P(fIP(i)) );
    end
    Y = sum./2;
end

% See Liu et al.
function newlindx = fIP(lindx)
    switch lindx
        case {0,8}
            i=1; j=2;
        case 1
            i=1; j=3;
        case 2
            i=2; j=3;
        case 3
            i=3; j=3;
        case 4
            i=3; j=2;
        case 5
            i=3; j=1;
        case 6
            i=2; j=1;
        case 7
            i=1; j=1;
    end
    newlindx = sub2ind([3 3],i,j);
end