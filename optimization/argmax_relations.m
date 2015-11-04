function ll = argmax_relations(lib,M,all_R,stroke_num)
% OPTIMIZE_RELATIONS
%
% For independent relations, it assumes that the global position
%   is locked to the token-level position
% For "attached" relations, it approximates the continuous
%   parameter by the discrete "motor" evaluations along the spline
%
% Input
%  M : Motor PRogram
%  all_R (optional): struct containing output of cache_enum_all_relations
%  stroke_num: vector of stroke indices that 
%              we want to optimize the relations for
%
% Output
%  ll : [length(stroke_num) x 1] all log-CPDs that depend directly
%    on each relation
%
    if ~exist('stroke_num','var')
       stroke_num = 1:M.ns;
    end
    
    PM = M.parameters;
    ll = zeros(length(stroke_num),1);
    for indx=1:length(stroke_num)
        s = stroke_num(indx);
        
        crit_pos = M.S{s}.pos_token;
        prev = 1:s-1;
        
        % caching to pre-compute the list of all relations,
        % and their prior under the type model
        if ~exist('all_R','var') || isempty(all_R)
            [list_R,list_score_type] = enum_all_relations(M.S(prev),lib);
        else
            list_R = all_R{s}.list_R;
            list_score_type = all_R{s}.list_score_type;
        end
        
        nR = length(list_R);
        scores = zeros(nR,1);
        
        % approximate the score for each of the relations
        for i=1:nR
            llextra = 0;
            switch list_R{i}.type
                case 'unihist'
                    gpos = crit_pos;
                    
                    % ensure that the global position is in bounds,
                    % otherwise we will get -inf score
                    if gpos(1) < 0, gpos(1) = 0; end
                    if gpos(1) > PM.imsize(1), gpos(1) = PM.imsize(1); end
                    if gpos(2) < -PM.imsize(2), gpos(2) = -PM.imsize(2); end
                    if gpos(2) > 0, gpos(2) = 0; end
                    
                    list_R{i}.gpos = gpos;
                case 'mid'
                    [~,eval_spot] = optimize_along_spline(crit_pos,list_R{i},M,lib);
                    list_R{i}.eval_spot_type = [];
                    list_R{i}.eval_spot_token = eval_spot;
                    llextra = CPD.score_relation_token_approx_marginalize(lib,eval_spot);
                                      
                otherwise
            end
            
            llprior = list_score_type(i);
            if isnan(llprior)
                llprior = CPD.score_relation_type(lib,list_R{i});
                assert(strcmp(list_R{i}.type,'unihist'));
            end            
            llfit = CPD.score_position(lib,crit_pos,list_R{i},M.S(prev));
            scores(i) = llprior + llfit + llextra;            
        end

        [ll(indx),windx] = max(scores);
        M.S{s}.R = list_R{windx};
    end

end

%
% Find the "evaluation spot" for a "mid" relation that
% is the most likely.
%
function [ll,eval_spot] = optimize_along_spline(crit_pos,R,M,lib)    
    myspline_eval = M.S{R.attach_spot}.motor{R.subid_spot};
    neval = size(myspline_eval,1);
    seval = bspline_gen_s(R.ncpt,neval);    
    dist = distfun(crit_pos,myspline_eval,lib);
    [mindist,windx] = min(dist);
    ll = -mindist - log(lib.rel.sigma_x) - log(lib.rel.sigma_y) - log(2*pi);
    eval_spot = seval(windx);
end

function dist = distfun(x,Y,lib)
    assert(isvector(x));
    x = x(:)';
    n = size(Y,1);
    dist = zeros(n,1);
    for i=1:n
        y = Y(i,:);
        dist(i) = 1./(2*lib.rel.sigma_x^2) .*(x(1)-y(1)).^2 ...
            + 1./(2*lib.rel.sigma_y^2) .*(x(2)-y(2)).^2;
    end
end