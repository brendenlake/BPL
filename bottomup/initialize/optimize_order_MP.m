%
% Optimize the order of strokes in a MotorProgram,
% before the relations are fixed
%
% Input
%  Q: reference to MotorProgram object
%
function optimize_order_MP(Q,lib)

    ps = defaultps_bottomup;

    % get all the permutations to try
    if Q.ns <= ps.max_ns_all_perm
        % try all possible combinations
        P = perms(1:Q.ns);
    else
        % try a subset of the permutations
        np = factorial(ps.max_ns_all_perm); %720        
        P = zeros(np,Q.ns);
        for i=1:np
            P(i,:) = randperm(Q.ns);
        end
        P = unique(P,'rows');
    end
        
    % score all the permutations    
    n = size(P,1);
    scores = zeros(n,1);
    for i=1:n
       QQ = Q.copy(); 
       perm = P(i,:);       
       QQ.S = QQ.S(perm);       
       scores(i) =  scoreMP_NoRel(QQ,lib,'type',true,'token',true,'stat',false,'image',false);
    end

    % pick the best order, and return Q without instantiating relations
    [~,windx] = randmax(scores);
    perm = P(windx,:);
    Q.S = Q.S(perm);
end