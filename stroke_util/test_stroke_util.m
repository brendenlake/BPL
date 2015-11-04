% Test code for the stroke utility procedures

% Fake character
stk1 = [1 2; 1 2; 1 2; 2 2];
stk2 = [1 2];
stk3 = [2 2];
stk4 = [0 0; 20 20];
stk1 = stk1 + randn(size(stk1));
stk2 = stk2 + randn(size(stk2));
stk3 = stk3 + randn(size(stk3));
stk4 = stk4 + randn(size(stk4));

vstk = [stk1; stk2; stk3; stk4];
S = {stk1; stk2; stk3; stk4};
S2 = {2*stk1; 2*stk2; 2*stk3; 2*stk4}; % rescale by 2
S_reverse = {stk1(end:-1:1,:); stk2; stk3; stk4(end:-1:1,:)};
com = mean(vstk);
rg = range(vstk,1);

% Check center of mass
com2 = com_char(S);
if aeq(com,com2)
   fprintf(1,'Center of mass test passed\n');
else
   error('COM failed');
end

% Check range of character
myrange = range_char(S);
if isequal(myrange,rg)
    fprintf(1,'Range test passed\n');
else
    error('Range test failed');
end

% Check removing strokes
minlen = 3;
mindist = 10;
RR = remove_short_stk(S,minlen,mindist);
if isequal(RR,S([1 4]))
   fprintf(1,'Remove stroke test passed\n');
else
   error('Remove stroke test failed');
end

% Check rescaling
R = apply_each_stroke(S,@(x)rescale_stk(x,[2 2]));
if isequal(R,S2)
   fprintf(1,'Apply each stroke test passed\n');
   fprintf(1,'Rescale test passed\n');
else
   error('Rescale test failed');
end

% Check reverse
RV = apply_each_stroke(S,@(x)reverse_stk(x));
if isequal(RV,S_reverse)
   fprintf(1,'Reverse test passed\n');
else
   error('Reverse test failed');
end