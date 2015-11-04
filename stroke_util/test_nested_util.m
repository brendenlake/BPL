% Test manipulation of nested cell arrays

% Test nested extraction
A = {{1 {2 3}} {4 5} {6 7 8} {9 10}};
C_real = vec(1:10);
C = flatten_nested(A);
for i=length(C_real)
   if C{i}~=C_real(i) 
      fprintf(1,'flatten_nested test failed\n');
      assert(false);
   end
end
fprintf(1,'flatten_nested test passed\n');

% Test nested element-wise modification
fnc = @(x)x + 1;
A = apply_to_nested(A,fnc);
C2 = flatten_nested(A);
for i=length(C_real)
   if C2{i}~=C_real(i)+1
      fprintf(1,'apply_to_nested test failed\n');
      assert(false);
   end
end
fprintf(1,'apply_to_nested test passed\n');

% MULTI-NESTED
% Test nested element-wise modification
A = { {{1} {2}} {{3} {4}} {{5} {6}} };
B = { {{-1} {-2}} {{-3} {-4}} {{-5} {-6}} };
C_real = zeros(6,1);
fnc = @(x,y) x + y;
A = multi_apply_to_nested(fnc,A,B);
C3 = flatten_nested(A);
for i=length(C_real)
   if C3{i}~=C_real(i)
      fprintf(1,'multi_apply_to_nested test failed\n');
      assert(false);
   end
end
fprintf(1,'multi_apply_to_nested test passed\n');

% MULTI-NESTED
% Test nested element-wise modification
A = { {{1} {2}} {{3} {4}} {{5} {6}} };
B = { {{-1} {-2}} {{-3} {-4}} {{-5} {-6}} };
C = { {{5} {10}} {{5} {10}} {{5} {12}} };
C_real = [5 10 5 10 5 12];
fnc = @(x,y,z) x + y + z;
A = multi_apply_to_nested(fnc,A,B,C);
C4 = flatten_nested(A);
for i=length(C_real)
   if C4{i}~=C_real(i)
      fprintf(1,'multi_apply_to_nested 3-arg test failed\n');
      assert(false);
   end
end
fprintf(1,'multi_apply_to_nested 3-arg test passed\n');