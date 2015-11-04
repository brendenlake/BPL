%
% Generate new tokens from a type using the prior. 
%
ps = defaultps;
nsamp = 9;
lib = loadlib;
motor_program = generate_character(lib);
samples = cell(nsamp,1);
for i=1:nsamp
    samples{i} = motor_program();
end

figure(1)
clf;
nrow = ceil(sqrt(nsamp));
for i=1:nsamp
    sptight(nrow,nrow,i);
    plot_image_only(samples{i}.pimg);
end