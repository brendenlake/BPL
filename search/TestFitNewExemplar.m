%
% Test the optimization algorithm that fits the parse
%  from one image of a character to another image "token" of that character
% 
nsamp = 10; % total number of MCMC samples
nsamp_mcmc = 10; % number of samples to keep
auto_affine = true; % start optimization after affine fitting
ns = 1; % number of strokes

ps = defaultps;
load(ps.libname,'lib');
motor_program = generate_character(lib,ns);
M1 = motor_program();        
M2 = motor_program();
while M1.ink_off_page || M2.ink_off_page   
   motor_program = generate_character(lib,ns);
   M1 = motor_program();        
   M2 = motor_program();      
end
M1.I = M1.pimg > 0.5;
M2.I = M2.pimg > 0.5;

% visualize
figure(1);
clf
subplot(1,3,1);
plot_motor_to_image(M1.I,M1.motor_warped);
title(['image score ',num2str(scoreMP(M1,lib,'image',true,'type',false,'token',false),3)]);
xlabel('example 1');
subplot(1,3,2);
plot_motor_to_image(M2.I,M2.motor_warped);
title(['image score ',num2str(scoreMP(M2,lib,'image',true,'type',false,'token',false),3)]);
xlabel('example 2');
pause(0.1);
drawnow

% compute samples for type-level variables
all_samples = mcmc_all(M1,lib,nsamp_mcmc,'type');
indx = round(linspace(nsamp,nsamp_mcmc,nsamp));
subset_samples = all_samples(indx);

% run the optimization
Mfit2 = FitNewExemplar(M2.I,subset_samples,lib,auto_affine);

subplot(1,3,3);
plot_motor_to_image(Mfit2.I,Mfit2.motor_warped);
title(['image score ',num2str(scoreMP(Mfit2,lib,'image',true,'type',false,'token',false),3)]);
xlabel('ex 1 fit to ex 2');