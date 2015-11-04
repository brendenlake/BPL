%
% Visualize the outcome of the skeleton extraction 
% process for a set of characters.
%
nimg = 9;
load('data_background_processed','D');

nrow = ceil(sqrt(nimg));
figure(1)
clf;
list_I = cell(nrow,1);
list_T = cell(nrow,1);
for i=1:nimg
    I = D.get('image','random','random','random');
    T = extract_skeleton(I);
    subplot(nrow,nrow,i);
    viz_skel.plot_trace(I,T);
    fprintf(1,'image %d\n',i);
    
    list_I{i} = I;
    list_T{i} = T;
end