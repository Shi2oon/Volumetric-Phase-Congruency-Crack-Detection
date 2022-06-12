%%
clear all; close all; clc
addpath('data'); 
% load head_ct; %V; 
load PCT
V =PCT;
%%
propagation_weight = 0.02; 
GAC_weight = .02; 
% g = ones(size(V)); % linear diffusion 
g = ac_gradient_map(V,1); 
delta_t = 1; 
mu = 0.2; 
fileHeader='B';
numZeros=5;
margin = 1; 
center = [43, 25, 20]; 
phi = zeros(size(V)); 
phi(center(1)-margin:center(1)+margin,...
    center(2)-margin:center(2)+margin,...
    center(3)-margin:center(3)+margin) = 1; 
%%
for i = 1:100
    figure(1);
    
    phi = ac_hybrid_model(V-mu, phi-.5, propagation_weight, GAC_weight, g, ...
        delta_t, 1); 
    if exist('h','var') && all(ishandle(h)), delete(h); end
    iso = isosurface(phi,0);
    h = patch(iso,'edgecolor','r','facecolor','w');  axis equal;  view(3); 
    set(gcf,'name', sprintf('#iters = %d',i));
    sprint_im=['%0' int2str(numZeros) 'd'];
%     fileName1=[fileHeader sprintf(sprint_im,i) '.tiff'];
 fileName1=[fileHeader sprintf(sprint_im,i)];
    drawnow; 
   
%     print('-f1',fileName1,'-dpng')
%     imwrite(f)
end

%%
figure;
slice = [1:6:60];
for i = 1:10
    subplot(2,5,i); imshow(V(:,:,slice(i)),[]); hold on; 
    c = contours(phi(:,:,slice(i)),[0,0]);
    zy_plot_contours(c,'linewidth',2);
end

Seg3D = zeros(size(V));
Seg3D(phi>0)=1;

