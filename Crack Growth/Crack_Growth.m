function [Crack] = Crack_Growth(Xdata,Ydata,Zdata,Timepoints,Di)
%Script to run crackfront program
%Ensure matlab data file is loaded, and mycolormap.mat is in the path.
% clc; clear
% load( [pwd '\input.mat'])
close all; clc
[Points,PointsSize]=setuprealall(Xdata,Ydata,Zdata,Timepoints);
% setuprealdemo;
Crack.data            = CrackFront8(Points, PointsSize, Timepoints, 4);
Crack.Points          = Crack.data{1}; 
Crack.growthvector    = Crack.data{2};
Crack.growth2d        = Crack.data{3};
Crack.length2d        = Crack.data{4};
Crack.growth2dnorm    = Crack.data{5};
Crack.longest         = Crack.data{6};
Crack.fastest         = Crack.data{7};
Crack.fastestperstep  = Crack.data{8};
Crack.meangrowth      = Crack.data{9};
Crack.meanlength      = Crack.data{10};

%%
mkdir(Di);      save([Di '\Data.mat'],'Crack');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  set(gcf,'position',[100,50,1400,920]); 
  savefig(FigHandle, fullfile(Di, [num2str(iFig), '.fig']));  
  saveas(FigHandle, fullfile(Di, [num2str(iFig), '.tif'])); close
end

%%
for iv=1:length(Crack.meangrowth);   Crack_Growth(iv) = sum(Crack.meangrowth(1:iv));   end
for iv=1:length(Crack.meanlength);   Crack_Length(iv) = sum(Crack.meanlength(1:iv));   end
Crack_Growth = [0 Crack_Growth];     Crack_Length     = [0 Crack_Length];
close;      plot(Timepoints, Crack_Growth,'--ok','MarkerFaceColor','k'); 
set(gcf,'position',[10 50 1100 950]);   box off;    
title('Mean growth rates'); xlabel('Cycles [K]');   ylabel('Growth rates')
saveas(gcf,[Di '\Crack Growth.fig']);   saveas(gcf,[Di '\Crack Growth.tif']); close
close;      plot(Timepoints, Crack_Length,'--ok','MarkerFaceColor','k');
set(gcf,'position',[10 50 1100 950]);   box off;    
xlabel('Cycles [K]');   ylabel('Length [mm]')
saveas(gcf,[Di '\Crack Length.fig']);   saveas(gcf,[Di '\Crack Length.tif']); close