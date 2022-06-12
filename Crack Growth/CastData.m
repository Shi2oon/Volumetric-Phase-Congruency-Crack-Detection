close all;      restoredefaultpath;  addpath(genpath(pwd));              
DS = com.mathworks.mde.desk.MLDesktop.getInstance();   	DS.closeGroup('Variables');
clc;clear;      set(0,'defaultAxesFontSize',35);set(0,'DefaultLineMarkerSize',15)   
Headers   = readtable('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Master_Excel v2.xlsx');
No_Cycles = table2array(Headers(:,5));           COM       = table2array(Headers(:,6));
Direct    = table2array(Headers(:,7));           PC_logic  = table2array(Headers(:,8)); 
%% Nor 1 to 58 .. Up 59 to 66 .. Down 67  to 69
started = 1;		ended = 58;         Zd = NaN(54,26);
counter = 0;        Yd = NaN(54,26);    Xd = NaN(54,26);
for ish = started:ended
  if PC_logic{ish} == 'Y'  && COM(ish) == 1
    Dir.inDir = [Direct{ish} '\Shoed'];
    Dir.Crack   = [Dir.inDir '\PC Analysis Uz'];
    switch Dir.Crack
        case 'P:\Shraddha\Displacement\450N-4500N 247k cycles_DVC_Vec_96x96x96_75ov=unknown\Shoed\PC Analysis Uz'
        otherwise
            load([Dir.Crack '\OutPut.mat'],'C','DU','crkTip','stepY','stepX');
            counter = counter+1;   
            Xdata(counter,1) = {crkTip(:,2).*stepX};
            Xd(1:length(crkTip(:,2)),counter) = squeeze(crkTip(:,2).*stepX);
            Ydata(counter,1) = {crkTip(:,1).*stepY};
            Yd(1:length(crkTip(:,2)),counter) = crkTip(:,1).*stepY;
            Zdata(counter,1) = {ones(size(Ydata{counter,1}))};
            Timepoints(1,counter) = str2num(cell2mat(No_Cycles(ish))).*1000;
            Zd(1:length(crkTip(:,2)),counter) = ones(size(Ydata{counter,1}))....
                                                *Timepoints(1,counter)./1000;
    end
    clearvars -except Xdata counter Timepoints Ydata Zdata ...
                      PC_logic Direct COM ish No_Cycles Xd Yd Zd
  end
end
%%
close all
s1=subplot(1,1,1); 	contourf(Xd,Yd,Zd,'LineStyle','none'); % brighten(0.5);
s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; %colormap jet;
c = colorbar;  c.Label.String = 'Cycles [K]'; xlabel('X [mm]');  ylabel('Y [mm]');
set(gcf,'position',[10 50 1100 950]); box off
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Length & Cycles.fig'); 
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Length & Cycles.tif'); close
%%
cd('C:\Users\scro3511\Downloads\codes')
Main
save('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Ka.mat','Crack','Xd','Yd','Zd','-append');
clear; clc; load('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Ka.mat');

%%
clc; clearvars Crack_Growth Crack_Length
for iv=1:length(Crack.meangrowth);   Crack_Growth(iv) = sum(Crack.meangrowth(1:iv));   end
for iv=1:length(Crack.meanlength);   Crack_Length(iv) = sum(Crack.meanlength(1:iv));   end
Crack_Growth = [0 Crack_Growth];     Crack_Length     = [0 Crack_Length];
NoCycles(isnan(NoCycles))=[];   NoCycles(NoCycles==247)=[];
close;      plot(NoCycles(1:26), Crack_Growth,'--ok','MarkerFaceColor','k'); 
set(gcf,'position',[10 50 1100 950]);   box off;    xlim([0 3e2])
title('Mean growth rates'); xlabel('Cycles [K]');   ylabel('Growth rates')
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Crack Growth.fig'); 
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Crack Growth.tif'); close
close;      plot(NoCycles(1:26), Crack_Length,'--ok','MarkerFaceColor','k');
set(gcf,'position',[10 50 1100 950]);   box off;    xlim([0 3e2])
xlabel('Cycles [K]');   ylabel('Length [mm]')
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Crack Length.fig'); 
saveas(gcf,'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Crack Length.tif'); close

%%
SIF = Ktrue;    SIF(isnan(SIF)) = [];   SIFe = Ketrue;    SIFe(isnan(SIFe)) = [];
SIF(24) = [];   SIFe(24) = [];
close all;  errorbar(SIF,Crack_Growth,SIFe,'horizontal','ok','MarkerSize',18,...
    'MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k'); hold off
KLO = SIF;    KLO(KLO<9.6)=NaN;     Klo = Crack_Growth; Klo(isnan(KLO))=[];
KLO(isnan(KLO))=[]; 
VB = fit(KLO', log10(Klo'),'c*(x)^b');
Klp = fit(KLO', Klo','exp1');
hold on;    Po = plot(Klp); hold off; legend off
Po.LineWidth = 5;
set(gcf,'WindowStyle','normal');    set(gcf,'position',[600,50,1000,920]);  
axis tight;     ylim([1e-6 3.3e-4]);  xlim([8.8 11.7]);
set(gca, 'YScale', 'log')
title({['da/dN = ' num2str(round(-VB.c,0)) '*(\DeltaK)^{' num2str(round(-VB.b,2)) ...
    '}, R = ' num2str(450/4500)]; ''});
ylabel('da/dN [mm/cycle]');       xlabel('\DeltaK [MPa m^{0.5}]');
saveas(gcf, 'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\SIF & Growth.fig');      
saveas(gcf, 'A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\SIF & Growth.tif'); close
save('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Ka.mat',...
    'Crack_Length','Crack_Growth','SIF','-append');
clear; load('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Ka.mat')
save('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Ka.mat')
load('A:\OneDrive - Nexus365\Work\Papers\Cast Iron\Crack Velocity\Ka.mat')