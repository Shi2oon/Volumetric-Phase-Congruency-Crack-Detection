function [Ncrackedslice]=SelectSlice(RawDisp,indexSD,Ncrackedslice)
% close all;

if Ncrackedslice == 999
    while Ncrackedslice == 999 
        H1 = figure;
        for iv =1: size(RawDisp,indexSD)
switch indexSD
    case 1
        MX=squeeze(RawDisp(iv,:,:));
        H1=imagesc(MX);     axis image;     colorbar
        set(gcf,'position',[600,100,1000,750])
        title (['Slice number ' num2str(iv) ' out of ' num2str(size(RawDisp,1))])
    case 2
        MX=squeeze(RawDisp(:,iv,:));
        imagesc(MX);     axis image;     colorbar
        set(gcf,'position',[600,100,1000,750])
        title (['Slice number ' num2str(iv) ' out of ' num2str(size(RawDisp,2))])
    case 3
        MX=squeeze(RawDisp(:,:,iv));
        imagesc(MX);     axis image;     colorbar
        set(gcf,'position',[600,100,1000,750])
        title (['Slice number ' num2str(iv) ' out of ' num2str(size(RawDisp,3))])
end 
%         waitforbuttonpress; 
        pause(0.25)
        end
    close 
    Ncrackedslice = input('So what is cracked slice number now! or type 999 : ');
    end
end 
if Ncrackedslice == 1 
    disp('WRONG CHOICE!! You cannot use 1st slice, we will use the 2nd now'); 
    Ncrackedslice = 2;      
end
end