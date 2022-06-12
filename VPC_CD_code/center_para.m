function center=center_para(RawDisp,indexSD,Ncrackedslice)
switch indexSD
    case 1
        MX=squeeze(RawDisp(Ncrackedslice,:,:));
        H1=imagesc(MX);axis image; colorbar
        set(gcf,'position',[600,100,1000,750])
        title ('Choose a point inside the crack')
        [ymask] = ginput(1);
        ymask=round(ymask);
        center=[Ncrackedslice ymask(2) ymask(1)];
    case 2
        MX=squeeze(RawDisp(:,Ncrackedslice,:));
        H1=imagesc(MX);axis image; colorbar
        set(gcf,'position',[600,100,1000,750])
        title ('Choose a point inside the crack')
        [ymask] = ginput(1);
        ymask=round(ymask);
        center=[ymask(2) Ncrackedslice ymask(1)];
    case 3
        MX=squeeze(RawDisp(:,:,Ncrackedslice));
        H1=imagesc(MX);axis image; colorbar
        set(gcf,'position',[600,100,1000,750])
        title ('Choose a point inside the crack')
        [ymask] = ginput(1);
        ymask=round(ymask);
        center=[ymask(2) ymask(1) Ncrackedslice];
end 
end