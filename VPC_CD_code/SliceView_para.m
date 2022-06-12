function SliceView=SliceView_para(RawDisp,indexSD)
siz = size(RawDisp);
switch indexSD
    case 1
        for j=1:siz(1)
            imagesc(squeeze(RawDisp(j,:,:)));colorbar;axis equal
            title(['Slice number ', num2str(j), ' (Note the number of a slice with a visible crack)']);
            pause
        end
    case 2
        for j=1:siz(2)
            imagesc(squeeze(RawDisp(:,j,:)));colorbar;axis equal
            title(['Slice number ', num2str(j), ' (Note the number of a slice with a visible crack)']);
            pause
        end
    case 3
        for j=1:siz(3)
            imagesc(squeeze(RawDisp(:,:,j)));colorbar;axis equal
            title(['Slice number ', num2str(j), ' (Note the number of a slice with a visible crack)']);
            pause
        end
end
end