% use this function is the crack is in a notched sample. the slices need to
% start from the crack tip rather than from sides.

function DO = ANG(filC,Kk)
close all;
disp('We will need to re-slice the 3D DVC by angles that starts from the notch');
disp('corner to be able to calculate the stress intesity factor correctly'); 
dot = importdata(filC);
DU.X   = reshape(dot(:,1),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
DU.Y   = reshape(dot(:,2),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
DU.Z   = reshape(dot(:,3),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
DU.Ux  = reshape(dot(:,4),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
DU.Uy  = reshape(dot(:,5),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
DU.Uz  = reshape(dot(:,6),length(unique(dot(:,1))),length(unique(dot(:,2))),length(unique(dot(:,3))));
contourf(squeeze(DU.X(:,:,Kk)),squeeze(DU.Y(:,:,Kk)),squeeze(DU.Uz(:,:,Kk))); axis image
DU.X = unique(DU.X);        DU.Y = unique(DU.Y);    DU.Z = unique(DU.Z);    
set(gcf,'position',[600,100,1000,750]);         colorbar; colormap jet
title('Click on the Corner where the cracks starts')
%         [xc,yc] = ginput(1);    title '';
        xc=4.208;       yc = 0.3729;
        if  abs(xc-max(DU.X)) > abs(xc-min(DU.X)) 
            x(1) = min(DU.X);        x(2) = max(DU.X);    
        else
            x(1) = max(DU.X);        x(2) = min(DU.X);    
        end
        if  abs(yc-max(DU.Y)) > abs(yc-min(DU.Y))
            y(1) = min(DU.Y);        y(2) = max(DU.Y);   
        else
            y(1) = max(DU.Y);        y(2) = min(DU.Y);    
        end

%%
Diag = rad2deg(atan((y(2)-y(1))/(x(2)-x(1))));
for Depth = 1:size(DU.Uz,3)
    slices = 0;           
    for vi = 1:91/length(DU.Y):91
        theta  = (vi-1)*sign(Diag);
        slices = slices+1;
        if      theta < Diag && sign(Diag) == 1% change if the crack is positive corner 
                Yc(Depth,slices,2) = (tan(deg2rad(theta))*(x(2)-x(1)))+y(1);
                Xc(Depth,slices,2) = x(2);
        elseif  theta > Diag && sign(Diag) == 1 % change if the crack is positive corner 
                Xc(Depth,slices,2) = (y(2)-y(1))/tan(deg2rad(theta))+x(1);
                Yc(Depth,slices,2) = y(2);        
        elseif  theta > Diag && sign(Diag) == -1 % change if the crack is positive corner 
                Yc(Depth,slices,2) = (tan(deg2rad(theta))*(x(2)-x(1)))+y(1);
                Xc(Depth,slices,2) = x(2);
        elseif  theta ==Diag
                Xc(Depth,slices,2) = x(2);          Yc(Depth,slices,2) = y(2);
        elseif  theta < Diag && sign(Diag) == -1 % change if the crack is positive corner 
                Xc(Depth,slices,2) = (y(2)-y(1))/tan(deg2rad(theta))+x(1);
                Yc(Depth,slices,2) = y(2);
        end
                Xc(Depth,slices,1) = x(1);          Yc(Depth,slices,1) = y(1);
                if Depth == round(size(DU.Uz,1)/2,0)
hold on; plot(squeeze(Xc(Depth,slices,:)),squeeze(Yc(Depth,slices,:)),'--w');  
                end
%%           colums row stacks and get data on the line
        [cx,cy,~] = improfile(squeeze(DU.Uz(:,:,Depth)),...
            squeeze(Xc(Depth,slices,:)),squeeze(Yc(Depth,slices,:)),length(DU.X));
             Ux = squeeze(DU.Ux(:,:,Depth));     
             Uy = squeeze(DU.Uy(:,:,Depth));      Uz = squeeze(DU.Uz(:,:,Depth));
        
        % re-arrange
    [~, index] = min(abs(DU.X-Xc(Depth,slices,1)));        Xc(Depth,slices,1) = DU.X(index);
    [~, index] = min(abs(DU.Y-Yc(Depth,slices,1)));        Yc(Depth,slices,1) = DU.Y(index);
    for i=1:length(cx)
        [~, index] = min(abs(DU.X-cx(i)));                 cx(i) = DU.X(index);
        [~, index] = min(abs(DU.Y-cy(i)));                 cy(i) = DU.Y(index);
        [indexy(i)] = ind2sub(size(DU.Y),find(DU.Y==cy(i)));
        [indexx(i)] = ind2sub(size(DU.X),find(DU.X==cx(i)));
        line2(i) = Ux(indexx(i),indexy(i)); 
        line3(i) = Uy(indexx(i),indexy(i));
        line4(i) = Uz(indexx(i),indexy(i)); 
    end
        DO.X(Depth,slices,:) = cx(:); % check how xes are arranged
        DO.Y(Depth,slices,:) = cy(:); % check how yes are arranged
        cz = ones(size(cy)).*DU.Z(Depth);
        DO.Z(Depth,slices,:) = cz(:);
        DO.Ux(Depth,slices,:)      = line2(:);
        DO.Uy(Depth,slices,:)      = line3(:);
        DO.Uz(Depth,slices,:)      = line4(:);
        
        if Depth == round(size(DU,1)/2,0)
            hold on;    plot(cx,cy,'-k');  
        end
    end
end
hold off;   

    saveas(gcf,[fileparts(filC) '\Angled_DVC.fig']);
    saveas(gcf,[fileparts(filC) '\Angled_DVC.tif']); close
 
T2 = table(DO.X(:),DO.Y(:),DO.Z(:),DO.Ux(:),DO.Uy(:),DO.Uz(:));
writetable(T2,[fileparts(filC) '\angeled.txt'],'Delimiter','\t','WriteRowNames',true);
