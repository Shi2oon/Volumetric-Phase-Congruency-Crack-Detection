function [DisplacementField1_C,VPC1,Ux1,Uy1,Uz1,DU]=...
        LoadData_para(PC1,PC2,Med_Circle,Med_Thresh,Dir,MainDispl,SliceDir)
temp1=importdata(Dir);
try
    D = temp1.data;
catch
    D=temp1;
    D(:,7) = ones(length(D),1);
end
if size(D,2)<=6
    [file,folder,ext]      = fileparts(Dir);
    [x,y,z,disX,disY,disZ] = filcompile (file,[folder ext]) ; 
    D  = [x(:) y(:) z(:) disX(:) disY(:) disZ(:)];
end
 
id = D(:,7)==0;
D(id,[4,5,6]) = nan;
[X, Y, Z, Ux1, Uy1, Uz1] = V2D(D);
% Ux1(Ux1==0) = nan;    Uy1(Uy1==0) = nan;      Uz1(Uz1==0) = nan;

switch MainDispl
    case 'Ux'
        DisplacementField1 = Ux1(:,:,:);
    case 'Uy'
        DisplacementField1 = Uy1(:,:,:);
    case 'Uz'
        DisplacementField1 = Uz1(:,:,:);
end

xaxis = 1:size(DisplacementField1,1);
yaxis = 1:size(DisplacementField1,2);
zaxis = 1:size(DisplacementField1,3);

switch SliceDir
    case 'Ux'
        for i = xaxis
             M1 = squeeze(DisplacementField1(i,:,:));
        % Get the min and max value of Displacement field
        % To get an harmonious colorbar for all graphs
        %    
        %     cmin = min(M1(~isnan(M1(:))));
        %     if ~isnan(cmin)
        %         cmax = max(M1(~isnan(M1(:))));
        %     else
        %         cmin = 0;
        %         cmax = 1;
        %     end
        %         
        %         %DEBUG
        %         subplot(1,5,1)
        %         imagesc(M1)
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement'])
        %         colorbar
        %         caxis([cmin cmax])
            % Module 1
            % Extract itterative slice, Cropped Volume
            M2(i,:,:) =inpaint_nans(M1,1);  %Replace the NaN value using interpolation 
        %         %DEBUG
        %         subplot(1,5,2)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement with NaN'])
        %         colorbar
        %         caxis([cmin cmax])    
            M_Del(i,:,:) = OutDel(squeeze(M2(i,:,:)),Med_Circle,Med_Thresh); %20, 0.001 %Smooth the Displacement field
        %         %DEBUG
        %         subplot(1,5,3)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement smoothed'])
        %         colorbar
        %         caxis([cmin cmax])
            M(i,:,:) =inpaint_nans(squeeze(M_Del(i,:,:)),1);    %Replace the NaN value using interpolation
        %         %DEBUG
        %         subplot(1,5,4)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement field'])
        %         colorbar
        %         caxis([cmin cmax])
            VPC1(i,:,:) = phasecongmono(squeeze(M(i,:,:)), 6,...
                PC1, 2, 0.55, 1.0,0.25,10,PC2,-1);
            %h=imshow(squeeze(VPC1(:,i,:)),[]); hold off;
        %         %DEBUG
        %         subplot(1,5,5)
        %         imagesc(squeeze(VPC1(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' PhaseCongruency'])
        %         colorbar
        end
        close all
        DisplacementField1_C = M;
    case 'Uy'
        for i = yaxis
             M1 = squeeze(DisplacementField1(:,i,:));
        % Get the min and max value of Displacement field
        % To get an harmonious colorbar for all graphs
        %    
        %     cmin = min(M1(~isnan(M1(:))));
        %     if ~isnan(cmin)
        %         cmax = max(M1(~isnan(M1(:))));
        %     else
        %         cmin = 0;
        %         cmax = 1;
        %     end
        %         
        %         %DEBUG
        %         subplot(1,5,1)
        %         imagesc(M1)
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement'])
        %         colorbar
        %         caxis([cmin cmax])
            % Module 1
            % Extract itterative slice, Cropped Volume
            M2(:,i,:) =inpaint_nans(M1,1);  %Replace the NaN value using interpolation 
        %         %DEBUG
        %         subplot(1,5,2)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement with NaN'])
        %         colorbar
        %         caxis([cmin cmax])    
            M_Del(:,i,:) = OutDel(squeeze(M2(:,i,:)),Med_Circle,Med_Thresh); %20, 0.001 %Smooth the Displacement field
        %         %DEBUG
        %         subplot(1,5,3)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement smoothed'])
        %         colorbar
        %         caxis([cmin cmax])
            M(:,i,:) =inpaint_nans(squeeze(M_Del(:,i,:)),1);    %Replace the NaN value using interpolation
        %         %DEBUG
        %         subplot(1,5,4)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement field'])
        %         colorbar
        %         caxis([cmin cmax])
            VPC1(:,i,:) = phasecongmono(squeeze(M(:,i,:)), 6,...
                PC1, 2, 0.55, 1.0,0.25,10,PC2,-1);
            %h=imshow(squeeze(VPC1(:,i,:)),[]); hold off;
        %         %DEBUG
        %         subplot(1,5,5)
        %         imagesc(squeeze(VPC1(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' PhaseCongruency'])
        %         colorbar
        end
        close all
        DisplacementField1_C = M;
    case 'Uz'
        for i = zaxis
             M1 = squeeze(DisplacementField1(:,:,i));
        % Get the min and max value of Displacement field
        % To get an harmonious colorbar for all graphs
        %    
        %     cmin = min(M1(~isnan(M1(:))));
        %     if ~isnan(cmin)
        %         cmax = max(M1(~isnan(M1(:))));
        %     else
        %         cmin = 0;
        %         cmax = 1;
        %     end
        %         
        %         %DEBUG
        %         subplot(1,5,1)
        %         imagesc(M1)
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement'])
        %         colorbar
        %         caxis([cmin cmax])
            % Module 1
            % Extract itterative slice, Cropped Volume
            M2(:,:,i) =inpaint_nans(M1,1);  %Replace the NaN value using interpolation 
        %         %DEBUG
        %         subplot(1,5,2)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement with NaN'])
        %         colorbar
        %         caxis([cmin cmax])    
            M_Del(:,:,i) = OutDel(squeeze(M2(:,:,i)),Med_Circle,Med_Thresh); %20, 0.001 %Smooth the Displacement field
        %         %DEBUG
        %         subplot(1,5,3)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement smoothed'])
        %         colorbar
        %         caxis([cmin cmax])
            M(:,:,i) =inpaint_nans(squeeze(M_Del(:,:,i)),1);    %Replace the NaN value using interpolation
        %         %DEBUG
        %         subplot(1,5,4)
        %         imagesc(squeeze(M2(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' z-displacement field'])
        %         colorbar
        %         caxis([cmin cmax])
            VPC1(:,:,i) = phasecongmono(squeeze(M(:,:,i)), 6,...
                PC1, 2, 0.55, 1.0,0.25,10,PC2,-1);
            %h=imshow(squeeze(VPC1(:,i,:)),[]); hold off;
        %         %DEBUG
        %         subplot(1,5,5)
        %         imagesc(squeeze(VPC1(:,i,:)))
        %         axis image
        %         xlabel('z')
        %         ylabel('x')
        %         title(['y-slice ' num2str(i) ' PhaseCongruency'])
        %         colorbar
        end
        close all
        DisplacementField1_C = M;         
end

DU.RawDisp = DisplacementField1_C;          DU.VPC = VPC1;   
DU.Ux = Ux1;            DU.Uy  = Uy1;       DU.Uz = Uz1;
DU.X  = X;              DU.Y   = Y;         DU.Z  = Z;
DU.MainDispl = MainDispl;
end