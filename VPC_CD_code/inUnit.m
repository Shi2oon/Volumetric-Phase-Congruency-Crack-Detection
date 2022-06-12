function [DU] = inUnit(DU,PixelSize,Overlap,Subset,Pixel_unit)
DU.Pixel2um = PixelSize*(1-0.01*Overlap)*Subset;
datum = DU.X.*DU.Pixel2um;
    
    %% Units + 
switch Pixel_unit
    case 'mm';              datum = datum.*1e-3;     offest2 = 1e-3;  % conver to m
    case '\mum';            datum = datum.*1e-6;     offest2 = 1e-6;  % conver to m
end

% and for some reason Abaqus is really bad in handling small dim.
if      mean(abs(datum(:,1)))<5e-5;     offset1 = 1e-6;       Pixel_unit = '\mum';
elseif  mean(abs(datum(:,1)))<5e-2;     offset1 = 1e-3;       Pixel_unit = 'mm'; 
else;   offset1 = 1;                     Pixel_unit = 'm';                end

DU.Pixel2um = DU.Pixel2um*offest2/offset1;
DU.Pixel_unit = Pixel_unit;

%% convert pixel to physical   
    DU.Phy.RawDisp  = DU.RawDisp.*DU.Pixel2um;
    DU.Phy.VPC      = DU.VPC.*DU.Pixel2um;
    DU.Phy.Ux       = DU.Ux.*DU.Pixel2um;
    DU.Phy.Uy       = DU.Uy.*DU.Pixel2um;
    DU.Phy.Uz       = DU.Uz.*DU.Pixel2um;
    DU.Phy.X        = DU.X.*DU.Pixel2um;
    DU.Phy.Y        = DU.Y.*DU.Pixel2um;
    DU.Phy.Z        = DU.Z.*DU.Pixel2um;
    DU.Phy.units    = DU.Pixel_unit;