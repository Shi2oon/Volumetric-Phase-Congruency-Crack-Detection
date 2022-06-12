function DU = re_arrangeData(indexSD,DU)
%% re-arrange data
    switch indexSD
    case 1;     DU.RawDisp = permute(DU.RawDisp, [1 2 3]);
                DU.VPC     = permute(DU.VPC, [1 2 3]);
                DU.Ux      = permute(DU.Ux, [1 2 3]);
                DU.Uy      = permute(DU.Uy, [1 2 3]);
                DU.Uz      = permute(DU.Uz, [1 2 3]);
                DU.X       = permute(DU.X, [1 2 3]);
                DU.Y       = permute(DU.Y, [1 2 3]);
                DU.Z       = permute(DU.Z, [1 2 3]);
            DU.Phy.RawDisp = permute(DU.Phy.RawDisp, [1 2 3]);
                DU.Phy.VPC = permute(DU.Phy.VPC,[1 2 3]);
                DU.Phy.Ux  = permute(DU.Phy.Ux,[1 2 3]);
                DU.Phy.Uy  = permute(DU.Phy.Uy,[1 2 3]);
                DU.Phy.Uz  = permute(DU.Phy.Uz,[1 2 3]);
                DU.Phy.X   = permute(DU.Phy.X,[1 2 3]);
                DU.Phy.Y   = permute(DU.Phy.Y,[1 2 3]);
                DU.Phy.Z   = permute(DU.Phy.Z,[1 2 3]);
                DU.Seg3D   = permute(DU.Seg3D,[1 2 3]);
                
    case 2;     DU.RawDisp = permute(DU.RawDisp, [2 1 3]);
                DU.VPC     = permute(DU.VPC, [2 1 3]);
                DU.Ux      = permute(DU.Ux, [2 1 3]);
                DU.Uy      = permute(DU.Uy, [2 1 3]);
                DU.Uz      = permute(DU.Uz, [2 1 3]);
                DU.Y1      = permute(DU.X, [2 1 3]);
                DU.X1      = permute(DU.Y, [2 1 3]);
                DU.Z       = permute(DU.Z, [2 1 3]);
            DU.Phy.RawDisp = permute(DU.Phy.RawDisp, [2 1 3]);
                DU.Phy.VPC = permute(DU.Phy.VPC,[2 1 3]);
                DU.Phy.Ux  = permute(DU.Phy.Ux,[2 1 3]);
                DU.Phy.Uy  = permute(DU.Phy.Uy,[2 1 3]);
                DU.Phy.Uz  = permute(DU.Phy.Uz,[2 1 3]);
                DU.Phy.Y1  = permute(DU.Phy.X,[2 1 3]);
                DU.Phy.X1  = permute(DU.Phy.Y,[2 1 3]);
                DU.Phy.Z   = permute(DU.Phy.Z,[2 1 3]);
                DU.Seg3D   = permute(DU.Seg3D,[2 1 3]);
                DU.Y     = DU.Y1;           DU.X = DU.X1;   
                DU.Phy.Y = DU.Phy.Y1;   DU.Phy.X = DU.Phy.X1;
                
    case 3;     DU.RawDisp = permute(DU.RawDisp, [3 2 1]);
                DU.VPC     = permute(DU.VPC, [3 2 1]);
                DU.Ux      = permute(DU.Ux, [3 2 1]);
                DU.Uy      = permute(DU.Uy, [3 2 1]);
                DU.Uz      = permute(DU.Uz, [3 2 1]);
                DU.Z1      = permute(DU.X, [3 2 1]);
                DU.Y       = permute(DU.Y, [3 2 1]);
                DU.X1      = permute(DU.Z, [3 2 1]);
            DU.Phy.RawDisp = permute(DU.Phy.RawDisp, [3 2 1]);
                DU.Phy.VPC = permute(DU.Phy.VPC,[3 2 1]);
                DU.Phy.Ux  = permute(DU.Phy.Ux,[3 2 1]);
                DU.Phy.Uy  = permute(DU.Phy.Uy,[3 2 1]);
                DU.Phy.Uz  = permute(DU.Phy.Uz,[3 2 1]);
                DU.Phy.Z1  = permute(DU.Phy.X,[3 2 1]);
                DU.Phy.Y   = permute(DU.Phy.Y,[3 2 1]);
                DU.Phy.X1  = permute(DU.Phy.Z,[3 2 1]);   
                DU.Seg3D   = permute(DU.Seg3D,[3 2 1]); 
                DU.Z     = DU.Z1;           DU.X = DU.X1;   
                DU.Phy.Z = DU.Phy.Z1;   DU.Phy.X = DU.Phy.X1;
    end
    
DU.Phy.Seg_Uz = DU.Phy.Uz;     DU.Phy.Seg_Uz(DU.Seg3D==1)=NaN;   
DU.Phy.Seg_Ux = DU.Phy.Ux;     DU.Phy.Seg_Ux(DU.Seg3D==1)=NaN;       	
DU.Phy.Seg_Uy = DU.Phy.Uy;     DU.Phy.Seg_Uy(DU.Seg3D==1)=NaN; 
DU.Phy.Seg_X  = DU.Phy.X;      DU.Phy.Seg_X(DU.Seg3D==1)=NaN;   
DU.Phy.Seg_Y  = DU.Phy.Y;      DU.Phy.Seg_Y(DU.Seg3D==1)=NaN;       	
DU.Phy.Seg_Z  = DU.Phy.Z;      DU.Phy.Seg_Z(DU.Seg3D==1)=NaN; 