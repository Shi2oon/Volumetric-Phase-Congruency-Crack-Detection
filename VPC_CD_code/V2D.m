function [X, Y, Z, Ux, Uy, Uz] = V2D(D)
% j=-1;
% k=1;
%   
% while j<0
%     if D(k,2)==D(k+1,2)
%         k=k+1;
%     else
%         lx=k;
%         j=1;
%     end
% end
% ly=size(D,1)/lx;
%     X=reshape(D(:,1),lx,ly);
%     Y=reshape(D(:,2),lx,ly);
%     Z=reshape(D(:,3),lx,ly);
%     
%     Ux=reshape(D(:,4),lx,ly);
%     Uy=reshape(D(:,5),lx,ly);
%     Uz=reshape(D(:,6),lx,ly);
% end

j=-1;
i=1;
while j<0
   if D(i,2)==D(i+1,2)
       i=i+1;
   else
       lx=i;
       j=1;
   end
end
lx;

tempZ=reshape(D(:,3),lx,size(D,1)/lx);
%%

j=-1;
i=1;
while j<0
   if tempZ(1,i)==tempZ(1,i+1)
       i=i+1;
   else
       ly=i;
       j=1;
   end
end
ly;
lz=size(D,1)/ly/lx;

    X=reshape(D(:,1),lx,ly,lz);
    Y=reshape(D(:,2),lx,ly,lz);
    Z=reshape(D(:,3),lx,ly,lz);
    
    Ux=reshape(D(:,4),lx,ly,lz);
    Uy=reshape(D(:,5),lx,ly,lz);
    Uz=reshape(D(:,6),lx,ly,lz);