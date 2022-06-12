function [cONa,bONa] = AspRatio_cylin(VxLstI,sizVOL,Ot,ex,ey,ez)


% from index list to coordinate list
[y,x,z] = ind2sub(sizVOL,VxLstI);

% from Cartesien coordinates to Cylindrical coordinates
[r,t,z] = XYZ2RTZ(Ot,ex,ey,ez,x,y,z);

% compute the mass center in the cylindrical coordinates
r0 = mean(r);
t0 = mean(t);
z0 = mean(z);

% compute the inertia tensor to the mass center
r = r-r0;
t = r0.*(t-t0);
z = z-z0;

Irr = sum(t.^2 + z.^2);
Itt = sum(r.^2 + z.^2);
Izz = sum(r.^2 + t.^2);
Itz = -sum(t.*z);
Irz = -sum(r.*z);
Irt = -sum(r.*t);

I = [Irr, Irt, Irz;
     Irt, Itt, Itz;
     Irz, Itz, Izz];
 
% compute the eigenvectors and eigenvalues
[dir,amp]=eig(I);
amp = [amp(1,1),amp(2,2),amp(3,3)];
 
cONa = sqrt((amp(1)+amp(2)-amp(3))/(amp(2)+amp(3)-amp(1)));
bONa = sqrt((amp(1)+amp(3)-amp(2))/(amp(2)+amp(3)-amp(1)));

