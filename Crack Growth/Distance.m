X = x00000;
Y = y00000;
Z = z00000;
Px= x00500;
Py= y00500;
Pz= z00500;
for i=1:23;
    for j=1:23;
        d(i,j)=sqrt((X(i)-Px(j)).^2+(Y(i)-Py(j)).^2+(Z(i)-Pz(j)).^2);     
    end
end
Dmin=min(d);
Dmax=max(d);
Shortest=min(Dmin);
Longest=max(Dmax);
disp(Shortest)
disp(Longest)
 