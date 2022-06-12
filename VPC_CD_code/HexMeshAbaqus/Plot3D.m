function Plot3D(M,X,Y,Z,units,tiz)
if isempty(X)
    [X,Y,Z] = meshgrid(1:size(M,2),1:size(M,1),1:size(M,3)); 
    units = 'Pixel';
    hs = slice(X,Y,Z,M,unique(X),unique(Y),unique(Z)) ;
else
    hs = slice(Y,X,Z,M,unique(X),unique(Y),unique(Z)) ;
end
shading interp;                             set(hs,'FaceAlpha',1);
set(gcf,'position',[30 50 1300 950]);       
y = unique(Y);      xlim([min(y) max(y)]);
x = unique(X);      ylim([min(x) max(x)]);
xlabel(['X [' units ']']);                
ylabel(['Y [' units ']']); 
zlabel(['Z [' units ']']);
axis image; 
view([45 45 15])
c = colorbar;       colormap jet;       c.Label.String = units;
title(tiz)
