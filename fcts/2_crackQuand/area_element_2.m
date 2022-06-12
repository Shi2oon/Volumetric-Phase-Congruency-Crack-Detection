function [dS] = area_element_2(N,G,xV,yV,zV)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%  Compute the intersection surface area between a plan and a cuboid
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	
%	[dS] = area_element_2(N,G,xV,yV,zV)
%	-----------------------------------
%
%    Inputs:
%    -------
%	>          N: the normal of the plan (matrix nx3)
%	>          G: the point that the plan is passing through (matrix nx3)
%	> [xV yV zV]: the 8 vertices of the cuboid (matrix nx8)
%
%    Output:
%    -------
%	> dS: the intersection area(s) (matrix nx1)
%
%	Note: n denotes the number of couples of cuboids and plans
%
%    written by Yang CHEN 02/06/2016
%

nbox = size(N,1);

%% find the intersection vertices
xI=zeros(nbox,12);    yI=xI;    zI=xI;
for i=1:12
    [n1,n2] = num2vertices(i); %each side line is associated to two vertices
    V1 = [xV(:,n1), yV(:,n1), zV(:,n1)];
    V2 = [xV(:,n2), yV(:,n2), zV(:,n2)];
    a = V2(:,1) - V1(:,1);
    b = V2(:,2) - V1(:,2);
    c = V2(:,3) - V1(:,3);
    tmp1 = ( N(:,1).*(G(:,1)-V1(:,1)) + ...
             N(:,2).*(G(:,2)-V1(:,2)) + ...
             N(:,3).*(G(:,3)-V1(:,3)) );
    tmp2 = ( a.*N(:,1) + b.*N(:,2) + c.*N(:,3) );
    
    xI(:,i) = V1(:,1) + a .* tmp1 ./ tmp2;
    yI(:,i) = V1(:,2) + b .* tmp1 ./ tmp2;
    zI(:,i) = V1(:,3) + c .* tmp1 ./ tmp2;
    
        % note: in the case that the S passes through the L (tmp2==0 <->
        % S//L;  tmp1==1 <-> V1 is in S), there will be no intersection
        % points for this being-passed L, however, the two vertices of the
        % being-passed L will be double-counted by the L's four
        % intersecting lines. Since the function "polyarea3D.m" is able to
        % provide the accurate area value even the polygon's vetices are
        % repeately presented, so this special case (S passing through L)
        % can be also correctly considered in this code. 
    
    % intersection points should be contained within the cuboid sides
    i0 = ( xI(:,i)-V1(:,1)>=0 & xI(:,i)-V2(:,1)<=0 ) & ...
         ( yI(:,i)-V1(:,2)>=0 & yI(:,i)-V2(:,2)<=0 ) & ...
         ( zI(:,i)-V1(:,3)>=0 & zI(:,i)-V2(:,3)<=0 );
    i0 = ~i0;
    xI(i0,i) = NaN;
    yI(i0,i) = NaN;
    zI(i0,i) = NaN;
    i0 = xI==Inf;   %plan//cuboid side&cuboid side does not belong to plan
    xI(i0,i) = NaN;
    yI(i0,i) = NaN;
    zI(i0,i) = NaN;
end


%% compute the intersection surface area from the intersection points
dS=zeros(nbox,1);
for i=1:nbox
    i0 = ~isnan(xI(i,:));
    x = xI(i,i0);
    y = yI(i,i0);
    z = zI(i,i0);
    if isempty(x); dS(i)=0; continue; end
    dS(i) = polyarea3D(x,y,z);
end



%%     Each side line of the cuboid is associated to two vertices
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n1,n2] = num2vertices(n)
switch n
    case 1
        n1 = 1;    n2 = 2;
    case 2
        n1 = 2;    n2 = 3;
    case 3
        n1 = 4;    n2 = 3;
    case 4
        n1 = 1;    n2 = 4;
    case 5
        n1 = 5;    n2 = 6;
    case 6
        n1 = 6;    n2 = 7;
    case 7
        n1 = 8;    n2 = 7;
    case 8
        n1 = 5;    n2 = 8;
    case 9
        n1 = 1;    n2 = 5;
    case 10
        n1 = 2;    n2 = 6;
    case 11
        n1 = 3;    n2 = 7;
    case 12
        n1 = 4;    n2 = 8;
end




