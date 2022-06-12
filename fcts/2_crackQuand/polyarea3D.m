function [area] = polyarea3D(x,y,z)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               compute the area of a 3D planar polygon
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Triangle decomposition
%
%       [area] = polyarea3D(x,y,z)
%       ===========================
%   
%   Inputs:
%   -------
%       > x,y,z: coordinates of the vertices of the polygon
%
%   Output:
%   -------
%       > area: the area of the polygon
%
%   written by Yang CHEN 14/06/2016
%

%% ---- compute the mass center of the planar polygon
G = [mean(x(:)) mean(y(:)) mean(z(:))];

%% ---- sort the vetices into counterclockwise order
% this is done over the coordinates of the polygon's plan
% 
% construct the local axes of the polygon's plan
GPi = [x(:)-G(1), y(:)-G(2), z(:)-G(3)];
ex = GPi(1,:)./norm(GPi(1,:));
ez = cross(GPi(1,:),GPi(2,:));    ez = ez./norm(ez);
ey = cross(ez,ex);

    % a=GPi_dot_ex;  b=GPi_dot_ey
a = GPi(:,1).*ex(1) + GPi(:,2).*ex(2) + GPi(:,3).*ex(3);
b = GPi(:,1).*ey(1) + GPi(:,2).*ey(2) + GPi(:,3).*ey(3);

    % Q1:   a>=0; b>0       => angle = asin(a/rho)
    % Q2:   a<0;  b>=0      => angle = asin(a/rho)
    % Q3:   a<+0; b<0       => angle = 2pi - asin(a/rho)
    % Q4:   a>0;  b<=0      => angle = 2pi - asin(a/rho)
    %                          with rho=amplitude_of_GPi
rho = sqrt(sum(GPi.*GPi,2));

i0 = b>=0;
i1 = b<0;

angle = zeros(length(x),1);
angle(i0) = acos(a(i0)./rho(i0));
angle(i1) = 2*pi - acos(a(i1)./rho(i1));

% sort the angles of each vertex
[angle,isort] = sort(angle,'ascend');

% sort the vertices
x = x(isort);
y = y(isort);
z = z(isort);

%% ---- compute the area of the polygon via trianglular decomposition
GPi = [x(:)-G(1), y(:)-G(2), z(:)-G(3)];
    % area = |a|*|b|* sin(theta)/2
    %      = |axb|/2
a_x = GPi(:,1);
a_y = GPi(:,2);
a_z = GPi(:,3);

GPii = [GPi(2:end,:);GPi(1,:)];
b_x = GPii(:,1);
b_y = GPii(:,2);
b_z = GPii(:,3);

crossProd = [a_y.*b_z - b_y.*a_z,...
             b_x.*a_z - a_x.*b_z,...
             a_x.*b_y - b_x.*a_y];

area = sum(sqrt(sum(crossProd.^2,2)))/2;

%% =======================================================================
% %% %%%% test
% x=[0 1 0 2 2 1]';
% y=[0 0.75 0.5 0 1 1]';
% z=[0 0 0 0 0 0]';
%     figure;plot3(x,y,z,'o');axis equal
%     xlabel('x');ylabel('y');zlabel('z');
% P = [x, y, z];
% % rotation 3D
% phi=pi/6;
% betaZ = pi/10;
% psi=0;
% [Ry,Rz1,Ry2] = RotMatrix(phi,betaZ,psi);
% P = P*Ry;
% P = P*Rz1;
% P = P*Ry2;
% 
% x = P(:,1);
% y = P(:,2);
% z = P(:,3);
%     hold on;plot3(x,y,z,'or')
% [area] = polyarea3D(x,y,z)
