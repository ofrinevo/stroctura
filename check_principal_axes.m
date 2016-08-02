%
clc 
clear all
close all

% 
n =8 % Circular Symmetry Degree
N = 20 % number of points in monomer

x = [rand; rand; rand]; % rotation axes

dx = [rand; rand; rand]*10; % shift

x = x./norm(x)*10

% points coordinates and masses
X = rand(N,1)*10;
Y = rand(N,1)*10;
Z = rand(N,1)*10;
M_1 = rand(N,1)*10;
M = M_1
% rotationn vector
rot_angle = 2*pi/n;
r = [x;rot_angle];

RotMatrix = vrrotvec2mat(r);

XYZ_rot = [X,Y,Z];
for in_rot = 2:n
    for in = 1:N
        XYZ_rot(in,:) = RotMatrix*XYZ_rot(in,:)';
    end
    X=[X,;XYZ_rot(:,1)];
    Y=[Y;XYZ_rot(:,2)];
    Z=[Z;XYZ_rot(:,3)];
    M=[M;M_1];
end

X =X+dx(1);
Y =Y+dx(2);
Z =Z+dx(3);

plot3(X,Y,Z,'o')
hold on
plot3([dx(1),dx(1)+x(1)],[dx(2),dx(2)+x(2)],[dx(3),dx(3)+x(3)],'r')



[Ig, PrinsIg, PrinsAx,V_cg] = calcPrincipalAxes(X,Y,Z,M)

%plot Principal Axes
plot3([0,PrinsAx(1,1)]+V_cg(1),[0,PrinsAx(2,1)]+V_cg(2),[0,PrinsAx(3,1)]+V_cg(3),'ko:')
plot3([0,PrinsAx(1,2)]+V_cg(1),[0,PrinsAx(2,2)]+V_cg(2),[0,PrinsAx(3,2)]+V_cg(3),'ko:')
plot3([0,PrinsAx(1,3)]+V_cg(1),[0,PrinsAx(2,3)]+V_cg(2),[0,PrinsAx(3,3)]+V_cg(3),'ko:')
axis equal
view(x)
