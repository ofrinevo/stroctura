function [Ig, PrinsIg, PrinsAx,V_cg] = calcPrincipalAxes(X,Y,Z,M)


mt = sum(M);           % The total mass of this system 

% Three components of the position vector of the center of mass are
Xcg =(1/mt)*sum(M.*X); 
Ycg =(1/mt)*sum(M.*Y);
Zcg =(1/mt)*sum(M.*Z); 
V_cg = [Xcg,Ycg,Zcg]   % [m]     
%%
Ig = [0.0,0.0,0.0;
      0.0,0.0,0.0;
      0.0,0.0,0.0];
%%  
% The total moment of inertia is the sum of moments of inertia for all 
% point masses in the system
for i =1:length(M)
    x = (X(i) - Xcg);
    y = (Y(i) - Ycg);
    z = (Z(i) - Zcg);
    m = M(i);
 Ig = Ig + [  m*(y^2 + z^2),    -m*x*y,         -m*x*z;
             -m*y*x,             m*(x^2 + z^2), -m*y*z;
             -m*x*z,            -m*y*z,          m*(y^2 + x^2)]; % [kg*m^2]
end


[PrinsAx,PrinsIg] = eig(Ig);

%%
