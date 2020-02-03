function [ux,uy]=point_force

%number of points = 20/resolution=20/0.1=200
GridNum = 101;    


%bulk and shear moduli 
v = 1.4;
%1/(2b)
b = 1/(4.0 *(1.0-v));

mu = 0.21;
%K = 2.0;

x =linspace(-4.,4.,GridNum);
[X,Y] = meshgrid(x,x);

%radius of the cell is considered 1; 

Mask = false(GridNum,GridNum);
Mask(X.^2+Y.^2>0.0)= true;

ux = zeros(GridNum, GridNum);
uy = zeros(GridNum, GridNum);
r =  zeros(GridNum, GridNum); 




Fxq = 1.0;
Fyq = 0.0;
xq=0.0;
yq=0.0;

% theta = linspace(0,2.0*pi*(1.0-1.0/GridNum),GridNum);
% deltal = 2.0 * sin(2.0*pi/GridNum);
% xq = sin(theta);
% yq = cos(theta);


r(Mask) = ((X(Mask)-xq).^2+(Y(Mask)-yq).^2).^0.5;



G11 = 1./r(Mask).*( (1-b)+b.* ((X(Mask)-xq).^2)./(r(Mask).^2));
G12 = 1./r(Mask).* b.* (X(Mask)-xq).*(Y(Mask)-yq)./(r(Mask).^2);
G21 = G12;
G22 = 1./r(Mask).*( (1-b)+b.* ((Y(Mask)-yq).^2)./(r(Mask).^2) );

    
    

ux(Mask) = ux(Mask)+1.0/(4.0*pi*mu) .*(Fxq*G11+ Fyq*G12);
                                                      
uy(Mask) = uy(Mask)+1.0/(4.0*pi*mu) .*(Fyq*G22+ Fxq*G21);


    
quiver(X(Mask),Y(Mask), ux(Mask), uy(Mask));
xlabel('x');
ylabel('y');
figure
View=(ux.^2+uy.^2).^0.5;
View(~Mask)=NaN;
surf(X,Y,View);shading flat
xlabel('x');
ylabel('y');


ux(~Mask)=NaN;
uy(~Mask)=NaN;
