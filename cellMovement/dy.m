function[duxy,duyy] = dy(sigma)

GridNum = 140;    

%constants to stress field
epsilon = -0.3;
B = -2.55; 
A = B/4.0;
radius = 1.0;

v0 = 0.003;
zeta = 300.0;
%bulk and shear moduli 
% mu = E/(4*(1+sigma));
%K= E/(3*(1-2*sigma));

%b= 2(1-sigma); 1/(2b)
b = 1/(4.0 *(1.0-sigma));

x =linspace(-4.0*radius,4.0*radius,GridNum);
[X,Y] = meshgrid(x,x);
delta_x= 2*4.0*radius/GridNum;
t = linspace(0, 2*pi*(1-1/GridNum), GridNum);


Cell = false(GridNum,GridNum);
Cell((X.^2+Y.^2) > radius)= true;
duxy = zeros(GridNum,GridNum);
duyy = zeros(GridNum,GridNum);

for m=1:GridNum
    tm = t(m);
    xc = zeros(GridNum,GridNum);
    xc(Cell) = X(Cell) - sin(tm);
    yc = zeros(GridNum,GridNum);
    yc(Cell) = Y(Cell)-cos(tm);
    r = zeros(GridNum,GridNum);
    r = (xc.^2 + yc.^2).^0.5;
    
    
    
    duxy(Cell) = duxy(Cell) + zeta* v0 * (1+epsilon+epsilon*cos(tm))./ (r(Cell).^5) ...
        .* (( (1+2*b)* xc(Cell).^2 +(1-b)*yc(Cell).^2)...
        * A*sin(tm)....
        + b*xc(Cell).*(xc(Cell).^2-2*yc(Cell).^2).*(1+B/2+B/2*cos(tm)));
    
    
    duyy(Cell) = duyy(Cell)-(1+epsilon+epsilon*cos(t(m)))...
                ./(xc(Cell).^2 + yc(Cell).^2).^2.5...
        .*((A*sin(t(m))*b*xc(Cell).*(xc(Cell).^2-2*yc(Cell).^2)...
        +(1+B/2+B/2*cos(t(m))).*yc(Cell).*((1-3*b)*xc(Cell).^2+yc(Cell).^2)));
        
    
end
SubMask = false(GridNum,GridNum);
SubMask(1:10:GridNum,1:10:GridNum)= true;

mymat = [X(Cell);Y(Cell);duxy(Cell);duyy(Cell)];


fid = fopen('duy.txt','a+');
fprintf(fid,'%2.6f\t %2.6f\t %2.6f\t %2.6f\n', mymat');
fclose(fid);


% figure    
%  quiver(X(SubMask&Cell),Y(SubMask&Cell),duxy(SubMask& Cell),duyy(SubMask&Cell),1,'k','Linewidth',1);
%Surface(X,Y, duxy);
% xlabel('x');
% ylabel('y');