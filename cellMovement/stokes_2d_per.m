function [vx,vy]=stokes_2d_per(theta)

N=length(theta);

h=1.0;

x=linspace(-pi,pi-2*pi/100,N);
y=linspace(0,3.0*h,N);

dx=x(2)-x(1);

for m=1:N
    for n=1:N
        for j=1:N
            
            dx1=(x(m)-x(j));
            dy1=(y(n)-h);
            dy2=(y(n)+h);
            denom1=cosh(dy1)-cos(dx1);
            denom2=cosh(dy2)-cos(dx1);
            
            if denom1==0
                denom1=1;
            end
            if denom2==0
                denom2=1;
            end
            
            vx1(n,m,j)=-theta(j)*(dy1*sin(dx1)/denom1/2-dy1*sin(dx1)/denom2/2 ...
                       -y(n)*sinh(dy2)*sin(dx1)/denom2^2);
            
            vy1(n,m,j)=-theta(j)*(0.5*log(denom2/denom1)+dy1*sinh(dy1)/denom1/2-dy2*sinh(dy2)/denom2/2 ...
                        +h*y(n)*(cosh(dy2)/denom2-sinh(dy2)^2/denom2^2));
                    
            %vx2(n,m,j)=-0.1*sin(x(j))*(0.5*log(denom2/denom1)-dy1*sinh(dy1)/denom1/2+dy2*sinh(dy2)/denom2/2 ...
            %            -h*y(n)*(cosh(dy2)/denom2-sinh(dy2)^2/denom2^2));
                    
            %vy2(n,m,j)=-0.1*sin(x(j))*(dy1*sin(dx1)/denom1/2-dy1*sin(dx1)/denom2/2 ...
            %           +y(n)*sinh(dy2)*sin(dx1)/denom2^2);
           
        end
        
        vx(n,m)=dx*sum(vx1(n,m,:));%+vx2(n,m,:));
        vy(n,m)=dx*sum(vy1(n,m,:));%+vy2(n,m,:));
    end
end


pcolor(vx),shading interp