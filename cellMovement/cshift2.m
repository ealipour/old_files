function [B] = cshift2(A,compass)

[W,L] = size(A);

if compass(1)<0
    
    Iy = [-compass(1)+1:W 1:-compass(1)];
    
elseif compass(1)>0
    
    Iy = [W-compass(1)+1:W 1:W-compass(1)];
    
else
    
    Iy = (1:W);
    
end

if compass(2)<0
    
    Ix = [-compass(2)+1:L 1:-compass(2)];
    
elseif compass(2)>0
    
    Ix = [L-compass(2)+1:L 1:L-compass(2)];
    
else
    
    Ix = (1:L);
    
end

B = A(Iy,Ix);