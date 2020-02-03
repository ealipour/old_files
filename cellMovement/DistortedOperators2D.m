function [LapXX,LapYY,LapXY,LapYX,GradX,GradY,GradXEdge,GradYEdge,VolMat,ControlVolume,EdgeLength,Link] ...
         = DistortedOperators2D(Cell,Edge,Distance,Dxx,Dxy,Dyx,Dyy,Nx,Ny,Dx,Dy)

%% Function outputs matrices that can be used
%% for finite volume algorithms.  The LapXX and LapYY are the x and y components of a
%% generalized Laplacian matrix with diffusion coefficients Dxx and Dyy.
%%  LapXY and LapYX are the XY second derivatives with diffusion (or shear) coefficients
%% Dxy and Dyx. These matrices compute the sum of fluxex
%% (assumed to be proportional to the
%% gradient of the function) through the faces that enclose the
%% ControlVolumes about the nodes.  GradX and GradY calculate
%% derivative matrices along the x and y directions for use when
%% calculating the divergence of a vector field (i.e., it is not the
%% derivatives of a function but the components of the integral of a
%% divergence).  The Edge components of the components of the gradients are not
%% included, but are computed with the matrices GradXEdge and GradYEdge.

%% define compass

  compass = [  1  0  1  0;
               0  1  0  1 ];
           
  DelX = [ 0   Dx   0 -Dx;
           Dx  Dx -Dx -Dx;
           Dx  0  -Dx  0 ];
       
  DelY = [ Dy  0  -Dy 0  ;
           Dy -Dy -Dy Dy ;
           0  -Dy  0  Dy ];

[Nl,Nw] = size(Cell);

%% define Link Matrix

    Unknowns = nnz(Cell);
    Link = zeros(Nl,Nw);
    Link(Cell) = 1:Unknowns;
    
%% Create Matrices
 
 LapXX = spalloc(Unknowns,Unknowns,0);
 LapXY = LapXX;
 LapYX = LapXX;
 LapYY = LapXX;
 GradX = LapXX;
 GradY = LapXX;
 GradXEdge = LapXX;
 GradYEdge = LapXX;
 VolMat = LapXX;
 ControlVolume = zeros(Unknowns,1);
 
 Fx = zeros(Nl,Nw);
 Fy = zeros(Nl,Nw);
 
 Fx(Edge) = -Distance(Edge).*Nx(Edge);
 Fy(Edge) = -Distance(Edge).*Ny(Edge);

     M1 = Cell                       ...
          & cshift2( Cell, [-1 0])  ...
          & cshift2( Cell, [0 -1]) ...
          & cshift2( Cell, [-1 -1]);
      
     A = cshift2(M1,[1 0]);
     B = cshift2(M1,[1 1]);
     C = cshift2(M1,[0 1]);
     
 %% Construct the connectivity matrices
     
Shift = makeSHIFT(Nl,Nw);

Pull = [ 1 8 7 6;
         1 6 5 4;
         1 4 3 2;
         1 2 9 8 ];
     
    Order2 = { M1, A, B, C };
     
 %% Define the Laplacian for each quarter face
        
 for j = 1:4
     
     DeltaX = DelX(1,j) + Fx(Shift(Order2{j},Pull(j,2))) - Fx(Shift(Order2{j},Pull(j,1)));
     DeltaY = DelY(1,j) + Fy(Shift(Order2{j},Pull(j,2))) - Fy(Shift(Order2{j},Pull(j,1)));
     
     DeltaX1 = DelX(3,j) + Fx(Shift(Order2{j},Pull(j,4))) - Fx(Shift(Order2{j},Pull(j,1)));
     DeltaY1 = DelY(3,j) + Fy(Shift(Order2{j},Pull(j,4))) - Fy(Shift(Order2{j},Pull(j,1)));
     
     DeltaX2 = 0.25.*( 2.*DelX(2,j) + Fx(Shift(Order2{j},Pull(j,4)))+ Fx(Shift(Order2{j},Pull(j,3))) ...
                                    + Fx(Shift(Order2{j},Pull(j,2))) - 3.*Fx(Shift(Order2{j},Pull(j,1))) );
     DeltaY2 = 0.25.*( 2.*DelY(2,j) + Fy(Shift(Order2{j},Pull(j,4)))+ Fy(Shift(Order2{j},Pull(j,3))) ...
                                    + Fy(Shift(Order2{j},Pull(j,2))) - 3.*Fy(Shift(Order2{j},Pull(j,1))) );
     
%% Define Masks

    Tail = Link(Shift(Order2{j},Pull(j,1)));
    Head = Link(Shift(Order2{j},Pull(j,2)));
    Cross = Link(Shift(Order2{j},Pull(j,3)));
    Over = Link(Shift(Order2{j},Pull(j,4)));
    
%% Find the normals to the faces

     npz = (  DeltaX1.*DeltaY2 - DeltaY1.*DeltaX2 );
        
     mag = sqrt(npz.^2);
     
     VolMat =   VolMat                                         ...
              + sparse(Tail,Tail,7.*mag./48,Unknowns,Unknowns) ...
              + sparse(Tail,Over,3.*mag./48,Unknowns,Unknowns) ...
              + sparse(Tail,Head,mag./48,Unknowns,Unknowns)    ...
              + sparse(Tail,Cross,mag./48,Unknowns,Unknowns);    
     
     ControlVolume(Tail) = ControlVolume(Tail) + mag./4;
                                 
     npz = (  DeltaX.*DeltaY2 - DeltaY.*DeltaX2 );
        
     mag = sqrt(npz.^2);
     
     VolMat =   VolMat                                         ...
              + sparse(Tail,Tail,7.*mag./48,Unknowns,Unknowns) ...
              + sparse(Tail,Head,3.*mag./48,Unknowns,Unknowns) ...
              + sparse(Tail,Over,mag./48,Unknowns,Unknowns)    ...
              + sparse(Tail,Cross,mag./48,Unknowns,Unknowns);
     
     ControlVolume(Tail) = ControlVolume(Tail) + mag./4;
 
 %% Normal to flux surface times the length
     
     nx = (DeltaY2 - DeltaY./2);
     ny = -(DeltaX2 - DeltaX./2);
     
     Len = sqrt(nx.^2 + ny.^2);
     
     nx = nx./Len;
     ny = ny./Len;

%% Calculate the Laplacian                                          
    
    Denom =   DeltaX2.*DeltaY - DeltaX.*DeltaY2;
    
%% compute flux terms propotional to the x component of the gradient

    Val = zeros(nnz(Tail),1);
    
    Val = Len.*(nx.*DeltaY2)./Denom;
    
    Dif = (3.*Dxx(Shift(Order2{j},Pull(j,1))) + 3.*Dxx(Shift(Order2{j},Pull(j,2))) ...
            + Dxx(Shift(Order2{j},Pull(j,3))) + Dxx(Shift(Order2{j},Pull(j,4))))./8;
    
    Val2 = zeros(nnz(Tail),1);
    
    Val2 = -Len.*(nx.*DeltaY)./Denom./4;
   
    LapXX =   LapXX ...
            - sparse(Tail,Tail,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Tail,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Over,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Over,Val2.*Dif,Unknowns,Unknowns);    
 
    Val = Len.*(ny.*DeltaY2)./Denom;
    Val2 = -Len.*(ny.*DeltaY)./Denom./4;
            
    Dif = (3.*Dyx(Shift(Order2{j},Pull(j,1))) + 3.*Dyx(Shift(Order2{j},Pull(j,2))) ...
            + Dyx(Shift(Order2{j},Pull(j,3))) + Dyx(Shift(Order2{j},Pull(j,4))))./8;
    
   
    LapYX =   LapYX ...
            - sparse(Tail,Tail,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Tail,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Over,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Over,Val2.*Dif,Unknowns,Unknowns);
        
%%  compute flux terms propotional to the y component of the gradient  
        
    Val = zeros(nnz(Tail),1);
    
    Val = -Len.*(ny.*DeltaX2)./Denom;
    
    Dif = (3.*Dyy(Shift(Order2{j},Pull(j,1))) + 3.*Dyy(Shift(Order2{j},Pull(j,2))) ...
            + Dyy(Shift(Order2{j},Pull(j,3))) + Dyy(Shift(Order2{j},Pull(j,4))))./8;
    
    
    Val2 = zeros(nnz(Tail),1);
    
    Val2 = Len.*(ny.*DeltaX)./Denom./4;
   
    LapYY =   LapYY ...
            - sparse(Tail,Tail,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Tail,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Over,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Over,Val2.*Dif,Unknowns,Unknowns);
        
    Val = -Len.*(nx.*DeltaX2)./Denom;
    
    Dif = (3.*Dxy(Shift(Order2{j},Pull(j,1))) + 3.*Dxy(Shift(Order2{j},Pull(j,2))) ...
            + Dxy(Shift(Order2{j},Pull(j,3))) + Dxy(Shift(Order2{j},Pull(j,4))))./8;
    
    
    Val2 = Len.*(nx.*DeltaX)./Denom./4;
   
    LapXY =   LapXY ...
            - sparse(Tail,Tail,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val.*Dif,Unknowns,Unknowns) ...
            - sparse(Tail,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Head,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Tail,Over,Val2.*Dif,Unknowns,Unknowns) ...
            + sparse(Head,Tail,3.*Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Head,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Cross,Val2.*Dif,Unknowns,Unknowns) ...
            - sparse(Head,Over,Val2.*Dif,Unknowns,Unknowns);
   
%% calculate the 1st derivative matrices
            
   Val = Len./8;

   GradX =   GradX ...
           - sparse(Tail,Tail,3.*Val.*nx,Unknowns,Unknowns)  ...
           - sparse(Tail,Head,3.*Val.*nx,Unknowns,Unknowns)  ...
           - sparse(Tail,Cross,Val.*nx,Unknowns,Unknowns)    ...
           - sparse(Tail,Over,Val.*nx,Unknowns,Unknowns)     ...
           + sparse(Head,Tail,3.*Val.*nx,Unknowns,Unknowns)  ...
           + sparse(Head,Head,3.*Val.*nx,Unknowns,Unknowns)  ...
           + sparse(Head,Cross,Val.*nx,Unknowns,Unknowns)    ...
           + sparse(Head,Over,Val.*nx,Unknowns,Unknowns);
       
  GradY =    GradY ...
           - sparse(Tail,Tail,3.*Val.*ny,Unknowns,Unknowns)  ...
           - sparse(Tail,Head,3.*Val.*ny,Unknowns,Unknowns)  ...
           - sparse(Tail,Cross,Val.*ny,Unknowns,Unknowns)    ...
           - sparse(Tail,Over,Val.*ny,Unknowns,Unknowns)     ...
           + sparse(Head,Tail,3.*Val.*ny,Unknowns,Unknowns)  ...
           + sparse(Head,Head,3.*Val.*ny,Unknowns,Unknowns)  ...
           + sparse(Head,Cross,Val.*ny,Unknowns,Unknowns)    ...
           + sparse(Head,Over,Val.*ny,Unknowns,Unknowns);

 end
 
 %% calculate edge lengths along outside edge

 EdgeLength = zeros(nnz(Cell),1);
 
 DelX = [ 0  Dx ];
 DelY = [ Dy 0 ];
 
 Pull = [ 1 4;
          1 2 ];
      
      compass2 = [  0  1 0 1 1 1 -1 -1; 
                   -1 -1 1 1 0 1  0  1 ];
 
      Core = Cell & ~Edge;
      
 for i = 1:2
     
     Head = Edge & cshift2(Edge,compass(:,i)) & ( cshift2(Core,compass2(:,4.*(i-1)+1)) | cshift2(Core,compass2(:,4.*(i-1)+2)) );
   
     Tail = Shift(Head,Pull(i,2));
     
     Len = sqrt(   (DelX(i) + Fx(Head)- Fx(Tail)).^2 ...
                 + (DelY(i) + Fy(Head)- Fy(Tail)).^2 );
             
     EdgeLength(Link(Head)) = EdgeLength(Link(Head)) + Len./2;
     EdgeLength(Link(Tail)) = EdgeLength(Link(Tail)) + Len./2;
     
     Val = -(DelY(i) + Fy(Head)- Fy(Tail))./8;
     GradXEdge =   GradXEdge ...
                 + sparse(Link(Head),Link(Head),3.*Val,Unknowns,Unknowns) ...
                 + sparse(Link(Head),Link(Tail),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Head),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Tail),3.*Val,Unknowns,Unknowns);
     
                 
     Val = (DelX(i) + Fx(Head)- Fx(Tail))./8;
     GradYEdge = GradYEdge ...
                 + sparse(Link(Head),Link(Head),3.*Val,Unknowns,Unknowns) ...
                 + sparse(Link(Head),Link(Tail),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Head),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Tail),3.*Val,Unknowns,Unknowns);

      Head = Edge & cshift2(Edge,compass(:,i)) & ( cshift2(Core,compass2(:,4.*(i-1)+3)) | cshift2(Core,compass2(:,4.*(i-1)+4)) );
   
     Tail = Shift(Head,Pull(i,2));
     
     Len = sqrt(   (DelX(i) + Fx(Head)- Fx(Tail)).^2 ...
                 + (DelY(i) + Fy(Head)- Fy(Tail)).^2 );
             
     EdgeLength(Link(Head)) = EdgeLength(Link(Head)) + Len./2;
     EdgeLength(Link(Tail)) = EdgeLength(Link(Tail)) + Len./2;
     
     Val = (DelY(i) + Fy(Head)- Fy(Tail))./8;
     GradXEdge =   GradXEdge ...
                 + sparse(Link(Head),Link(Head),3.*Val,Unknowns,Unknowns) ...
                 + sparse(Link(Head),Link(Tail),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Head),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Tail),3.*Val,Unknowns,Unknowns);
     
                 
     Val = -(DelX(i) + Fx(Head)- Fx(Tail))./8;
     GradYEdge = GradYEdge ...
                 + sparse(Link(Head),Link(Head),3.*Val,Unknowns,Unknowns) ...
                 + sparse(Link(Head),Link(Tail),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Head),Val,Unknowns,Unknowns) ...
                 + sparse(Link(Tail),Link(Tail),3.*Val,Unknowns,Unknowns);
                        
 end
 
 EdgeLength = EdgeLength(Link(Edge));
 