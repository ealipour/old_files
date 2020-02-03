function Shift = makeSHIFT(Tall,Wide)
  
  Area = Tall.* Wide;
        
%  S(:,5) = [ 3 6 9    S(:,4) = [ 3 6 9    S(:,3) = [ 9 3 6
%             1 4 7               1 4 7               7 1 4
%             2 5 8 ]             2 5 8 ]             8 2 5 ]
%     
%  S(:,6) = [ 4 7 1    S(:,1) = [ 1 4 7    S(:,2) = [ 7 1 4
%             5 8 2               2 5 8               8 2 5
%             6 9 3 ]             3 6 9 ]             9 3 6 ]
%   
%  S(:,7) = [ 5 8 2    S(:,8) = [ 2 5 8    S(:,9) = [ 8 2 5
%             6 9 3               3 6 9               9 3 6
%             4 7 1 ]             1 4 7 ]             7 1 4 ] 
    
  Shift      = zeros(Area,9);                        % None
  Shift(:,1) =     1:Area   ;
  
  Shift([(1+Tall):Area 1:Tall],2) = Shift(:,1);      % East
  
  Shift(    [2:Area 1],4) = Shift(             :,1); % North
  Shift(   1:Tall:Area,4) = Shift(Tall:Tall:Area,1);  

  Shift(:,6) = Shift([(1+Tall):Area 1:Tall],1);      % West
  
  Shift(             :,8) = Shift(    [2:Area 1],1); % South
  Shift(Tall:Tall:Area,8) = Shift(   1:Tall:Area,1);
  
  Shift(:,3) = Shift(Shift(:,4),2); % North-East
  Shift(:,5) = Shift(Shift(:,6),4); % North-West
  Shift(:,7) = Shift(Shift(:,8),6); % South-West
  Shift(:,9) = Shift(Shift(:,2),8); % South-East
