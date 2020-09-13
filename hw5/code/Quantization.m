function q = Quantization(r,Q)
% Q = 2
if Q == 1
    r(r>=0) = 1;
    r(r<0) = -1;
    q = r;
    
% Q = 4
elseif Q  == 2
    r(r>=0.75) = 1;
    r(r<0.75 & r>=0) = 0.5 ;
    r(r<0 & r>=-0.75 ) = -0.5;
    r(r<-0.75) = -1;
    q = r;
    
% Q = 8
elseif Q == 3
    r = r;
    r(r>=0.875) = 1;
    r(r<0.875 & r>=0.625) = 0.75;
    r(r<0.625 & r>=0.375) = 0.5;
    r(r<0.375 & r>=0) = 0.25;
    r(r<0 & r>=-0.375) = -0.25;
    r(r<-0.375 & r>=-0.625) = -0.5;
    r(r<-0.625 & r>=-0.875) = -0.75;
    r(r<-0.875) = -1;
    q = r;
    
else
    fprintf('Given the wrong quantization level.')
end
end