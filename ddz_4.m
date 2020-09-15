function D=ddz_4(z)
% compute first derivative using 4th order centred 
% difference scheme.
% z is evenly spaced.

% check for equal spacing
if abs(std(diff(z))/mean(diff(z)))>.000001
    disp(['ddz: values not evenly spaced!'])
    d = NaN;
    return
end

del = z(2)-z(1);
N = length(z);

D = zeros(N,N);

for k=3:N-2
    D(k,k-2) = 1/12;
    D(k,k-1) = -2/3;
    D(k,k)   = 0;
    D(k,k+1) = 2/3;
    D(k,k+2) = -1/12;
end 
 
   D(1,1) = -3/2; D(1,2) = 2; D(1,3) = -1/2;
   D(2,2) = -3/2; D(2,3) = 2; D(2,4) = -1/2;
   
   D(N,N)     = 3/2; D(N,N-1)   = -2;  D(N,N-2)   = 1/2;
   D(N-1,N-1) = 3/2; D(N-1,N-2) = -2;  D(N-1,N-3) = 1/2;

   D = D/del;
end


