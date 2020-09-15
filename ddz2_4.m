function D=ddz2_4(z)
% compute second derivative using 4th order centred 
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
    D(k,k-2) = -1/12;
    D(k,k-1) = 4/3;
    D(k,k)   = -5/2;
    D(k,k+1) = 4/3;
    D(k,k+2) = -1/12;
end 
 
  D(1,1) = 2; D(1,2) = -5; D(1,3) = 4; D(1,4) = -1;
  D(2,2) = 2; D(2,3) = -5; D(2,4) = 4; D(2,5) = -1;

  D(N,N)   = 2;  D(N,N-1)   = -5;  D(N,N-2)   = 4; D(N,N-3)   = -1;
  D(N-1,N-1) = 2;  D(N-1,N-2) = -5;  D(N-1,N-3) = 4; D(N-1,N-4) = -1;

  D = D/del^2;     

end
