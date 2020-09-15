function D=ddz4_4(z)
% compute fourth derivative using 4th order centred 
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

for k=4:N-3
    D(k,k-3) = -1/6;
    D(k,k-2) = 2;
    D(k,k-1) = -13/2;
    D(k,k)   = 28/3;
    D(k,k+1) = -13/2; 
    D(k,k+2) = 2;
    D(k,k+3) = -1/6;
end 
  
  D(1,1) = 3; D(1,2) = -14;  D(1,3) = 26;  D(1,4) = -24;  D(1,5) = 11;  D(1,6) = -2;
  D(2,2) = 3; D(2,3) = -14;  D(2,4) = 26;  D(2,5) = -24;  D(2,6) = 11;  D(2,7) = -2;
  D(3,3) = 3; D(3,4) = -14;  D(3,5) = 26;  D(3,6) = -24;  D(3,7) = 11;  D(3,8) = -2;
  
  D(N,N)     = 3; D(N,N-1)   = -14; D(N,N-2)   = 26; D(N,N-3)   = -24; D(N,N-4)   = 11; D(N,N-5)   = -2;
  D(N-1,N-1) = 3; D(N-1,N-2) = -14; D(N-1,N-3) = 26; D(N-1,N-4) = -24; D(N-1,N-5) = 11; D(N-1,N-6) = -2;
  D(N-2,N-2) = 3; D(N-2,N-3) = -14; D(N-2,N-4) = 26; D(N-2,N-5) = -24; D(N-2,N-6) = 11; D(N-2,N-7) = -2;

  D = D/del^4;
end
