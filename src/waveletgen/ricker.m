% Ricker is the 2nd derivative of a gaussian
%  f:   output Ricker wavelet
%  n :  number of data points returned
%  fpeak: peak frequency (khz). --> sqrt(2)/tau for Ricker
%  dt:  time sampling rate

function f = ricker(fpeak, dt)

fc = fpeak;
timeshift = sqrt(2);
nts = floor(timeshift/fc/dt+0.5);

% the total length of f[i] should be 2*nts+1
b = (3.1415926*fc*1000)*(3.1415926*fc*1000);
a = 2.*b;

for i=1:nts+1
    t = (i-1)*dt*0.001;
    t2 = t*t;
    f(nts+i) = (1.-a*t2)*exp(-b*t2); % nts.... 2*nts */
    f(nts+2-i) = f(nts+i);                % nts ...0-->nts....2*nts */
end

