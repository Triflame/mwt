clear;clc;
path(path,'/uufs/geophys.utah.edu/common/tomofs/chaiwoot/codes/lib/SegyMAT');

[w10,h1,h2] = ReadSu('src10.su');
[w5,h1,h2] = ReadSu('src5.su');
[w5_hamming,h1,h2] = ReadSu('src5_hamming.su');
[w5_wiener,h1,h2] = ReadSu('src5_wiener.su');

nt = length(w10);
t = (0:nt-1)*0.001;
% set(gcf,'Position',[600, 100, 400, 700]);
subplot(221);plot(t,w10);xlabel('Time (seconds)');ylabel('Amplitude');
title('(a) 10-Hz Ricker Wavelet');axis([0 1.1 -1 1]);
subplot(222);plot(t,w5);xlabel('Time (seconds)');ylabel('Amplitude');
title('(b) Targeted 5-Hz Ricker Wavelet');axis([0 1.1 -1 1]);
subplot(223);plot(t,w5_hamming/max(w5_hamming));
xlabel('Time (seconds)');ylabel('Amplitude');
title('(c) 5-Hz Ricker Wavelet using Hamming-Window Filter');axis([0 1.1 -1 1]);
subplot(224);plot(t,w5_wiener);
title('(d) 5-Hz Ricker Wavelet using Wiener Filter');axis([0 1.1 -1 1]);
xlabel('Time (seconds)');ylabel('Amplitude');

% set(gcf,'PaperPositionMode','auto');
print -deps wavelet.eps
% print -dpng wavelet.png -r300
