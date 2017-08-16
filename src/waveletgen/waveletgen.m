clear;clc

fpeak = 10;
dt = 0.001;
nt = 2000;
f = ricker(fpeak,dt);
f = [f(:); zeros(nt-length(f),1)];
plot(f);

