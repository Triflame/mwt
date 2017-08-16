clear;clc

nx = 499;
nz = 101;
dx = 5;

fid = fopen('vel1','r');
v = fread(fid, nx*nz,'float32');
fclose(fid);
v = reshape(v, nz, nx);

imagesc(v);
