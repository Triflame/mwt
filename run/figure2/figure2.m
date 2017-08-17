clear;clc;

nx = 101;
nz = 101;
dx = 5;

x = (0:nx-1)*dx;
z = (0:nz-1)*dx;

% true velocity model
fid = fopen('vel_true.dat','r');
v = fread(fid,nx*nz,'float32');
fclose(fid);
vtrue = reshape(v, nz, nx);

% initial velocity model
fid = fopen('vel_grad.dat','r');
v = fread(fid,nx*nz,'float32');
fclose(fid);
vgrad = reshape(v, nz, nx);

nx = 501;

% traveltime tomogram
fid = fopen('vel_trt.dat','r');
v = fread(fid,nx*nz,'float32');
fclose(fid);
vtrt = reshape(v, nz, nx);
vtrt = vtrt(:,201:301);

% early-arrival waveform tomogram
fid = fopen('vel_ewt.dat','r');
v = fread(fid,nx*nz,'float32');
fclose(fid);
vewt = reshape(v, nz, nx);
vewt = vewt(:,201:301);

% multiscale waveform tomogram
fid = fopen('vel_mwt.dat','r');
v = fread(fid,nx*nz,'float32');
fclose(fid);
vmwt = reshape(v, nz, nx);
vmwt = vmwt(:,201:301);

colormap(gray)
left = 650;
top = -50;
subplot(321);imagesc(x,z,vtrue,[1500 3000]);xlabel('X (m)');ylabel('Z (m)');colorbar;
title('(a) True Velocity Model');text(left,top,'m/s');
subplot(322);imagesc(x,z,vgrad,[1500 3000]);xlabel('X (m)');ylabel('Z (m)');colorbar;
title('(b) Initial Velocity Model');text(left,top,'m/s');
subplot(323);imagesc(x,z,vtrt,[1500 3000]);xlabel('X (m)');ylabel('Z (m)');colorbar;
title('(c) TRT Tomogram');text(left,top,'m/s');
subplot(324);imagesc(x,z,vewt,[1500 3000]);xlabel('X (m)');ylabel('Z (m)');colorbar;
title('(d) EWT Tomogram');text(left,top,'m/s');
subplot(325);imagesc(x,z,vmwt,[1500 3000]);xlabel('X (m)');ylabel('Z (m)');colorbar;
title('(e) MWT Tomogram');text(left,top,'m/s');

print -dpng figure2.png

exit

