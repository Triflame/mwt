clear;clc

nx = 499;
nz = 101;
dx = 5;

fid = fopen('vel4','r');
v = fread(fid, nx*nz,'float32');
fclose(fid);
v = reshape(v, nz, nx);
v = v(:,199:299);

nx = 101;
nz = 101;
x = (0:nx-1)*dx;
z = (0:nz-1)*dx;
colormap(gray);imagesc(x,z,v);colorbar;title('(c) TRT Tomogram');
xlabel('X (m)');ylabel('Z (m)');

print -dpng figure2c.png
exit

