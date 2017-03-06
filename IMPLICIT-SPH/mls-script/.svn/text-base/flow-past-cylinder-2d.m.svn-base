clear all;
close all;

fid = fopen("flow-past-cylinder-2d.out");
tmp = fscanf(fid, '%f');
fclose(fid); 

N = tmp(1);
lx = tmp(2);
ly = tmp(3);

r   = tmp(4);
off = tmp(5);

dx = ly/N;

xhi =  lx+off
xlo = -lx+off
yhi =  ly+dx*0.5;
ylo = -ly-dx*0.5;
zhi =  dx*0.5;
zlo = -dx*0.5;

bx = [xlo+dx/2 : dx : xhi]';
n = size(bx,1);

by = ly*ones(n,1); 
bz = zeros(n,1);

dtheta = dx/r/2;
ntheta = 2*pi/dtheta;
theta = linspace(-pi, pi, ntheta)'; %%[-pi : dtheta : pi-dtheta/2]';
theta(end) = [];
n = size(theta,1);

cx = cos(theta)*r;
cy = sin(theta)*r;
cz = zeros(n,1);

xx = [bx;bx;cx];
yy = [-by;by;cy];
zz = [bz;bz;cz];

%xx = [bx;bx];
%yy = [-by;by];
%zz = [bz;bz];

%xx = [cx];
%yy = [cy];
%zz = [cz];

natoms = size(xx,1);

tag = (1:natoms)';

type = [2*ones(2*size(bx,1),1); 3*ones(size(cx,1),1)];

density = 1.0*ones(natoms,1);                                       
viscosity = 0.1*ones(natoms,1);  
pressure = zeros(natoms,1);
psi = zeros(natoms,1);
eps = ones(natoms,1);

Atoms = [tag, type, density, viscosity, pressure, psi, eps, xx, yy, zz]; 

fid = fopen('flow-past-cylinder-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Flow Past a Cylinder 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', 3);

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

