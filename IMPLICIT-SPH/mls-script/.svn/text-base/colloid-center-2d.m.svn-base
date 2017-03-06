clear all;
close all;

% box geometry and put a sphere in the middle 
box_half_length = pi;
box_half_width  = pi;
sphere_radius   = 1.5;

% discretization
N = 128;

Nx = round(N*(box_half_width/box_half_length));
Ny = N;
                                              
box_half_width = Nx*box_half_length/N;

dx = (box_half_width*2)/Nx;
dy = (box_half_length*2)/Ny;

xhi =  box_half_width;
xlo = -box_half_width;
yhi =  box_half_length;
ylo = -box_half_length;
zhi =  0.5;
zlo = -0.5;

x = xlo+dx/2 : dx : xhi;
y = ylo+dy/2 : dy : yhi;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 

tag = (1:natoms)';
type = ones(natoms,1);

% exclude solid types
type(X.^2 + Y.^2 < (sphere_radius + dy/3)^2) = 2;
cond = (type(:) == 2);

X(cond,:) = [];
Y(cond,:) = [];
Z(cond,:) = [];
type(cond,:) = [];

delta = dy/sphere_radius;
theta = [0.0 : delta : 2*pi - delta/2]';
nboundary = size(theta,1);

BX = sphere_radius*cos(theta);
BY = sphere_radius*sin(theta);
BZ = 0.0 * ones(nboundary,1);
BT = 2 * ones(nboundary,1);

X = [ X ; BX ];
Y = [ Y ; BY ];
Z = [ Z ; BZ ];

type = [ type ; BT ];

natoms = size(X,1); 
printf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = density;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('colloid-center-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for box flow 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %.8f #--datafile\' colloid-center-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
