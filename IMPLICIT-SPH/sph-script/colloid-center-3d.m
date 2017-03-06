clear all;
close all;

% box geometry and put a sphere in the middle 
box_half_length = pi;
box_half_width  = pi;
box_half_height = pi;

sphere_radius   = 1.5;

% discretization
N = 64;

Nx = N*(box_half_width/box_half_length);
Ny = N*(box_half_length/box_half_length);
Nz = N*(box_half_height/box_half_length);

natoms = Nx*Ny*Nz;
printf('effective number of atoms used in the box = %d\n', natoms);

dx = (box_half_width*2)/Nx;
dy = (box_half_length*2)/Ny;
dz = (box_half_height*2)/Nz;

printf('dx, dy, dz = %f, %f, %f\n', dx,dy,dz);

xhi =  box_half_width;
xlo = -box_half_width;
yhi =  box_half_length;
ylo = -box_half_length;
zhi =  box_half_height;
zlo = -box_half_height;

x = xlo : dx : xhi-dx/2;
y = ylo : dy : yhi-dy/2;
z = zlo : dz : zhi-dz/2;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 
printf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 1.0 * ones(natoms,1);

% solid type
type(X.^2 + Y.^2 + Z.^2 < sphere_radius^2) = 2;

density(type==2) = 1000.0;
psi(type==2) = 0.0;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('colloid-center-3d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for box flow 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo, xhi);
fprintf(fid,'%f %f ylo yhi\n', ylo, yhi);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %f #--datafile\' colloid-center-2d.lmp", max([dx,dy,dz]));
[r,s] = system(cmd);
disp(s);
