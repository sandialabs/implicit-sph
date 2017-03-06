clear all;
close all;

% box geometry and put a sphere in the middle 
ref_length      = 1;
box_half_length = 1.0*ref_length;
box_half_width  = 1.0*ref_length;
sphere_radius   = 0.6*ref_length;

% discretization
N = 128;

Nx = round(N*(box_half_width/box_half_length)); 
Ny = N;

box_half_width = Nx*box_half_length/N; 

natoms = Nx*Ny;
printf('effective number of atoms used in the box = %d\n', natoms);

dx = (box_half_width*2)/Nx;
dy = (box_half_length*2)/Ny;

printf('dx, dy = %.8f, %.8f\n', dx,dy);

xhi =  box_half_width;
xlo = -box_half_width;
yhi =  box_half_length;
ylo = -box_half_length;
zhi =  0.5;
zlo = -0.5;

x = xlo : dx : xhi-dx/2;
y = ylo : dy : yhi-dy/2;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 
printf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

water = struct('density',1.0e-21, 'viscosity',1.0e10, 'psi', 0.0);

scaphi = 4.142e-21/1.6e-19;
%solid = struct('density',1.0e-18, 'viscosity',0.0,'psi',-50e-3/scaphi);
solid = struct('density',1.0e-21, 'viscosity',1.0e10, 'psi',-50e-3/scaphi);

% array
density = water.density * ones(natoms,1);                                       
viscosity = water.viscosity * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = water.psi * ones(natoms,1);
eps = 6.949320148e-10*ones(natoms,1);

% solid type
type((X-xlo).^2 + (Y-ylo).^2 < sphere_radius^2) = 2;
type((X-xhi).^2 + (Y-ylo).^2 < sphere_radius^2) = 2;
type((X-xhi).^2 + (Y-yhi).^2 < sphere_radius^2) = 2;
type((X-xlo).^2 + (Y-yhi).^2 < sphere_radius^2) = 2;

density(type==2) = solid.density;
psi(type==2) = solid.psi;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('charged-membrane-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Pressure driven flow through charged membrane 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%8d %2d % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e \n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %.8f #--datafile\' charged-membrane-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
