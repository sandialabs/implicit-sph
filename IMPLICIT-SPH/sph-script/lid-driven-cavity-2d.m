clear all;
close all;

box_half_x = 0.5;
box_half_y = 0.5;

% discretization
N = 128;

Nx = N*(box_half_x/box_half_y);
Ny = N;

natoms = Nx*Ny;
printf('effective number of atoms used in the cavity = %d\n', natoms);

dx = (box_half_x*2)/Nx;
dy = (box_half_y*2)/Ny;

printf('dx, dy = %f, %f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
delta = max([dx,dy]);                                       
nn = 10;  

wall_half_x = delta*nn + box_half_x;
wall_half_y = delta*nn + box_half_y;

xhi =  wall_half_x;
xlo = -wall_half_x;
yhi =  wall_half_y;
ylo = -wall_half_y;
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
printf('actual number of atoms used in the cavity including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.01 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 1.0 * ones(natoms,1);

% solid type
type(X.^2 > box_half_x^2) = 2;
type(Y.^2 > box_half_y^2) = 2;

tmp = type;

type(X.^2 <= (box_half_x - 4*delta)^2 & Y > box_half_y & tmp == 2) = 3;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('lid-driven-cavity-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for lid driven cavity\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo, xhi);
fprintf(fid,'%f %f ylo yhi\n', ylo, yhi);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile dx/c\\variable           dx equal %f #--datafile dx\' lid-driven-cavity-2d.lmp", delta);
[r,s] = system(cmd);
disp(s);

cmd=sprintf("sed -i \'/--datafile xx/c\\variable           xx equal %f #--datafile xx\' lid-driven-cavity-2d.lmp", box_half_x);
[r,s] = system(cmd);
disp(s);

cmd=sprintf("sed -i \'/--datafile yy/c\\variable           yy equal %f #--datafile yy\' lid-driven-cavity-2d.lmp", box_half_y);
[r,s] = system(cmd);
disp(s);

