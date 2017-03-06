clear all;
close all;

box_length = 2.0*pi;
box_height = 2.0*pi;

% discretization
% ----------------------------------------------------------------
N = 128;

Nx = round(N*(box_length/box_height));
Ny = N;

dx = (box_length)/Nx;
dy = (box_height)/Ny;

fprintf('dx, dy = %.8f, %.8f\n', dx,dy);

offset = 4;

xhi =  box_length + offset*dx;
xlo =  0.0 - offset*dx;
yhi =  box_height + offset*dy;
ylo =  0.0 - offset*dy;
zhi =  0.5;
zlo = -0.5;

Nx = Nx + 2*offset;
Ny = Ny + 2*offset;
% ----------------------------------------------------------------
x = linspace(xlo, xhi, Nx);
y = linspace(ylo, yhi, Ny);
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

natoms = size(X,1); 
type = ones(natoms,1);

eps = dx/1000;

type(X < 0.0 | X > box_length) = 2;
type(Y < 0.0 | Y > box_height) = 2;

% ----------------------------------------------------------------
natoms = size(X,1);
fprintf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 1.0 * ones(natoms,1);

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('shift-test-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for shift-test 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo-dy, xhi+dy);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo-dy, yhi+dy);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf('sed -i \''/--datafile/c\\variable           dx equal %.8f #--datafile\'' shift-test-2d.lmp', max(dx,dy));
[r,s] = system(cmd);
disp(s);
