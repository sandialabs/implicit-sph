clear all;
close all;

% box
box_length_left = pi;
box_length_right = pi;

box_length = box_length_left + box_length_right;

box_height_lower = pi/2.0;
box_height_upper = pi/2.0;

box_height = box_height_lower + box_height_upper;

% cylinder
cyl_radius = 1.02;

% discretization
% ----------------------------------------------------------------
N = 40;

Nx = round(N*(box_length/box_height));
Ny = N;

box_length = Nx*box_height/N;
box_length_left = box_length - box_length_right;

natoms = Nx*Ny;

dx = (box_length)/Nx;
dy = (box_height)/Ny;

fprintf('dx, dy = %.8f, %.8f\n', dx,dy);

xhi =  box_length_right;
xlo = -box_length_left;
yhi =  box_height_upper;
ylo = -box_height_lower;
zhi =  0.5;
zlo = -0.5;

% ----------------------------------------------------------------
x = xlo : dx : xhi;
y = ylo : dy : yhi;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

natoms = size(X,1); 
type = ones(natoms,1);

type(X.^2 + Y.^2 < (cyl_radius + dy/2)^2) = 2; 
cond = (type(:) == 2);

X(cond,:) = [];
Y(cond,:) = [];
Z(cond,:) = [];
type(cond,:) = [];

% ----------------------------------------------------------------
delta = dy/cyl_radius;                                       
theta = [0.0 : delta : 2*pi - delta/2]';        
nboundary = size(theta,1);                                       
                                       
BX = cyl_radius*cos(theta);                                       
BY = cyl_radius*sin(theta);                                       
BZ = 0.0 * ones(nboundary,1);                                       
btype = 2 * ones(nboundary,1);

% ----------------------------------------------------------------
X = [ X ; BX ];
Y = [ Y ; BY ];
Z = [ Z ; BZ ];

type = [ type ; btype ];

% ----------------------------------------------------------------
type(X >=  box_length_right) = 2;
type(X <= -box_length_left) = 2;

type(Y >=  box_height_upper) = 2;
type(Y <= -box_height_lower) = 2;

% ----------------------------------------------------------------
natoms = size(X,1);
fprintf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = density;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('unit-test-compact-poisson-boundary-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Unit test compact poisson 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo-dx/2, xhi+dx/2);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo-dy/2, yhi+dy/2);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf('sed -i \''/--datafile/c\\variable           dx equal %.8f #--datafile\'' unit-test-compact-poisson-boundary-2d.lmp', max(dx,dy));
[r,s] = system(cmd);
disp(s);
