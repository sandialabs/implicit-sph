clear all;
close all;

% channel
channel_length = 1.5;
channel_height = 0.5;

% spinner
spinner_length = 0.1;
spinner_height = 0.25;

% discretization
% ----------------------------------------------------------------
N = 64;

Nx = round(N*(channel_length/channel_height));
Ny = N;

channel_length = Nx*channel_height/N;

natoms = Nx*Ny;

dx = (channel_length*2)/Nx;
dy = (channel_height*2)/Ny;

fprintf('dx, dy = %.8f, %.8f\n', dx,dy);

xhi =  channel_length*1.7;
xlo = -channel_length*0.3;
yhi =  channel_height;
ylo = -channel_height;
zhi =  0.5;
zlo = -0.5;

% ----------------------------------------------------------------
x = xlo+dx/2 : dx : xhi;
y = ylo      : dy : yhi;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

natoms = size(X,1); 
type = ones(natoms,1);

type(X.^2 < (spinner_length + dx/2)^2 & Y.^2 < (spinner_height + dy/2)^2) = 2;
cond = (type(:) == 2);

X(cond,:) = [];
Y(cond,:) = [];
Z(cond,:) = [];
type(cond,:) = [];

% ----------------------------------------------------------------
bx = linspace(-spinner_length, spinner_length, 2*spinner_length/dx);
by = linspace(-spinner_height, spinner_height, 2*spinner_height/dy);
bz = 0;

[BX, BY, BZ] = meshgrid(bx,by,bz);
BX = reshape(BX, numel(BX), 1);
BY = reshape(BY, numel(BY), 1);
BZ = reshape(BZ, numel(BZ), 1);

nboundary = size(BX,1);
btype = 2*ones(nboundary,1);

btype(BX.^2 < spinner_length^2 & BY.^2 < spinner_height^2) = 3;
btype((BX.^2 + BY.^2) >= (spinner_length^2 + spinner_height^2)) = 3;
cond = (btype(:) == 3);

BX(cond,:) = [];
BY(cond,:) = [];
BZ(cond,:) = [];
btype(cond,:) = [];

% ----------------------------------------------------------------
X = [ X ; BX ];
Y = [ Y ; BY ];
Z = [ Z ; BZ ];

type = [ type ; btype ];

type(Y.^2 >= channel_height^2) = 2;

% ----------------------------------------------------------------
natoms = size(X,1);
fprintf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.005 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = density;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('spinner-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Spinner 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo-dy/2, yhi+dy/2);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf('sed -i \''/--datafile/c\\variable           dx equal %.8f #--datafile\'' spinner-2d.lmp', max(dx,dy));
[r,s] = system(cmd);
disp(s);
