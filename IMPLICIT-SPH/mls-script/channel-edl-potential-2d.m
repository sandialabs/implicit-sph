clear all;
close all;

% geometry (symmetry w.r.t. y axis)
%          |------------>| channel_length
%          |------->|      potential_zone
% ------------------------
%                          channel_radius
%       orign-------------
%           
% ------------------------
channel_length = 0.2;
channel_radius = 1.0;
potential_zone = 1.2*channel_length;

% discretization
N = 256;

Nx = round(N*(channel_length/channel_radius));
Ny = N;

channel_length = Nx*channel_radius/N; 

dx = (channel_length*2)/Nx;
dy = (channel_radius*2)/Ny;

printf('dx, dy = %.8f, %.8f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
nlayer = 10;

xhi =  channel_length;
xlo = -channel_length;;
yhi =  channel_radius;
ylo = -channel_radius;
zhi =  0.5;
zlo = -0.5;

x = xlo+dx/2 : dx : xhi;
y = ylo      : dy : yhi;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 
printf('actual number of atoms used in the channel including boundary = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 1.0 * ones(natoms,1);

% boundary type
type(Y.^2 >= channel_radius^2) = 2;

psi(type==2 & (X.^2 < potential_zone^2)) = 1.0;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('channel-edl-potential-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for channel flow 2D with PoissonBoltzmann\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo-dy, yhi+dy);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %.8f #--datafile\' channel-edl-potential-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
