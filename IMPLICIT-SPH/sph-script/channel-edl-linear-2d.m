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

Nx = N*(channel_length/channel_radius);
Ny = N;

%natoms = Nx*Ny;
%printf('effective number of atoms used in the channel = %d\n', natoms);

dx = (channel_length*2)/Nx;
dy = (channel_radius*2)/Ny;

delta = max(dx,dy);

printf('dx, dy = %f, %f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
nlayer = 10;
outer_wall_radius = delta*nlayer + channel_radius;;

xhi =  channel_length;
xlo = -channel_length;;
yhi =  outer_wall_radius;
ylo = -outer_wall_radius;
zhi =  0.5;
zlo = -0.5;

x = xlo : dx : xhi-dx/2;

use_regular = 1;
if use_regular == 1
  % regular distribution
  y = ylo : dy : yhi-dy/2;
else
  % anisotropic distribution
  nlayer = 20;
  delta2 = delta/10;
  ya = ylo:delta2:ylo+delta*nlayer;
  yb = yhi-delta*nlayer:delta2:yhi;
  y =unique([ya ylo+delta*nlayer:delta:yhi-delta*nlayer yb]);
end

z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 
printf('actual number of atoms used in the channel including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 0.02 * ones(natoms,1);

% solid type
type(Y.^2 > channel_radius^2) = 2;
psi(type==2 & (X.^2 < potential_zone^2)) = 1.0;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('channel-edl-linear-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for channel flow 2D with PoissonBoltzmann\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo, xhi);
fprintf(fid,'%f %f ylo yhi\n', ylo, yhi);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %f #--datafile\' channel-edl-linear-2d.lmp", delta);
[r,s] = system(cmd);
disp(s);
