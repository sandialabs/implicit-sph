clear all;
close all;

channel_half_length = 0.5;
channel_radius = 0.5;

% discretization
N = 128;

Nx = round(N*(channel_half_length/channel_radius));
Ny = N;

channel_length = Nx*channel_radius/N;

%natoms = Nx*Ny;
%printf('effective number of atoms used in the channel = %d\n', natoms);

dx = (channel_half_length*2)/Nx;
dy = (channel_radius*2)/Ny;

delta = max(dx,dy);

printf('dx, dy = %f, %f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
nlayer = 10;
outer_wall_radius = delta*nlayer + channel_radius;;

xhi =  2*channel_half_length;
xlo =  0.0;
yhi =  outer_wall_radius;
ylo = -outer_wall_radius;
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

use_case_one = 1;
if use_case_one == 1
  %% case 1 
  psi(type==2 & (X.^2 <  (channel_half_length)^2)) = 1.0;
  psi(type==2 & (X.^2 >= (channel_half_length)^2)) = -1.0;
else
  %% case 2 
  psi(type==2 & (X.^2 <  (channel_half_length)^2) & Y>0) =  1.0;
  psi(type==2 & (X.^2 <  (channel_half_length)^2) & Y<0) = -1.0;
  psi(type==2 & (X.^2 >= (channel_half_length)^2) & Y>0) = -1.0;
  psi(type==2 & (X.^2 >= (channel_half_length)^2) & Y<0) =  1.0;
end

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('channel-edl-alternate-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for channel flow 2D with PoissonBoltzmann\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo, xhi);
fprintf(fid,'%f %f ylo yhi\n', ylo, yhi);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %f #--datafile\' channel-edl-alternate-2d.lmp", delta);
[r,s] = system(cmd);
disp(s);
