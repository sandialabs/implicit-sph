clear all;
close all;

channel_half_x = 0.1;
channel_half_y = 1.0;
channel_half_z = 1.0;

% discretization
N = 128;

Nx = N*(channel_half_x/channel_half_y);
Ny = N*(channel_half_y/channel_half_y);
Nz = N*(channel_half_z/channel_half_y);

natoms = Nx*Ny*Nz;
printf('effective number of atoms used in the channel = %d\n', natoms);

dx = (channel_half_x*2)/Nx;
dy = (channel_half_y*2)/Ny;
dz = (channel_half_z*2)/Nz;

printf('dx, dy, dz = %f, %f, %f\n', dx,dy,dz);

% solid zone needs at least 5*dy (solid zone)
wall_half_y = max([dx,dy,dz])*10 + channel_half_y;
wall_half_z = max([dx,dy,dz])*10 + channel_half_z;

xhi =  channel_half_x;
xlo = -channel_half_x;
yhi =  wall_half_y;
ylo = -wall_half_y;
zhi =  wall_half_z;
zlo = -wall_half_z;

x = xlo : dx : xhi-dx/2;
y = ylo : dy : yhi-dy/2;
z = zlo : dz : zhi-dz/2;


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
eps = 1.0 * ones(natoms,1);

% solid type
type(Y.^2 > channel_half_y^2) = 2;
type(Z.^2 > channel_half_z^2) = 2;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('channel-moving-wall-3d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for channel flow 3D with moving wall\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo, xhi);
fprintf(fid,'%f %f ylo yhi\n', ylo, yhi);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %f #--datafile\' channel-moving-wall-3d.lmp", max([dx,dy,dz]));
[r,s] = system(cmd);
disp(s);
