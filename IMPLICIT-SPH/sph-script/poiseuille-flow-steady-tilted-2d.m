clear all;
close all;

channel_length = 0.2;
channel_radius = 0.4;

theta = pi/6;

Q = [cos(theta) -sin(theta); sin(theta), cos(theta)];

domain_max = max((Q* [channel_length, channel_length, -channel_length,  -channel_length; -channel_radius, channel_radius, channel_radius, -channel_radius])');
x_max= domain_max(1);
y_max=domain_max(2);

% discretization
N = 64;

Nx = round(N*(x_max/y_max));
Ny = N;

x_max = Nx*y_max/N;

natoms = Nx*Ny;
printf('effective number of atoms used in the channel = %d\n', natoms);

dx = (x_max*2)/Nx;
dy = (y_max*2)/Ny;

printf('dx, dy = %f, %f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
xhi = max(dx,dy)*10 + x_max;
yhi = max(dx,dy)*10 + y_max;
xlo = -xhi;
ylo = -yhi;
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
viscosity = 1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = 1.0 * ones(natoms,1);

% solid type
type(abs(Q(1,1)*X + Q(2,1)*Y) > channel_length) = 2;
 type(abs(Q(1,2)*X + Q(2,2)*Y) > channel_radius) = 2;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('poiseuille-flow-steady-tilted-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Poiseuille Flow 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %.8f #--datafile\' poiseuille-flow-steady-tilted-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
