clear all;
close all;

width = 2*pi;
height = 2*pi;

% discretization
N = 96;

Nx = N;
Ny = N;

natoms = Nx*Ny;
printf('effective number of atoms used in the channel = %d\n', natoms);

dx = (width)/Nx;
dy = (height)/Ny;

printf('dx, dy = %f, %f\n', dx,dy);


xhi =  width/2;
xlo = -width/2;
yhi =  height/2;
ylo = -height/2;
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
eps = sqrt(1+X.^2+Y.^2);


% solid type
% type(Y.^2 > channel_radius^2) = 2;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('poisson-boltzmann-dielectric-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file poisson-boltzmann-dielectric\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %.8f #--datafile\' poisson-boltzmann-dielectric-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
