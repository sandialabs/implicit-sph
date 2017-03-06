clear all;
close all;

% channel
channel_width = 1.5;
channel_height = 0.5;

% spinner
spinner_width = 0.25;
spinner_height = 0.12;

% discretization
N = 128;

Nx = round(N*(channel_width/channel_height));
Ny = N;

channel_width = Nx*channel_height/N;

natoms = Nx*Ny;
fprintf('effective number of atoms used in the channel = %d\n', natoms);

dx = (channel_width*2)/Nx;
dy = (channel_height*2)/Ny;

fprintf('dx, dy = %f, %f\n', dx,dy);

% solid zone needs at least 5*dy (solid zone)
outer_wall_height = max(dx,dy)*10 + channel_height;;

xhi =  channel_width;
xlo = -channel_width;;
yhi =  outer_wall_height;
ylo = -outer_wall_height;
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
fprintf('actual number of atoms used in the channel including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);

% solid type
type(Y.^2 > channel_height^2) = 2;
type(X.^2 < spinner_width^2 & Y.^2 < spinner_height^2) = 2;

Atoms = [tag, type, density, viscosity, pressure, psi, X, Y, Z]; 

fid = fopen('spinner-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for Spinner 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%.8f %.8f xlo xhi\n', xlo, xhi);
fprintf(fid,'%.8f %.8f ylo yhi\n', ylo, yhi);
fprintf(fid,'%.8f %.8f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n', Atoms'); 

fclose(fid);

cmd=sprintf('sed -i \''/--datafile/c\\variable           dx equal %.8f #--datafile\'' spinner-2d.lmp', max(dx,dy));
[r,s] = system(cmd);
disp(s);
