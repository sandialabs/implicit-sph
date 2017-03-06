clear all;
close all;

% box geometry and put a sphere in the middle 
square_box = pi *2.0;

% discretization
N = 64;

%natoms = (N-1)*(N-1);
%printf('effective number of atoms used in the box = %d\n', natoms);

dx = (square_box)/(N-1);
dy = (square_box)/(N-1);

xhi =  square_box;
xlo =  0.0;
yhi =  square_box;
ylo =  0.0;
zhi =  0.5;
zlo = -0.5;

x = xlo : dx : xhi;
y = ylo : dy : yhi;
z = 0.;

[X, Y, Z] = meshgrid(x,y,z);
X = reshape(X, numel(X), 1);
Y = reshape(Y, numel(Y), 1);
Z = reshape(Z, numel(Z), 1);

% actual number of atoms used for this model
natoms = size(X,1); 
printf('actual number of atoms used in the box including solids = %d\n', natoms);

tag = (1:natoms)';
type = ones(natoms,1);

density = 1.0 * ones(natoms,1);                                       
viscosity = 0.1 * ones(natoms,1);  
pressure = 0.0 * ones(natoms,1);
psi = 0.0 * ones(natoms,1);
eps = density;

% solid type
type((X==xlo)|(X==xhi)|(Y==ylo)|(Y==yhi))=2;

Atoms = [tag, type, density, viscosity, pressure, psi, eps, X, Y, Z]; 

fid = fopen('poisson-boundary-2d.data', 'wt' );
fprintf(fid,'LAMMPS ISPH data file for box flow 2D\n');
fprintf(fid,'%d atoms\n', natoms);
fprintf(fid,'%d atom types\n', size(unique(type),1));

fprintf(fid,'%f %f xlo xhi\n', xlo-dx/2, xhi+dx/2);
fprintf(fid,'%f %f ylo yhi\n', ylo-dy/2, yhi+dy/2);
fprintf(fid,'%f %f zlo zhi\n', zlo, zhi);

fprintf(fid,'\nAtoms\n\n');

fprintf(fid,'%d %d %f %f %f %f %f %f %f %f\n', Atoms'); 

fclose(fid);

cmd=sprintf("sed -i \'/--datafile/c\\variable           dx equal %f #--datafile\' poisson-boundary-2d.lmp", max(dx,dy));
[r,s] = system(cmd);
disp(s);
