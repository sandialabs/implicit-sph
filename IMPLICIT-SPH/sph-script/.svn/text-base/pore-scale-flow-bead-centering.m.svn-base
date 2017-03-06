clear all;
close all;

%%%%%%%%%%%%%%%%
% use SI units %
%%%%%%%%%%%%%%%%

% bead centeroids file open
file_id = fopen('pore-scale-flow-bead-centeroids-3d.dat','r');
bead_centeroids = fscanf(file_id,'%f %f %f',[3 Inf]);
fclose(file_id);

% outline geometry
bead_position_scale_to_si = 0.002;
bead_radius = 2.5e-4; %% 2.5e-4;

% scaling geometry
bead_centeroids = bead_centeroids * bead_position_scale_to_si;

% % original problem
% %%%%%%%%%%%%%%%%%%
n = 96; %% 192, 384;

cylinder_eff_half_length = 12.8e-3/2;
cylinder_half_buffer = 1.58e-3/2;

cylinder_length = cylinder_eff_half_length*2 + cylinder_half_buffer*2;

cylinder_radius = 4.4e-3;
cylinder_eff_radius = cylinder_radius; 

% discretization
nx = n;
nz = n;
ny = n*(cylinder_length/2/cylinder_radius);

natoms = nx*ny*nz;
fprintf('effective number of atoms used in the box = %d\n', natoms);

dx = cylinder_radius*2/nx;
dz = cylinder_radius*2/nz;
dy = cylinder_length  /ny;

fprintf('dx, dy, dz = %f, %f, %f\n', dx,dy,dz);

wall = 6*max([dx dy dz]);

% put the geometry in the center
for i=1:3
    min_bead_centeroids(i) = min(bead_centeroids(i,:));
    max_bead_centeroids(i) = max(bead_centeroids(i,:));
end
offsets = (min_bead_centeroids + max_bead_centeroids)/2;

fprintf('min bead centeroids = %f %f %f\n', min_bead_centeroids-offsets);
fprintf('max bead centeroids = %f %f %f\n', max_bead_centeroids-offsets);

%% from alex
offsets = [ (7.8498e-4 + 9.6050e-3)/2 ...
            (7.8000e-4 + 1.3600e-2)/2 - cylinder_half_buffer ...
            (5.7998e-4 + 9.3700e-3)/2 ];
for i=1:3
    bead_centeroids(i,:) = bead_centeroids(i,:) - offsets(i);
end

% lammps box 
cylinder_half_buffer
wall 
xhi =  cylinder_radius
xlo = -cylinder_radius
zhi =  xhi
zlo =  xlo
yhi =  cylinder_length/2
ylo = -cylinder_length/2

format long;
bead_centeroids'