% Hyposearch input job file
%
% Job template: 02: Use flat bathymetry at 3000 meters.

event   = 'FLAT_SYNTH';
jobname = 'FS1';

% Specify which phases to solve for which station,
% same order as stations.coor. Phases in order:
% P SP M MM

solvephases = [1 1 1 1;
               1 1 1 1;
               1 1 1 1;
               1 1 1 1];

% Specify which phases to use in RMS calculation
usephases   = [1 1 1 1;
               1 1 1 1;
               1 1 1 1;
               1 1 1 1];

% Specify margins of model in meters (on UPS grid) from the corners of 
% the array triangle, format: top, down, left, right.
margins     = [1000 1000 1000 1000];


% Velocity and depth specifications
velp = [1500; 5800];
vels = [   0; 3200];

% Interfaces, specify to use fixed seafloor (0), IBCAO bathymetry (1) or a custom
% seafloor (2).
seafloortype = 0;

% fixed seafloor
seafloor_depth = 3000; % positive depth


% custom, create a triangulated surface which extends the full bounds of
% the final model (coordinates should be UPS), vertices:
I = [   0     0   3000;
        0  2000   3000;
      500     0   3000;
      500  2000   3000;
     1000     0   5000;
     1000  2000   5000;
     1500     0   5000;
     1500  2000   5000;
     2000     0   3000;
     2000  2000   3000;
     2500     0   3000;
     2500  2000   3000;
     ];
   
x = I(:,1);
y = I(:,2);
z = I(:,3);
dt = delaunay(x, y);
seafloor = TriRep (dt, x, y, z);
clear I x y z dt;

% Region of interest (optional, may be automatically calculated)
region_auto = true;

% if manual: boundaries of box of region of interest, will be gridded
% according to makegrid.
% rxmin = -10000; rxmax = 10000;
% rymin = -10000; rymax = 10000;
% rzmin = 4000; rzmax = 7000;

% number of grid points
region_size_auto = false;
rnx = 20;
rny = 20;
rnz = 10;

%% Interfaces
ninterfaces = 1; % including seafloor
interfaces_add = false;
% interface_tris_add = cell(ninterfaces-1,1);

