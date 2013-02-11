% Hyposearch input job file
%
% Job template: 01: use IBCAO bathymetry for seafloor

event   = 'BATH_MULTI_LAYER';
jobname = 'BM1';

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
%
% Will be ignored when region_auto is set
margins     = [60000 60000 60000 60000];


% Velocity and depth specifications
velp = [1500; 3500; 5400; 7400; 8100; 8250];
vels = [   0; 2020; 3120; 4500; 4680; 4760];


% Interfaces, specify to use fixed seafloor (0), IBCAO bathymetry (1) or a custom
% seafloor (2).
seafloortype = 1;

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
region_auto = false;

% if manual: boundaries of box of region of interest, will be gridded
% according to makegrid.
% load nodes
quakes = load ('quakes.coor');
stations = load ('stations.coor');
nodes = cat(1, stations, quakes);
margins = 10000;
l = min(nodes(:,1));
r = max(nodes(:,1));
u = max(nodes(:,2));
d = min(nodes(:,2));

rxmin = l - margins;
rxmax = r + margins;
rymin = u + margins;
rymax = d - margins;
rzmin = 5000; rzmax = 7500;

% number of grid points
region_size_auto = false;
rnx = 20;
rny = 20;
rnz = 10;

%% Interfaces
ninterfaces     = 4; % including seafloor
interfaces_add  = true;
interface_tris_add = cell(ninterfaces-1,1);

depths=[4300 5650 7000 10000 ]; %40000];
for d=2:length(depths)
  xmin = rxmin; xmax = rxmax;
  ymin = rymin; ymax = rymax;
  
  xx = [xmin xmax];
  yy = [ymin ymax];
  zz = -depths(d) .* ones(2);
  [yg, xg] = meshgrid (yy, xx);
  yg = yg(:); xg = xg(:);
  dt = delaunay(xg, yg);

  int = TriRep (dt, xg, yg, -zz(:));
  interface_tris_add{d-1} = int;
end
