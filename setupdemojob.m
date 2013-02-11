% demo job for hyposearch
xmin = 0; xmax = 2500;
ymin = 0; ymax = 2000;
zmin = 0; zmax = 7000;

bounds = [xmin xmax; ymin ymax; zmin zmax];

velp     = [1500; 5000];  % velocity of layers: water, bedrock, even denser
vels     = [   0; 4000];  % S - velocity


% Seafloor vertices [x, y, z]
I = [ 0    0   3000;
      0  2000   3000;
     500    0   3000;
     500  2000   3000;
     1000    0   5000;
     1000  2000   5000;
     1500   0   5000;
     1500 2000   5000;
     2000   0   3000;
     2000 2000   3000;
     2500   0   3000;
     2500 2000   3000;
     ];
   

x = I(:,1);
y = I(:,2);
z = I(:,3);
dt = delaunay(x, y);
bathymetry = TriRep (dt, x, y, z);


stations = [1200 700 0];


fprintf ('hs: set up region of interest box..\n');
rxmin = 0; rxmax = 2500;
rymin = 0; rymax = 2000;
rzmin = 4000; rzmax = 7000;

% Create box around region of interest
[bxg, byg, bzg] = meshgrid ([rxmin rxmax], [rymin rymax], [rzmin rzmax]);
bxg = bxg(:); byg = byg(:); bzg = bzg(:);

ch  = convhull (bxg, byg, bzg);
R_tri = TriRep (ch, bxg, byg, bzg);
R_layer = 2; % Target box is located in this layer
clear bxg byg bzg box;