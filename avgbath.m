function [a, s] = avgbath (nodes)
% get avarage depth of bath within convex hull of nodes

ncid = netcdf.open ('ibcao_det.grd', 'NC_NOWRITE'); % will be created by makehyposearch.sh
x = netcdf.getVar(ncid, 0); % longitude
y = netcdf.getVar(ncid, 1); % latitude
z = netcdf.getVar(ncid, 2); % bathymetry (depth is negative)
z = -double(z);

dx = 500; dy = 500; % meters

netcdf.close (ncid);

[yg, xg] = meshgrid(y, x);
dt = delaunay (xg(:), yg(:));
bath = TriRep (dt, xg(:), yg(:), z(:));

nodes = nodes(:,1:2);
orignodes = nodes;

margin = 5000;
mx = mean (nodes(:,1));
my = mean (nodes(:,2));
nodes(nodes(:,1)>mx,1) = nodes(nodes(:,1)>mx,1) + margin;
nodes(nodes(:,1)<mx,1) = nodes(nodes(:,1)<mx,1) - margin;
nodes(nodes(:,2)>my,2) = nodes(nodes(:,2)>my,2) + margin;
nodes(nodes(:,2)<my,2) = nodes(nodes(:,2)<my,2) - margin;

X = bath.X(:,1:2);

dt = DelaunayTri (nodes);

triids = pointLocation (dt, X);
bx = bath.X(~isnan(triids),:);
dtb = delaunay (bx(:,1:2));
dt  = TriRep (dtb, bx);

doplot = true;
if (doplot)
  figure(1); clf('reset');

  trisurf (dt);
  xlabel ('x [m]');
  ylabel ('y [m]');
  zlabel ('Depth [m]');
  set (gca, 'ZDir', 'reverse');
  title ('Sample bathymetry used to get avarage depth');
  hold on;
  plot3 (orignodes(:,1), orignodes(:,2), repmat([0], size(nodes,1), 1), 'rx');
  plot3 (nodes(:,1), nodes(:,2), repmat([0], size(nodes,1), 1), 'kx');
  legend ('Bathymetry', 'Orignal nodes', 'Nodes with margin');
end

a = mean(bx(:,3));
s = std(bx(:,3));

fprintf ('Average: %g [m], std. deviation: %g [m], max %g [m], min: %g [m]\n', a, s, max(bx(:,3)), min(bx(:,3)));

end