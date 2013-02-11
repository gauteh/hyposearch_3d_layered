function [rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz)
  
%fprintf ('sg: setting up target grid..\n');
%% Grid size, accuracy of arrival rays

%% Build grid
% margin to match a grid point, 1.0 is full cell
ep_x = 1.0;
ep_y = 1.0;
ep_z = 1.0;

% create background mesh and crop it to the hull of R
rminx = min(R_tri.X(:,1));
rmaxx = max(R_tri.X(:,1));
rminy = min(R_tri.X(:,2));
rmaxy = max(R_tri.X(:,2));
rminz = min(R_tri.X(:,3));
rmaxz = max(R_tri.X(:,3));

rdx = (rmaxx - rminx) / (rnx);
rdy = (rmaxy - rminy) / (rny);
rdz = (rmaxz - rminz) / (rnz);

rgxv  = linspace (rminx + rdx/2, rmaxx - rdx/2, rnx);
rgyv  = linspace (rminy + rdy/2, rmaxy - rdy/2, rny);
rgzv  = linspace (rminz + rdz/2, rmaxz - rdz/2, rnz);

[rgx, rgy, rgz] = meshgrid (rgxv, rgyv, rgzv); % grid points, centers

end