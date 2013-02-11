function [rays, r_ray_p] = solvegrid (stations, bounds, velp, vels, interface_tris, R_tri, rnx, rny, rnz, ninterfaces, nlayers, phase, fig)
%
% solvegrid: solve two point problem of stations <-> grid points of region
% of interest (R_tri).
%
% Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-02-11
%
% Output:
%
%   rays:     all rays traced
%   r_ray_p:  grid points with ray index

% Problem:
%
% Find ray path between station and all grid points of interest
% - For all stations
%   - For all relevant phases
%     - For all grid cells within region of interest
%
% Approach (for each station):
%
% 1. Solve initial probe rays for first (simplest) phase
% 2. Identify grid cells solved sufficiently good (ray intersected grid
%    cell close enough to center)
%
% If unsolved cells:
% 3. Interpolate take off conditions for missing cells
% 4. Solve interpolated take off conditions
% 5. Go to 2.
%
% If all solved:
% Solve next phase
%
% When done:
% - Calculate travel time for all phases for all grid cells
%
%

%% Plot model
doplot = true;

if (doplot)
  figure(fig); clf('reset');
  for i = 1:ninterfaces
    trisurf (interface_tris{i});
    hold on;
  end

  set(gca, 'ZDir', 'reverse');
  title (sprintf('Model and rays, Phase: %s', phase));
  xlabel ('x [m]');
  ylabel ('y [m]');
  zlabel ('Depth [m]');

  % Plot region of interest
  trisurf (R_tri, 'FaceAlpha', 0.3, 'FaceColor', 'r');

  % Plot stations
  plot3 (stations(:,1), stations(:,2), stations(:,3), 'k*');
  zlim ([bounds(3,1) bounds(3,2)]);
  
end
%%

%% Initial shots
initmode = 0; % 0 = equidistant arrival points, 1 = take off range
if (initmode == 0)
% calculate equidistant arrival points at seafloor
if (strcmp(phase, 'P') || strcmp(phase, 'SP'))
  nax = 200; nay = 200; % number of points
elseif (strcmp(phase, 'PIPBP'))
  %nax = 50; nay = 50; % number of points
  nax = 200; nay = 200; % number of points
elseif (strcmp(phase, 'PIPBPIPBP'))
  %nax = 50; nay = 50; % number of points
  nax = 200; nay = 200; % number of points
end
%nax = 20; nay = 20;
%nax = 20; nay = 20;

fprintf ('sg: initial shots: cartesian: nax=%d, nay=%d\n', nax, nay);
ax = linspace (bounds(1,1), bounds(1,2), nax);
ay = linspace (bounds(2,1), bounds(2,2), nay);
[axg, ayg] = meshgrid (ax, ay);

Fa = TriScatteredInterp (interface_tris{1}.X(:,1), interface_tris{1}.X(:,2), interface_tris{1}.X(:,3));
azg = Fa(axg, ayg);
axg = axg(:); ayg = ayg(:); azg = azg(:);

if (doplot)
  plot3 (axg, ayg, azg, 'rx');
end

% make direction vector (z = 0)
axg = axg - stations(1,1);
ayg = ayg - stations(1,2);

aa  = [axg ayg azg];
aan = sqrt (sum (abs (aa).^2, 2)); % norm
aa  = aa ./ repmat(aan,1,3);       % normalize

% calculate thetas
thetas   = acos (aa(:,3));
phis     = atan2(aa(:,2), aa(:,1));

elseif (initmode == 1)
  %fprintf ('sg: do initial shots..\n');
  thetaend = (pi/2);
  if (strcmp(phase, 'PIPBP'))
    nth      = 50;  % no of thetas
    nphi     = 360; % no of phis
    thetaend = (pi/2);
  elseif (strcmp(phase, 'PIPBPIPBP'))
    nth  = 50;
    nphi = 360;
    thetaend = (pi/2)/6;
  else
    nth  = 40;
    nphi = 60;
    thetaend = (pi/2)/12*11;
  end

  fprintf ('sg: initial shots: radial: nth=%d, nphi=%d\n', nth, nphi);
  thetas  = linspace (thetaend/nth, thetaend , nth-1);
  phis    = linspace (0, 2*pi, nphi);
  
  [tt, pp] = meshgrid (thetas, phis);
  
  tx = cos(pp) .* sin(tt);
  ty = sin(pp) .* sin(tt);
  tz = cos(tt);

  % add theta = 0 for phi = 0
  thetas = [0; tt(:)];
  phis   = [0; pp(:)];
  tx = [0; tx(:)];
  ty = [0; ty(:)];
  tz = [1; tz(:)];
  aa = [tx ty tz];
end



%% Create box around model
[bxg, byg, bzg] = meshgrid (bounds(1,:), bounds(2,:), bounds(3,:));
bxg = bxg(:); byg = byg(:); bzg = bzg(:);

ch  = convhull (bxg, byg, bzg);
box = TriRep (ch, bxg, byg, bzg);

bvert0 = box.X (box.Triangulation(:,1),:);
bvert1 = box.X (box.Triangulation(:,2),:);
bvert2 = box.X (box.Triangulation(:,3),:);
bnvert = size (bvert0,1);
%trisurf (box, 'FaceAlpha', 0.2);
clear bxg byg bzg box;

% Create surface
[bxg, byg, bzg] = meshgrid (bounds(1,:), bounds(2,:), bounds(3,1));
bxg = bxg(:); byg = byg(:); bzg = bzg(:);

dt      = delaunay (bxg, byg);
box     = TriRep (dt, bxg, byg, bzg);

svert0 = box.X (box.Triangulation(:,1),:);
svert1 = box.X (box.Triangulation(:,2),:);
svert2 = box.X (box.Triangulation(:,3),:);
snvert = size (svert0,1);
%trisurf (box, 'FaceAlpha', 0.8, 'FaceColor', 'y');
clear bxg byg bzg box;
%%

tri_opt = struct('triangle', 'two sided', 'border', 'inclusive', 'eps', 1e-5);

t0 = tic;

%% Set up target grid
[rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);

%inhull  = ~isnan(tsearchn (R_tri.X, delaunay(R_tri.X), [rgx(:) rgy(:) rgz(:)]));
%rb      = rb(inhull,:); % grid points to solve for
%plot3 (rb(:,1), rb(:,2), rb(:,3), 'bx'); drawnow;
gp      = size(rgx(:),1);  % no of grid points

gvertices0 = nan(12 * gp, 3);
gvertices1 = nan(12 * gp, 3);
gvertices2 = nan(12 * gp, 3);
ngvertices = 12 * gp;

% create box for each grid point (duplicate standard box)
xx = [-rdx/2 rdx/2] .* ep_x;
yy = [-rdy/2 rdy/2] .* ep_y;
zz = [-rdz/2 rdz/2] .* ep_z;
[gxx, gyy, gzz] = meshgrid (xx, yy, zz);
gxx = gxx(:); gyy = gyy(:); gzz = gzz(:);
gch = convhull (gxx, gyy, gzz); % keep this (is common for all grid boxes)

for k=1:gp
  % move box to each grid point
  [x, y, z] = ind2sub(size(rgx), k);
  gxx_ = gxx + rgx(x,y,z);
  gyy_ = gyy + rgy(x,y,z);
  gzz_ = gzz + rgz(x,y,z);

  % vertices of convex hull of grid box
  X_ = [gxx_ gyy_ gzz_];
  
  gvertices0((k-1) * 12 + (1:12),:) = X_(gch(:,1),:);
  gvertices1((k-1) * 12 + (1:12),:) = X_(gch(:,2),:);
  gvertices2((k-1) * 12 + (1:12),:) = X_(gch(:,3),:);
  
 % trisurf (gch, gxx_, gyy_, gzz_, 'FaceAlpha', 0.3); % plot box
end

%% Set up layer boxes (not including waterlayer, but bottom layer beneath last interface)
nlayers = ninterfaces;
layers = cell(nlayers,1);
min_r_layer = NaN;
max_r_layer = NaN;

% create artifical bottom interface
bottom = interface_tris{end}.X;
bottom(:,3) = bounds(3,2);
bottom = TriRep (delaunay(bottom(:,1:2)), bottom);

for k=1:ninterfaces
  if (k==ninterfaces)
    l = [interface_tris{k}.X; bottom.X;];
  else
    l   = [interface_tris{k}.X; interface_tris{k+1}.X;];
  end
  dt  = convhull(l);
  L   = TriRep (dt, l);
  layers{k} = L;
  
  % figure out if this is maximum layer
  % if any part of this layer is inside R_tri maximum layer is
  % equal or greater
  inhull = any(~isnan(tsearchn(R_tri.X, delaunay(R_tri.X), L.X))) || ...
           any(~isnan(tsearchn(L.X, delaunay(L.X), R_tri.X)));
  if (inhull)
    max_r_layer = k +1;
    
    % first layer which is at least partly inside R_tri is min_r_layer
    if isnan(min_r_layer)
      min_r_layer = k + 1;
    end
  end

end

assert (size(velp,1) == size(vels,1), 'sg: velp and vels are not the same length');
assert (size(velp,1) >= max_r_layer, 'sg: velocity not specified for all layers in target grid');

%fprintf ('sg: min r layer: %d, max r layer: %d\n', min_r_layer, max_r_layer);

%% Ray points (phase dependandt)
points_offset = 0;
if strcmp(phase, 'P')
  % sea floor, layer
  ray_points = max_r_layer;
elseif strcmp(phase, 'SP')
  ray_points = max_r_layer;
elseif strcmp(phase, 'PIPBP')
  % sea floor, surface, sea floor, through layers..
  ray_points = max_r_layer + 2;
  points_offset = 2;
elseif (strcmp(phase, 'PIPBPIPBP'))
  % sea floor, surface, sea floor, surface, sea floor, through layers
  ray_points = max_r_layer + 4;
  points_offset = 4;
else
  error ('sg: unknown phase.');
end

%% Ray result vector: Z direction is coordinate of interface hit
cpoint = 1; % current layer
rays = nan(size(aa,1), 13, ray_points+1);
rays(1:end,1:3, cpoint) = [ aa(:,1), aa(:,2), aa(:,3) ]; % take off vector (normalized)
rays(:,4, cpoint) = stations(1);                    % take off coordinates
rays(:,5, cpoint) = stations(2);
rays(:,6, cpoint) = stations(3);
rays(:,7, cpoint) = 0;                              % length since last intersection
rays(:,8, cpoint) = 0;                              % type: 0 = original shot, 1 = transmitted, 2 = reflected
rays(:,9, cpoint) = 0;                              % termination:
                                                    %   0 = unsolved
                                                    %   1 = hit interface
                                                    %   2 = hit box (went out of model)
                                                    %   3 = hit surface
                                                    %   4 = error (ie bad direction vector)
rays(:,10,cpoint) = 0;                              % time elapsed (total)
rays(:,11,cpoint) = false;                          % matched to grid point
rays(:,12:13,cpoint) = [thetas phis];               % take off angle

clear tt pp tx ty tz;
drawnow;

r_ray_p = nan(rnx, rny, rnz, 3); % ray matching grid point
                                 % [ 
                                 %   ray number of rays(..)
                                 %   distance to point from latest intersection
                                 %   time of point intersection
                                 % ]
%%

%% Solve rays untill all grid points have been solved for this phase
allgridpoints_solved  = false; gridpoint_iter = 1; 
gridpoint_maxiter = 3; % max iterations for grid point

% determine 
while (~allgridpoints_solved)
  %fprintf ('sg: grid iteration: %i\n', gridpoint_iter);
  grt0 = tic;
  
  cpoint = 1;
  cinter = 1;

  while cpoint <= ray_points
    plt0 = tic;
    fprintf ('sg: solving for phase: %s, iteration %i, point: %i', phase, gridpoint_iter, cpoint);
    drawnow('update'); % flush output
    
    if (strcmp(phase, 'PIPBP')) % 1st multiple
      if (cpoint == 1)
        ccinter = 1;    % sea floor
      elseif (cpoint == 2)
        ccinter = NaN;  % surface
        hasinter = true;
        vert0 = svert0;
        vert1 = svert1;
        vert2 = svert2;
        nvert = snvert;
      elseif (cpoint == 3)
        ccinter = 1;    % sea floor
      else
        ccinter = cpoint - 2; % generic
      end
    elseif (strcmp(phase, 'PIPBPIPBP'))
      if (cpoint == 1)
        ccinter = 1; % sea floor
      elseif (cpoint == 2)
        ccinter = NaN;  % surface
        hasinter = true;
        vert0 = svert0;
        vert1 = svert1;
        vert2 = svert2;
        nvert = snvert;
      elseif (cpoint == 3)
        ccinter = 1; % sea floor
      elseif (cpoint == 4)
        ccinter = NaN;  % surface
        hasinter = true;
        vert0 = svert0;
        vert1 = svert1;
        vert2 = svert2;
        nvert = snvert;
      elseif (cpoint == 5)
        ccinter = 1; % sea floor
      else
        ccinter = cpoint - 4; % generic
      end
    elseif (strcmp(phase, 'SP')) % S converted to P
      ccinter = cpoint; % generic
    elseif (strcmp(phase, 'P'))
      ccinter = cpoint; % generic
    end
    
    % generic interface
    if (~isnan(ccinter))
      if (ccinter > ninterfaces)
        hasinter = false;
      else
        intf = interface_tris{ccinter};
        vert0 = intf.X(intf.Triangulation(:,1),:);
        vert1 = intf.X(intf.Triangulation(:,2),:);
        vert2 = intf.X(intf.Triangulation(:,3),:);
        nvert = size(vert0, 1);
        hasinter = true;
      end
    end
    
    if (cpoint == 1),
      rayc = 'b';
    elseif (cpoint == 2),
      rayc = 'g';
    elseif (cpoint == 3),
      rayc = 'y';
    elseif (cpoint == 4),
      rayc = 'c';
    elseif (cpoint == 5),
      rayc = 'm';
    else
      rayc = 'k';
    end

    % solve all orignal shots or transmitted rays for next layer
    % todo: this can be vectorized by repeating the intersections for each
    %       ray
    ri = find ((rays(:,9,cpoint) == 0));

    fprintf (', rays: %i..', length(ri));
    drawnow('update'); % flush output

    for k=ri'
      p = rays(k, 1:3, cpoint); % direction vector
      x = rays(k, 4:6, cpoint); % p0

      if (norm(p)>(1+10*eps) || (norm (p) < (1 - 10*eps)) || any(isnan(p))),
        erreps = ((norm(p)-1)/eps);
        %fprintf ('sg: direction vector is not normalized or is NaN, relative eps: %f, continuing..\n', erreps);
        rays (k, 9, cpoint)     = 4; % set ray to terminated and continue
        rays (k, 9, cpoint + 1) = 4;
        
        continue;
      end
      % find intersecting tetrahedron
      if (hasinter)
        tic;
        [intersect d u vv] = TriangleRayIntersection (repmat(x, nvert, 1), repmat(p, nvert, 1), vert0, vert1, ...
                                                      vert2, tri_opt);

        % Only keep rays with positive time (discard backwards or p0 intersects)
        ii        = (d > eps);
        ii_i      = find(ii==1);     % indexes
        intersect = intersect(ii);
        d         = d(ii);
        u         = u(ii);
        vv        = vv(ii);

      else
        intersect = 0;
      end

      termination = NaN;

      if any(intersect),
  %       fprintf('sg: number of: faces=%i, points=%i, intersections=%i; time=%f sec\n', ...
  %               size(intf.X,1), nvert, sum(intersect), toc);
  
        if (strcmp(phase, 'PIPBP'))
          if (cpoint == 2)
            termination = 3; % hit surface
          else
            termination = 1; % other interface
          end
        elseif (strcmp(phase, 'PIPBPIPBP'))
          if (cpoint == 2 || cpoint == 4)
            termination = 3; % hit surface
          else
            termination = 1; % other interface
          end
        else
          termination = 1; % hit interface
        end
      else
        % If no intersect at interface, find intersect at
        % model boundary
  %       fprintf ('sg: ray did not hit interface, tracing to box boundary.\n');
        tic;
        [intersect d u vv] = TriangleRayIntersection (repmat(x, bnvert, 1), repmat(p, bnvert, 1), bvert0, ...
                                                      bvert1, bvert2, tri_opt);

        % Only keep rays with positive time (discard backwards or p0 intersects)
        ii        = (d > eps);
        ii_i      = find(ii==1);   % indexes
        intersect = intersect(ii);
        d         = d(ii);
        u         = u(ii);
        vv        = vv(ii);

  %       fprintf('sg: box: number of: faces=%i, points=%i, intersections=%i; time=%f sec\n', ...
  %               size(intf.X,1), nvert, sum(intersect), toc);

        termination = 2; % hit box
      end



      if (~any(intersect))
        v = [p; x];
        l = [0 1; 200 1] * v;
        if (doplot)
          plot3(l(:,1), l(:,2), l(:,3), 'r');
        end

%         fprintf ('sg: error: No intersections could be found: not even box, continuing..\n');
        rays (k, 9, cpoint)     = 4; % set ray to terminated and continue
        rays (k, 9, cpoint + 1) = 4;
        
        continue;
      end

      % if more than one intersect is found the intersection is exactly at the
      % border between two tetrahedra of the interface.

      j = find(intersect == 1);
      if (length (j) ~= 1),
        %fprintf ('sg: error: more than one intersect solutions for ray, could have hit exactly between two tetrahedra, using first.\n')
        j = j(1);
      end

  %     fprintf ('sg: intersect, t = %f, u = %f\n', t(j), u(j));
      v = [p; x];
      l = [0 1; d(j) 1] * v;

      % only plot first shot
      %if (cpoint == 1)
      if (doplot)
        plot3 (l(:,1), l(:,2), l(:,3), rayc);
      end
      %end

      rays (k, 4:6, cpoint + 1) = l(2,:); % save intersection point in ray
      rays (k, 7, cpoint + 1)   = d(j);   % length of this segment

      rays (k, 9, cpoint)   = termination; % termination reason
      rays (k, 9, cpoint + 1) = termination; % will be changed below if continue trace
      rays (k, 11, cpoint + 1) = false;    % not matched
      if (strcmp(phase, 'P'))
        rays (k, 10, cpoint + 1)  = rays(k,10,cpoint) + d(j) / velp(cpoint); % time
      elseif (strcmp(phase, 'SP'))
        if (cpoint > 1)
          rays (k, 10, cpoint + 1)  = rays(k,10,cpoint) + d(j) / vels(cpoint); % time
        else
          rays (k, 10, cpoint + 1)  = rays(k,10,cpoint) + d(j) / velp(cpoint); % time
        end
      elseif (strcmp(phase, 'PIPBP'))
        if (cpoint <= 2)
          vcpoint = 1;
        else
          vcpoint = cpoint -2;
        end
        rays (k, 10, cpoint + 1)  = rays(k,10,cpoint) + d(j) / velp(vcpoint); % time
      elseif (strcmp(phase, 'PIPBPIPBP'))
        if (cpoint <= 4)
          vcpoint = 1;
        else
          vcpoint = cpoint -4;
        end
        rays (k, 10, cpoint + 1)  = rays(k,10,cpoint) + d(j) / velp(vcpoint); % time
      end
          

      if (termination == 1 || termination == 3),
        transmit = false;
        reflect  = false;
        
        if (strcmp(phase, 'PIPBP')),
          if (cpoint <= 2)
            reflect = true;
          else
            if (termination == 1)
              transmit = true;
            end
          end
        elseif (strcmp(phase, 'PIPBPIPBP'))
          if (cpoint <= 4)
            reflect = true;
          else
            if (termination == 1)
              transmit = true;
            end
          end
        else
          if (termination == 1)
            transmit = true;
          end
        end
        
        if (reflect),
          rays (k, 8, cpoint + 1) = 2; % reflect ray
          rays (k, 9, cpoint + 1) = 0; % unsolved

          vertex_i = ii_i(j);

          % calculate normal to triangle
          v0 = vert1(vertex_i,:)' - vert0(vertex_i,:)';
          v1 = vert2(vertex_i,:)' - vert0(vertex_i,:)';

          n = cross(v1, v0);
          n = n ./ norm(n);
          
          % calculate reflected ray: Snells law, Keers et. al. 97, p. 18, B9
          % http://en.wikipedia.org/wiki/Snells_law#Vector_form
          
          pn = 2 * dot(-p, n) * n + p';
          rays (k, 1:3, cpoint + 1) = pn';
        end
        if (transmit)
          rays (k, 8, cpoint + 1) = 1; % transmit ray
          rays (k, 9, cpoint + 1) = 0; % unsolved

          vertex_i = ii_i(j);

          % calculate normal to triangle
          v0 = vert1(vertex_i,:)' - vert0(vertex_i,:)';
          v1 = vert2(vertex_i,:)' - vert0(vertex_i,:)';

  %       plot3(vert0(vertex_i,1), vert0(vertex_i,2), vert0(vertex_i,3), 'kx');
  %       plot3(vert1(vertex_i,1), vert1(vertex_i,2), vert1(vertex_i,3), 'bx');
  %       plot3(vert2(vertex_i,1), vert2(vertex_i,2), vert2(vertex_i,3), 'gx');

          n = cross(v1, v0);
          n = n ./ norm(n);

          % calculate transmitted ray: Snells law, Keers et. al. 97, p. 18, B9
          % http://en.wikipedia.org/wiki/Snells_law#Vector_form
          if (strcmp(phase, 'P'))
            n1    = 1./velp(cpoint);
            n2    = 1./velp(cpoint + 1);
          elseif (strcmp(phase, 'SP'))
            if (cpoint == 1)
              n1    = 1./velp(cpoint);
              n2    = 1./vels(cpoint + 1);
            else
              n1    = 1./vels(cpoint);
              n2    = 1./vels(cpoint + 1);            
            end
          elseif (strcmp(phase, 'PIPBP'))
            % We get to seafloor at third point
            n1 = 1./velp(cpoint-2);
            n2 = 1./velp(cpoint-1);
          elseif (strcmp(phase, 'PIPBPIPBP'))
            % We got to seafloor at fifth point
            n1 = 1./velp(cpoint-4);
            n2 = 1./velp(cpoint-3);            
          end
          
%           cost1 = dot (n, p);
%           radicand = 1 - (n1/n2)^2 * (1 - cost1^2);
%           assert (radicand > 0, 'sg: error: negative radicand: total internal reflection');
%           cost2 = sqrt ( radicand );
%           pn    = (n1/n2) * p' + sign(cost1) * (- n1/n2 * cost1 + cost2) * n;

          % From Degreve, 2006:
          N     = n1 / n2;
          assert (N > 0, 'sg: error: currently no scenario with rays going from higher to lower velocities.');
          
          cosI  = -dot(n, p);
          sinT2 = N.^2 * (1.0 - cosI.^2);
          
          %assert (sinT2 < 1.0, 'sg: error: TIR.');
          
          if (sinT2 >= 1.0)
            % Total internal refraction
            %fprintf ('sg: TIR\n');
            rays (k, 9, cpoint)     = 4; % set ray to terminated and continue
            rays (k, 9, cpoint + 1) = 4;
          else
     
            cosT  = sqrt (1.0 - sinT2);

            pn    = N * p + (N * cosI - cosT) * n'; % p is normalized so this should also be normalized
            
            %assert (pn(3) <= p(3), 'sg: error: ray should not be refracted to a smaller angle in this scenario.');
            
            %assert (cost1 > 0, 'sg: error cost1 negative');
            rays (k, 1:3, cpoint + 1) = pn;
          end
          %quiver3 (l(2,1), l(2,2), l(2,3), n(1), n(2), n(3), 'b');
          %quiver3 (l(2,1), l(2,2), l(2,3), pn(1), pn(2), pn(3), 'r');
        end
      end

    end  

    % solve rays for next layer
    cpoint = cpoint + 1;
    cinter = cinter + 1;

    fprintf ('done: %f secs\n', toc(plt0));
    drawnow('update'); % flush output
  end
  fprintf ('sg: solved all rays, rays: %i, time: %f sec\n', size(rays,1), toc(t0));
  drawnow;

  %% match rays to grid points (each layer separately since the ray path changes for each layer)
  mtt0 = tic;
  gridmatches = 0;
  for kl=min_r_layer:max_r_layer
    % 1. find rays that hit convex hull of r (search rays in region layer)
    % 2. figure out which hit any grid point boxes and match
    
    % TODO: Don't check for intersects to grid cells not in current layer
    
    R_layer = kl + points_offset;

    layer_X   = layers{kl-1}.X;
    layer_dt  = delaunay (layer_X);
    
    mt0 = tic;

    % solved rays in region layer that have not been matched before
    rri = find((rays(:,9,R_layer-1) == 1) & ...
               (rays(:,9,R_layer) > 0)    & ...
               (rays(:,11,R_layer) == false));


    rvert0 = R_tri.X(R_tri.Triangulation(:,1),:);
    rvert1 = R_tri.X(R_tri.Triangulation(:,2),:);
    rvert2 = R_tri.X(R_tri.Triangulation(:,3),:);
    rnvert = size(rvert0, 1);

    fprintf ('sg: matching rays to grid points.. (layer: %i, rays: %i)..', kl, length(rri));
    drawnow('update'); % flush output
    for rk=rri'
      % check if ray intersects convex hull of R
      p = rays(rk, 1:3, R_layer); % direction vector
      x = rays(rk, 4:6, R_layer); % p0

      % mark ray matched
      rays(rk,11,R_layer) = true;

      v = [p; x];
      l = [0 1; rays(rk,7,R_layer+1) 1] * v;

  %     if (gridpoint_iter > 1)
  %       plot3 (l(:,1), l(:,2), l(:,3), 'y'); drawnow;
  %     end

      if (norm(p)>(1+10*eps) || (norm (p) < (1 - 10*eps)) || any(isnan(p))),
          erreps = ((norm(p)-1)/eps);
          %fprintf ('sg: direction vector is not normalized or is NaN, relative eps: %f, continuing..\n', erreps);

          continue;
        end

      % find intersecting tetrahedron

      it0 = tic;
      [intersect d u vv] = TriangleRayIntersection (repmat(x, rnvert, 1), repmat(p, rnvert, 1), rvert0, rvert1, ...
                                                    rvert2, tri_opt);

      % Only keep rays with positive time (discard backwards or p0 intersects)
      ii        = (d > eps);
      ii_i      = find(ii==1);     % indexes
      intersect = intersect(ii);
      d         = d(ii);
      u         = u(ii);
      vv        = vv(ii);


      %% Match intersecting ray to grid boxes
      if any(intersect)

        % test against inside boxes
  %       fprintf ('sg: found match (%f secs), determining grid point..\n', toc(it0));
        bt0 = tic;
        [intersect d u vv] = TriangleRayIntersection (repmat(x, ngvertices, 1), repmat(p, ngvertices, 1), gvertices0, gvertices1, ...
                                                      gvertices2, tri_opt);

  %       fprintf ('sg: found %i matching triangles (%f secs), determining grid point\n', length(intersect(intersect==1)), toc(bt0));

        ii_i = find((intersect==1) & (d>eps));
        gpi  = floor (ii_i ./ 12) +1;
        [~, ui] = unique (gpi);
        uii = ones(size(gpi));
        uii(ui) = 0; % mask points without two intersects
        gpi_u = gpi(uii==1);
        %gpi_u = unique (gpi);

        %fprintf ('sg: %i grid points matched\n', length(gpi_u));

        % plot matching line
        if (doplot)
          plot3 (l(:,1), l(:,2), l(:,3), 'r');
        end
        
        % mask out grid point matches outside current layer
        if (~isempty(gpi_u))
          inhull = ~isnan(tsearchn(layer_X, layer_dt, [rgx(gpi_u) rgy(gpi_u) rgz(gpi_u)]));
          gpi_u = gpi_u(inhull);
        end

        % pick match if better than any previous match and plot boxes of all solutions
        for b=gpi_u'

          % index
          [xi, yi, zi] = ind2sub(size(rgx), b);

          gxx_ = gxx + rgx(xi,yi,zi);
          gyy_ = gyy + rgy(xi,yi,zi);
          gzz_ = gzz + rgz(xi,yi,zi);

          if (doplot)
            trisurf (gch, gxx_, gyy_, gzz_, 'FaceColor', 'b', 'FaceAlpha', 0.7);
          end

          gridmatches = gridmatches + 1;


          % NOT USED: calculate distance from intersection (closest point on ray to
          %                                       grid center)
          % http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
          gp0 = [rgx(xi,yi,zi) rgy(xi,yi,zi) rgz(xi,yi,zi)]; % center of grid point
          %gpl     = (x - gp0) - (dot((x - gp0), gp0)) .* p;

          gpd = norm(gp0 - x); % distance from intersection point (interface)
          dgp = gpd + sum(rays(rk, 7, 1:R_layer));

          % calculate time of intersection
          if (strcmp(phase, 'P'))
            dgt = gpd / velp(R_layer) + rays(rk, 10, R_layer);
          elseif (strcmp(phase, 'PIPBP'))
            dgt = gpd / velp(R_layer-points_offset) + rays(rk, 10, R_layer);
          elseif (strcmp(phase, 'PIPBPIPBP'))
            dgt = gpd / velp(R_layer-points_offset) + rays(rk, 10, R_layer);          
          elseif (strcmp(phase, 'SP'))
            dgt = gpd / vels(R_layer) + rays(rk, 10, R_layer);
          end

          % store grid point <-> ray solution 
          if (isnan(r_ray_p(xi, yi, zi, 3)) || ...
              (dgt > r_ray_p(xi, yi, zi, 3)))
            %fprintf ('sg: found new solution for grid point: %i, %i, %i\n', rbg(b,1), ...
            %         rbg(b,2), rbg(b,3));
            r_ray_p(xi, yi, zi, 1) = rk;  % ray no
            r_ray_p(xi, yi, zi, 2) = dgp; % distance
            r_ray_p(xi, yi, zi, 3) = dgt; % time              
          end
        end
      end
    end
    
    fprintf ('done: %f secs.\n', toc(mt0));
  end

  fprintf ('sg: matching done: %i matches, %f secs, ', gridmatches, toc(mtt0));

  % check if all grid points have been solved
  unsolved_gridpoints = isnan(r_ray_p(:,:,:,1));
  fprintf ('unsolved gridpoints: %i\n', length(unsolved_gridpoints(unsolved_gridpoints==1)));
  drawnow('update'); % flush output
  
  if (~isempty(unsolved_gridpoints))  
    %% Interpolate new rays for unsolved grid points
    cpoint = 1;
    cinter = 1;
    
    % grid points to find
    urb  = [rgx(unsolved_gridpoints) rgy(unsolved_gridpoints) rgz(unsolved_gridpoints)];

    solved_gridpoints = ~unsolved_gridpoints;
    r_ray_p_t = r_ray_p(:,:,:,1);
    solvedrays = r_ray_p_t(solved_gridpoints);
    %solvedrays = find(rays(:,9,R_layer+1)==2 | ...
    %                  rays(:,9,R_layer+1)==1); % rays that have been solved in
                                               % layer of region of interest
    fth  = rays(solvedrays,12,1); % take off angle
    fphi = rays(solvedrays,13,1); % take off angle
    %gr   = rays(solvedrays,4:6,R_layer+1);
    gr = [rgx(solved_gridpoints) rgy(solved_gridpoints) rgz(solved_gridpoints)];
    
    % TriScatterInterp only interpolates within the convex hull of its
    % function
        
    Fth  = TriScatteredInterp (gr, fth);
    Fphi = TriScatteredInterp (gr, fphi);

    thetas  = Fth(urb);
    phis    = Fphi(urb);

%     thetas = griddata3(gr(:,1), gr(:,2), gr(:,3), fth, urb(:,1), urb(:,2), urb(:,3));
%     phis   = griddata3(gr(:,1), gr(:,2), gr(:,3), fphi, urb(:,1), urb(:,2), urb(:,3));

    
    % remove NaNs
    inan = ~isnan(thetas);
    if (any(inan))
      fprintf ('sg: warning: nan output from new angle interpolation.\n');
    end

    thetas = thetas(inan);
    phis   = phis(inan);

    % remove duplicates from new interpolated rays
%     [~, tii, ~] = unique (thetas);
%     [~, pii, ~] = unique (phis);
%     uii = zeros(size(thetas));
%     uii(tii)  = 1;
%     uii(pii)  = 1;
%     uii       = logical(uii); % make logical array (not index)
%     thetas    = thetas(uii);
%     phis      = phis(uii);
    
    targ = [thetas phis];
    [~, uii, ~] = unique (targ, 'rows');
    thetas    = thetas(uii);
    phis      = phis(uii);

    % remove duplicates in new rays and old rays
    tii = ismember (thetas, rays(:,12,1));
    pii = ismember (phis, rays(:,13,1));
    uii = (~pii | ~tii);
    thetas = thetas(uii);
    phis   = phis(uii);    

    if (~isempty(thetas))
      tx = cos(phis) .* sin(thetas);
      ty = sin(phis) .* sin(thetas);
      tz = cos(thetas);

      nrays = nan(length(tx(:)), 13, ray_points+1);
      nrays(1:end,1:3, cpoint) = [ tx(:), ty(:), tz(:) ];  % take off vector (normalized)
      nrays(:,4, cpoint) = stations(1);              % take off coordinates
      nrays(:,5, cpoint) = stations(2);
      nrays(:,6, cpoint) = stations(3);
      nrays(:,7, cpoint) = 0;                              % length since last intersection
      nrays(:,8, cpoint) = 0;                              % type: 0 = original shot, 1 = transmitted, 2 = reflected
      nrays(:,9, cpoint) = 0;                              % termination:
                                                           %   0 = unsolved
                                                           %   1 = hit interface
                                                           %   2 = hit box (went out of model)
                                                           %   3 = hit surface
                                                           %   4 = bad termination
      nrays(:,10,cpoint) = 0;                              % time elapsed (total)
      nrays(:,11,cpoint) = false;                          % matched to grid point
      nrays(:,12:13,cpoint) = [thetas phis];               % take off angle

      rays = cat(1, rays, nrays);
    else
      fprintf ('sg: no new rays could be calculated, giving up.\n');
      allgridpoints_solved = true;
    end
    %%
  else
    fprintf ('sg: all grid points solved.\n');
    allgridpoints_solved = true;
  end

  if (gridpoint_iter >= gridpoint_maxiter)
    allgridpoints_solved = true;
  end

  gridpoint_iter = gridpoint_iter + 1;
  fprintf ('sg: grid iteration done: %f secs\n', toc(grt0));
end

%%

fprintf ('sg: done: %f secs\n', toc(t0));

end
