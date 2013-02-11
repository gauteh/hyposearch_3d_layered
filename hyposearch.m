function [traveltimes] = hyposearch (bounds, interface_tris, stations, velp, vels, R_tri, rnx, rny, rnz)
  % Hyposearch
  %
  %
  % Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-02-11
  %
  %

  fprintf ('hs: hyposearch\n');
  
  ht0 = tic;

  %% Set up model
  fprintf ('hs: setting up models..\n');
  
  ninterfaces = size(interface_tris,1); % interfaces in model
  nlayers     = ninterfaces + 1;        % and layers (between interfaces)

  % figure(2); clf('reset');
  % l = 0:ninterfaces;
  % [haxes, hl1, hl2] = plotyy (l, [velp vels], l(2:end), tcoeff, @stairs);
  % set (haxes(:), 'XTick', l);
  % ylabel (haxes(1), 'Velocity [m/s]')
  % ylabel (haxes(2), 'Transmission coeff')
  % set (hl1, 'Marker', 'x');
  % set (hl2, 'Marker', 'x');
  % xlabel ('Layer');
  % title ('Velocity model parameters');
  % legend ('P velocity', 'S velocity', 'Transmission coefficient');

  %% Interfaces

  % Interface:
  %  - Function z of x and y defed along arbitrary points
  %  - Transmission and reflection coefficients
  %
  % Must be defined on entire model within boundaries (var: bounds)
  %interface_tris = cell(ninterfaces, 1);


  % Find Delaunay triangulation for interface
  %fprintf ('hs: setting up interfaces..\n');
  %interface_tris{1} = interfaces;

  % Lower velocity increase (dont track reflections)
  % I = [ 0     0    6000;
  %       0   2000    6000;
  %       2500   0    6000;
  %       2500 2000    6000;
  %      ];
  % interfaces{2} = I;

  %%
  
  %% Solve all phases for all stations
  nstations = size(stations,1);
  sst0 = tic;
  
  % Cell array to store traveltime results
  traveltimes = cell (nstations, 8);
  rays    = cell(nstations, 4);
  matches = cell(nstations, 4);
  
  % grid points of region of interest
  %[~, ~, ~, rb, ~, ~, ~, ~, ~, ~, ~] = makegrid (R_tri);
  
  %solvephases = { 'P', 'SP', 'PIPBP', 'PIPBPIPBP' };
  njobs = nstations * 4;
  
  parfor ns=1:njobs
    [ks, ps] = ind2sub([nstations 4], ns);
    
    st0 = tic;
    
    station = stations(ks, :);
    
    phf = ps; % figure
    ph  = 0;
    
    if (ps == 1)
      ph = 'P';
    elseif (ps == 2)
      ph = 'SP';
    elseif (ps == 3)
      ph = 'PIPBP';
    elseif (ps == 4)
      ph = 'PIPBPIPBP';
    end
    
    fprintf ('hs: solving for station %i, phase: %s..\n', stations(ks,3), ph);
    
    [r, r_r] = solvegrid (station, bounds, velp, vels, interface_tris, R_tri, ...
                                   rnx, rny, rnz, ninterfaces, nlayers, ph, phf);
    
    rays(ns)    = { r };
    matches(ns) = { r_r };
    %% Solve for grid in region of interest for phase P
%    [rays_p, r_ray_p] = solvegrid (station, bounds, velp, vels, interface_tris, R_tri, ...
 %                                  R_layer, rnx, rny, rnz, ninterfaces, nlayers, 'P', 1);

    
    %% Mask out unused rays and re-index solved grid points
%     [gnx, gny, gnz, ~] = size(r_ray_p);
%     matched_rays    = r_ray_p(:,:,:,1);
%     solved          = ~isnan(matched_rays);
%     f_solved        = solved(:);
%     f_matched_rays  = matched_rays(solved); % flatten
% 
%     rays_p = rays_p(f_matched_rays,:,:); % mask out unused rays
% 
%     % Re-index matched rays in gridpoint ray match
%     f_matched_rays = nan(size(f_solved));
%     f_matched_rays(f_solved) = 1:size(rays_p,1);
%     matched_rays = reshape(f_matched_rays, gnx, gny, gnz);
% 
%     r_ray_p(:,:,:,1) = matched_rays;

    % Store
    %traveltimes(ks, :) = {rays_p, r_ray_p};
%     traveltimes(ks, 2) = {r_ray_p};
    
    %% Solve for grid in region of interest for phase SP (converted at seafloor)
%    [rays_s, r_ray_s] = solvegrid (station, bounds, velp, vels, interface_tris, R_tri, ...
 %                                  R_layer, rnx, rny, rnz, ninterfaces, nlayers, 'SP', 2);


    %% Mask out unused rays and re-index solved grid points
%     [gnx, gny, gnz, ~] = size(r_ray_s);
%     matched_rays    = r_ray_s(:,:,:,1);
%     solved          = ~isnan(matched_rays);
%     f_solved        = solved(:);
%     f_matched_rays  = matched_rays(solved); % flatten
% 
%     rays_s = rays_s(f_matched_rays,:,:); % mask out unused rays
% 
%     % Re-index matched rays in gridpoint ray match
%     f_matched_rays = nan(size(f_solved));
%     f_matched_rays(f_solved) = 1:size(rays_s,1);
%     matched_rays = reshape(f_matched_rays, gnx, gny, gnz);
% 
%     r_ray_s(:,:,:,1) = matched_rays;
    
    % Store
    %traveltimes(ks, 3) = {rays_s};
    %traveltimes(ks, 4) = {r_ray_s};


  %   phases:
  %   P = transmitted P ray all the way to the ice, surface
  %   S = transmitted S ray to the sea floor, B
  %   SP = converted S to P at the sea floor
  %   B = bottom, sea floor reflection
  %   I = ice, surface reflection

    %% Solve for grid in region of interest for phase PIPBP(first multiple)
%    [rays_m, r_ray_m] = solvegrid (station, bounds, velp, vels, interface_tris, R_tri, ...
 %                                 R_layer, rnx, rny, rnz, ninterfaces, nlayers, 'PIPBP', 3);


    %% Mask out unused rays and re-index solved grid points
%     [gnx, gny, gnz, ~] = size(r_ray_m);
%     matched_rays    = r_ray_m(:,:,:,1);
%     solved          = ~isnan(matched_rays);
%     f_solved        = solved(:);
%     f_matched_rays  = matched_rays(solved); % flatten
% 
%     rays_m = rays_m(f_matched_rays,:,:); % mask out unused rays
% 
%     % Re-index matched rays in gridpoint ray match
%     f_matched_rays = nan(size(f_solved));
%     f_matched_rays(f_solved) = 1:size(rays_m,1);
%     matched_rays = reshape(f_matched_rays, gnx, gny, gnz);
% 
%     r_ray_m(:,:,:,1) = matched_rays;
    
    % Store
    %traveltimes(ks, 5) = {rays_m};
    %traveltimes(ks, 6) = {r_ray_m};


    %% Solve for grid in region of interest for phase PIPBPIPBP (second multiple)
%    [rays_mm, r_ray_mm] = solvegrid (station, bounds, velp, vels, interface_tris, R_tri, ...
 %                                  R_layer, rnx, rny, rnz, ninterfaces, nlayers, 'PIPBPIPBP', 4);


    %% Mask out unused rays and re-index solved grid points
%     [gnx, gny, gnz, ~] = size(r_ray_mm);
%     matched_rays    = r_ray_mm(:,:,:,1);
%     solved          = ~isnan(matched_rays);
%     f_solved        = solved(:);
%     f_matched_rays  = matched_rays(solved); % flatten
% 
%     rays_mm = rays_mm(f_matched_rays,:,:); % mask out unused rays
% 
%     % Re-index matched rays in gridpoint ray match
%     f_matched_rays = nan(size(f_solved));
%     f_matched_rays(f_solved) = 1:size(rays_mm,1);
%     matched_rays = reshape(f_matched_rays, gnx, gny, gnz);
% 
%     r_ray_mm(:,:,:,1) = matched_rays;
    
    % Store
    %traveltimes(ks, 7) = {rays_mm};
    %traveltimes(ks, 8) = {r_ray_mm};
    
    %traveltimes(ks,:) = {rays_p, r_ray_p, rays_s, r_ray_s, rays_m, r_ray_m, rays_mm, r_ray_mm};
    %traveltimes(ns:(ns+1)) = rays_p;
    %traveltimes(ns + 1) = r_ray_p;
    fprintf ('hs: station %i, phase %s done: %f secs\n', stations(ks,3), ph, toc(st0));
  end
  
  fprintf ('hs: all stations done: %f secs.\n', toc(sst0));
  
  % Combine to traveltimes
  traveltimes(:,[1 3 5 7]) = rays;
  traveltimes(:,[2 4 6 8]) = matches;
  
  fprintf ('hs: done: %f secs.\n', toc(ht0));
end                           

