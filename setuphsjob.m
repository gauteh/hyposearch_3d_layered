function [bounds, interface_tris, stations, quakes, phases, velp, vels, R_tri, usephases, event, jobname, rnx, rny, rnz] = setuphsjob (defaults)
% setuphsjob: load job from current directory and prepare for solving
%
% author: Gaute Hope <eg@gaute.vetsj.com> / 2013-02-18
%
% defaults specify whether to use default settings or to load settings from
% file hs_job.m

  fprintf ('hs: setup: loading batyhymetry and job data..\n');
  
  if (~defaults)
    fprintf ('hs: loading hs_job.m..');
    hs_job;
  else
    fprintf ('hs: loading defaults..');
    margins       = 10000 .* ones(1,4);
    seafloortype  = 1;
    solvephases   = ones(4,4);
    usephases     = ones(4,4);
    velp          = [1500; 5800];
    vels          = [   0; 3200];
    
    event         = 'Default event';
    jobname       = 'Default job';
  end
  
  fprintf ('event: %s, job: %s\n', event, jobname);
  
  figure(1); clf('reset');
  
  if (seafloortype == 1)
    % load ibcao of this region
    ncid = netcdf.open ('ibcao_det.grd', 'NC_NOWRITE');
    x = netcdf.getVar(ncid, 0); % longitude
    y = netcdf.getVar(ncid, 1); % latitude
    z = netcdf.getVar(ncid, 2); % bathymetry (depth is negative)
    z = double(z);

    dx = 500; dy = 500; % meters

    netcdf.close (ncid);

    subplot(3,1, 1);
    %surf (y, x, z);
    imagesc(y, x, z); colormap (jet);
    title ('IBCAO (general extracted region)');
    hold all;

    % load region
    region = load('ibcao_det.coor');

  end
  
  % load stations
  stations = load('stations.coor');
  
  if (seafloortype == 1 || seafloortype == 2)
    cc = ['g', 'y', 'k'];
    h = [];
    for i=1:3
      h(i) = scatter(stations(i,2), stations(i,1),  sprintf('%co', cc(stations(i,3)+1)), 'filled');
    end

    legend (h', 'GAK2', 'GAK3', 'GAK4');
  end
  
  % load quake
  quakes = load ('quakes.coor');
  cc = ['r', 'b'];
  
  if (seafloortype == 1 || seafloortype == 2)
    h = [];
    for i=1:size(quakes,1)
      h(i) = scatter(quakes(i,2), quakes(i,1), sprintf('%co', cc(quakes(i,3)+1)), 'filled');
    end
  end

  %legend (h, 'HYPOSAT solutions');

  % margins
  nodes = cat(1, stations, quakes);
  if (~region_auto)
    nodes = cat(1, nodes, [rxmin rymin rzmin]);
    nodes = cat(1, nodes, [rxmax rymax rzmax]);
  end
  
  if (interfaces_add)
    for k=1:(ninterfaces-1)
      nodes = cat(1, nodes, interface_tris_add{k}.X);
    end
  end
  
  l = min(nodes(:,1));
  r = max(nodes(:,1));
  u = max(nodes(:,2));
  d = min(nodes(:,2));

  if (region_auto)
    l = l - margins(3);
    r = r + margins(4);
    u = u + margins(1);
    d = d - margins(2);
  end
  % 
  % [yy, xx] = meshgrid(y, x);
  % 
  % iyy = (yy<u) & (yy>d);
  % ixx = (xx<r) & (xx>l);
  % 
  % ii  = ixx & iyy;

  if (seafloortype == 1)
    % constrain IBCAO
    xmin = find(x>l, 1, 'first');
    xmax = find(x<r, 1, 'last');
    ymin = find(y>d, 1, 'first');
    ymax = find(y<u, 1, 'last');

    clear l r u d m ncid;

    xx = x(xmin:xmax);
    yy = y(ymin:ymax);
    zz = z(xmin:xmax,ymin:ymax);

    subplot(3,1,2);
    %figure(7); clf('reset');
    imagesc(yy, xx, zz);
    colormap (jet);
    hold on;
    title ('IBCAO (constrained region for solving)');
    
  elseif (seafloortype == 0)
    % fixed depth
    xx = [l r];
    yy = [d u];
    zz = -seafloor_depth .* ones(2);
    
    subplot(2,1,1);
    %figure(7); clf('reset');
    imagesc(yy, xx, zz);
    colormap (jet);
    hold on;
    title ('Fixed depth');
  end
  
  % load stations
  cc = ['g', 'y', 'k'];
  h = [];
  for i=1:3
    h(i) = scatter(stations(i,2), stations(i,1),  sprintf('%co', cc(stations(i,3)+1)), 'filled');
  end

  legend (h, 'GAK2', 'GAK3', 'GAK4');

  % load quake
  cc = ['r', 'b'];
  h = [];
  for i=1:size(quakes,1)
    h(i) = scatter(quakes(i,2), quakes(i,1), sprintf('%co', cc(quakes(i,3)+1)), 'filled');
  end

  % Triangulate surface
  if (seafloortype == 2)
    bathymetry = seafloor;
  else
    [yg, xg] = meshgrid (yy, xx);
    yg = yg(:); xg = xg(:);
    dt = delaunay(xg, yg);

    bathymetry = TriRep (dt, xg, yg, -zz(:));
  end
  
  trisurf (bathymetry);
  title ('Bathymetry triangulation (depth positive, constrained region)');
  set(gca, 'ZDir', 'reverse');
  colorbar;

  bounds = [xx(1) xx(end); 
            yy(1) yy(end); 
            0     7000];
          
  
  deepest = bounds(3,2); % deepest so far


  %% Interfaces
  interface_tris    = cell(ninterfaces,1);
  interface_tris{1} = bathymetry;
  deepest = max([deepest bathymetry.X(:,3)']);
  fprintf ('hs: setup: interfaces: %d\n', ninterfaces);
  if (interfaces_add)
    for k=1:(ninterfaces-1)
      interface_tris{k+1} = interface_tris_add{k};
      trisurf (interface_tris_add{k}, 'FaceAlpha', 0.6);
      
      deepest = max([deepest interface_tris{k+1}.X(:,3)']);
    end
  end

  %% load phasereadings
  fid = fopen ('phases.tt');
  phasedata = textscan(fid, '%s%s%d%d%d%d%d%f');
  fclose(fid);

  station_n = phasedata{1};
  phase_n   = phasedata{2};

  year    = double(phasedata{3});
  month   = double(phasedata{4});
  day     = double(phasedata{5});
  hour    = double(phasedata{6});
  minute  = double(phasedata{7});
  second  = double(phasedata{8});


  % gak2 = 0, gak3 = 1, gak4 = 2, gaks = 3
  station = zeros(length(station_n), 1);
  station(strcmp(station_n, 'GAK2')) = 0;
  station(strcmp(station_n, 'GAK3')) = 1;
  station(strcmp(station_n, 'GAK4')) = 2;
  station(strcmp(station_n, 'GAKS')) = 3;

  % phase: P = 1, S = 2, M = 3, MM = 4
  phases = zeros(length(phase_n), 2);
  phases(strcmp(phase_n, 'P_'),2) = 1;
  phases(strcmp(phase_n, 'S_'),2) = 2;
  phases(strcmp(phase_n, 'M_'),2) = 3;
  phases(strcmp(phase_n, 'MM'),2) = 4;
  time    = datenum (year, month, day, hour, minute, second);
  time    = time .* 24 .* 60 .* 60; % change to seconds
  phases(:,3) = time;
  phases(:,1) = station;
  
  % check that we got all required phases
  nstations = size(stations,1);
  for ks=1:nstations
    s = stations(ks,3);
    assert (length(find( (phases(:,1)==s) & (phases(:,2) == 1) )) == 1, 'No or too many P phase for station: %d', s);
    assert (length(find( (phases(:,1)==s) & (phases(:,2) == 2) )) == 1, 'No or too many S phase for station: %d', s);
    assert (length(find( (phases(:,1)==s) & (phases(:,2) == 3) )) < 2, 'Too many M phase for station: %d', s);
    assert (length(find( (phases(:,1)==s) & (phases(:,2) == 4) )) < 2, 'Too many MM phase for station: %d', s);
    
    if (isempty(find( (phases(:,1)==s) & (phases(:,2) == 3), 1 )))
      fprintf ('hs: Empty phase M for station: %d.\n', s);
      phases = cat(1, phases, [s 3 NaN]);
    end
    
    if (isempty(find( (phases(:,1)==s) & (phases(:,2) == 4), 1 )))
      fprintf ('hs: Empty phase MM for station: %d.\n', s);
      phases = cat(1, phases, [s 4 NaN]);
    end
      
  end

  %% Region of interest box
  fprintf ('hs: setup: region of interest..\n');
  if (region_auto)
    rxmin = xx(1); rxmax = xx(end);
    rymin = yy(1); rymax = yy(end);
    rzmin = 4000; rzmax = 8000;
  else
    fprintf ('hs: setup: using manually defined bounds of region.\n');
  end
  
  deepest = max([deepest rzmax]);
  bounds(3,2) = deepest + 1000; % with some margin
  
  % TODO: This layer does not have a velocity, when R_tri goes into it it
  % will become a problem: should in that case be set in hs_job.
  
  if (region_size_auto)
    rnx = 10;
    rny = 10;
    rnz = 10;
  else
    fprintf ('hs: setup: using manually defined grid size.\n');
  end

  % Create box around region of interest
  [bxg, byg, bzg] = meshgrid ([rxmin rxmax], [rymin rymax], [rzmin rzmax]);
  bxg = bxg(:); byg = byg(:); bzg = bzg(:);

  ch  = convhull (bxg, byg, bzg);
  R_tri = TriRep (ch, bxg, byg, bzg);
  trisurf (R_tri, 'FaceAlpha', 0.7, 'FaceColor', 'r');
  clear bxg byg bzg box;


  %% Job ready
  fprintf ('hs: setup: job ready.\n');
end

