function [mindt_i, s, t_dterr, dterrs, tttime, oterrs, otime, t_otime, t_otimeerr] = searchgrid (R_tri, rnx, rny, rnz, stations, quakes, phasereadings, traveltimes, dt, fig, figtitle)

% searchgrid: Search and match grid <-> ray pairs and calculate rms value
%
% Author: Gaute Hope <eg@gaute.vetsj.com> / 2013-02-18
%
% input:
% 
%   R_tri:          Triangulation of region of interest
%   stations:       Stations
%   quakes:         Hyposat quakes
%   phasereadings:  Absolute time phase readings.
%   traveltimes:    Calculated traveltimes with rays and matches
%   ot:             Which phases to use
%   fig:            Which figures to use for plotting
%   figtitle:       Title of figure
%
% output:
%
%   i:              Best match index in rb
%   s:              Values of best match in rb
%   t_dterr:        Total rms
%   dterrs:         Station rms values

if (length(fig) < 2)
  error ('sr: need 2 figure handles\n');
end

figure(fig(1)); clf('reset');

fprintf ('sr: calculate rms of origin time and differences..\n');

%% Setup
nstations = size(stations,1);
nphases   = 4;

[rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);

dterrs    = zeros(rnx, rny, rnz, nstations); % phase differences error
oterrs    = zeros(rnx, rny, rnz, nstations); % origin time error
tttime    = zeros(rnx, rny, rnz, nstations, nphases); % theoretical travel times for all stations and all phases
otime     = zeros(rnx, rny, rnz, nstations); % origin time
t_otime   = zeros(rnx, rny, rnz);            % avarage origin time
t_otimeerr = zeros(rnx, rny, rnz);           % origin time error to avarage time

t_dterr   = zeros(rnx, rny, rnz); % total differences error on all stations
%t_ttderr  = zeros(rnx, rny, rnz);

origtime = zeros(nstations, nphases); % original arrival time (readings from input)

for ks=1:nstations
  pr = phasereadings(phasereadings(:,1)==stations(ks,3), 2:3);
  
  %% Phase P
  pt    = pr(pr(:,1)==1,2);
  ray_p = traveltimes{ks,1};
  rp    = traveltimes{ks,2};
  

  dtp = rp(:,:,:,3);  % travel time
  otp = pt - dtp;     % origin time
  origtime(ks, 1) = pt;
  %%

  %% Phase S
  st    = pr(pr(:,1)==2,2);
  ray_s = traveltimes{ks,3};
  rs    = traveltimes{ks,4};

  dts = rs(:,:,:,3);  % travel time
  ots = st - dts;     % origin time
  origtime(ks, 2) = st;
  %%

  %% Phase M
  mt    = pr(pr(:,1)==3,2);
  ray_m = traveltimes{ks,5};
  rm    = traveltimes{ks,6};

  dtm = rm(:,:,:,3);  % travel time
  otm = mt - dtm;     % origin time
  origtime(ks, 3) = mt;
  %%

  %% Phase MM
  mmt    = pr(pr(:,1)==4,2);
  ray_mm = traveltimes{ks,7};
  rmm    = traveltimes{ks,8};

  dtmm = rmm(:,:,:,3);  % travel time
  otmm = mmt - dtmm;    % origin time
  origtime(ks, 4) = mmt;
  %%

  %% Difference time error
  
  % unweighted
  dterr = zeros(size(otp));
  if (dt(1,2))
    dterr = ((st - pt) - (dts - dtp)).^2;
  end
  
  if (dt(1,3))
    dterr = dterr + ((mt - pt) - (dtm - dtp)).^2;
  end
  
  if (dt(1,4))
    dterr = dterr + ((mmt - pt) - (dtmm - dtp)).^2;
  end
  
  if (dt(2,3))
    dterr = dterr + ((mt - st) - (dtm - dts)).^2;
  end
  
  if (dt(2,4))
    dterr = dterr + ((mmt - st) - (dtmm - dts)).^2;
  end
  
  if (dt(3,4))
    dterr = dterr + ((mmt - mt) - (dtmm - dtm));
  end
  
  n = dt(1,2) + dt(1,3) + dt(1,4) + dt(2,3) + dt(2,4) + dt(3,4);

  %dterr   = sqrt(dterr ./ n);
  t_dterr = t_dterr + dterr;
  dterrs(:,:,:,ks) = dterr;
  
  %% Plot result
  az = 0; el = 90;
  
  subplot(nstations,1,ks);
%   slice (oterr, 1:5,1:5,1:5);
  S = [rgx(:) rgy(:) rgz(:) dterr(:)];
  
  Sn = S(isnan(S(:,4)),:);
  S  = S(~isnan(S(:,4)),:);
  hold on;
  
  % plot HYPOSAT quakes
  %scatter3(quakes(:,1), quakes(:,2), 6000 .* ones(size(quakes,1),1), 100, 'ob', 'filled');
  
  h = scatter3(S(:,1), S(:,2), S(:,3), 100, S(:,4), 'filled', 'MarkerEdgeColor', 'k');
  %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'MarkerEdgeColor', 'k');
  %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'rx');
  
  xlabel ('x [m] UPS');
  ylabel ('y [m] UPS');
  zlabel ('depth [m]');
  colormap (hot); colorbar;
  set(gca, 'ZDir', 'reverse');
  view (az, el);
  
  title (sprintf('%s: RMS: Differences, station: %i', figtitle, stations(ks,3)));

%   subplot (nstations,2,(ks-1)*2 +2);
% %   slice (ttderr, 3,2,3);
%   S = [rb ttderr(:)];
%   
%   Sn = S(isnan(S(:,4)),:);
%   S  = S(~isnan(S(:,4)),:);
%   hold on;
%   
%   % plot HYPOSAT quakes
%   scatter3(quakes(:,1), quakes(:,2), [6000 6000], 100, 'ob', 'filled');
%   
%   scatter3(S(:,1), S(:,2), S(:,3), 100, S(:,4), 'filled', 'MarkerEdgeColor', 'k');
%   %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'MarkerEdgeColor', 'k');
%   %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'rx');
%   
% 
%   
%   xlabel ('x [m] UPS');
%   ylabel ('y [m] UPS');
%   zlabel ('depth [m]');
%   colormap (hot); colorbar;
%   set(gca, 'ZDir', 'reverse');
%   view (az, el);
%   title (sprintf('RMS: Differences, station: %i', stations(ks,3)));

  [mindt, i_mindt] = min(dterr(:));
%   minttd  = min(min(min(ttderr)));
  fprintf ('sr: station %i, minimum rms, differences: %f.\n', stations(ks,3), mindt);
  
  %% calculate origin time
  wp = [1.0 1.0 0.0 0.0]; % weights for phases used in origin time calculation
  
  % otp, ots, otm, otmm holds origin time for phase reading for each grid
  % point
  
  % find weighted avarage origin time for each grid point
  ot  = zeros(size(otp));
  ote = zeros(size(otp));
  if (dt(1,1) && wp(1) ~= 0) 
    ot = ot + otp .* wp(1);

  end
  if (dt(1,2) && wp(2) ~= 0)
    ot = ot + ots .* wp(2);

  end
  if (dt(1,3) && wp(3) ~= 0)
    ot = ot + otm .* wp(3);

  end
  if (dt(1,4) && wp(4) ~= 0)
    ot = ot + otmm .* wp(4);
  end
  
  n  = dt(1,1:4) * wp';
  ot = ot ./ n;

  otime(:,:,:,ks) = ot;
  
  % find weighted error for origin time and phase arrival
  we = wp; % weights for error calculation (currently the same as wp, could be tweaked)
  if (dt(1,1) && we(1) ~= 0) 
    ote = ote + we(1) .* (otp - ot).^2;
  end
  if (dt(1,2) && we(2) ~= 0)
    ote = ote + we(2) .* (ots - ot).^2;
  end
  if (dt(1,3) && we(3) ~= 0)
    ote = ote + we(3) .* (otm - ot).^2;
  end
  if (dt(1,4) && we(4) ~= 0)
    ote = ote + we(4) .* (otmm - ot).^2;
  end
  
  ote = sqrt(ote ./ sum(we));
  oterrs(:,:,:,ks) = ote; % will be equal to dterr
  t_otime = t_otime + ot;
  
  % Store travel times
  tttime(:,:,:,ks,1) = dtp;
  tttime(:,:,:,ks,2) = dts;
  tttime(:,:,:,ks,3) = dtm;
  tttime(:,:,:,ks,4) = dtmm;
  
  [minot, i_minot] = min(ote(:));
  fprintf ('sr: station %i, minimum weighted rms, origin time: %f.\n', stations(ks,3), minot);
  
  tot_err = sqrt (ote.^2 + dterr.^2);
  [mintt_, i_tot] = min(tot_err(:));
  fprintf ('sr: station %i, minimum total rms: %f.\n', stations(ks,3), mintt_);
  % plot best fit
%   plot3 (rb(i_mindt,1), rb(i_mindt,2), rb(i_mindt,3),  'xy');
  plot3 (rgx(i_tot), rgy(i_tot), rgz(i_tot),  'xr');
end

%% Combine total errors
dtn = dt(1,2) + dt(1,3) + dt(1,4) + dt(2,3) + dt(2,4) + dt(3,4);
t_dterr = sqrt (t_dterr ./ (3.*dtn));

% Avarage origin time for all grid points for all stations (wp is weight to
% phases).
t_otime     = t_otime ./ 3;

% Find traveltime error for all grid points for all stations and all phases
wp = [1.0 1.0 1.0 1.0];
n  = 0;
if (dt(1,1) && wp(1) ~= 0)
  ph = 1; % phase
  t_otimeerr  = t_otimeerr + ...
                ((origtime(1,ph) - t_otime) - tttime(:,:,:,1,ph)).^2 + ...
                ((origtime(2,ph) - t_otime) - tttime(:,:,:,2,ph)).^2 + ...
                ((origtime(3,ph) - t_otime) - tttime(:,:,:,3,ph)).^2;
  n = n + 3 .* wp(ph);
end
if (dt(1,2) && wp(2) ~= 0)
  ph = 2; % phase
  t_otimeerr  = t_otimeerr + ...
                ((origtime(1,ph) - t_otime) - tttime(:,:,:,1,ph)).^2 + ...
                ((origtime(2,ph) - t_otime) - tttime(:,:,:,2,ph)).^2 + ...
                ((origtime(3,ph) - t_otime) - tttime(:,:,:,3,ph)).^2;
  n = n + 3 .* wp(ph);
end
if (dt(1,3) && wp(3) ~= 0)
  ph = 3; % phase
  t_otimeerr  = t_otimeerr + ...
                ((origtime(1,ph) - t_otime) - tttime(:,:,:,1,ph)).^2 + ...
                ((origtime(2,ph) - t_otime) - tttime(:,:,:,2,ph)).^2 + ...
                ((origtime(3,ph) - t_otime) - tttime(:,:,:,3,ph)).^2;
  n = n + 3 .* wp(ph);
end
if (dt(1,4) && wp(4) ~= 0)
  ph = 4; % phase
  t_otimeerr  = t_otimeerr + ...
                ((origtime(1,ph) - t_otime) - tttime(:,:,:,1,ph)).^2 + ...
                ((origtime(2,ph) - t_otime) - tttime(:,:,:,2,ph)).^2 + ...
                ((origtime(3,ph) - t_otime) - tttime(:,:,:,3,ph)).^2;
  n = n + 3 .* wp(ph);
end

t_otimeerr = sqrt(t_otimeerr ./ n);

%t_otimeerr  = sqrt ( ((t_otime - otime(:,:,:,1)).^2 + (t_otime - otime(:,:,:,2)).^2 + (t_otime - otime(:,:,:,3)).^2) ./ 3 );




% subplot (1,2,2);
% %   slice (ttderr, 3,2,3);
% S = [rb t_ttderr(:)];
% 
% Sn = S(isnan(S(:,4)),:);
% S  = S(~isnan(S(:,4)),:);
% hold on;
% 
% % plot HYPOSAT quakes
% scatter3(quakes(:,1), quakes(:,2), [6000 6000], 100, 'ob', 'filled');
% 
% scatter3(S(:,1), S(:,2), S(:,3), 100, S(:,4), 'filled', 'MarkerEdgeColor', 'k');
% %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'MarkerEdgeColor', 'k');
% %scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'rx');
% 
% 
% 
% xlabel ('x [m] UPS');
% ylabel ('y [m] UPS');
% zlabel ('depth [m]');
% colormap (hot); colorbar;
% set(gca, 'ZDir', 'reverse');
% view (az, el);
% title ('Combined RMS: Differences');

mindt   = min(min(min(t_dterr)));
minot   = min(min(min(t_otimeerr)));
% minttd  = min(min(min(t_ttderr)));
fprintf ('sr: combined minimum rms, differences: %f.\n', mindt);
fprintf ('sr: combined minimum rms, origin: %f.\n', minot);

wt = [1.0 1.0]; % weights for origin and difference errors (todo: certify rms calc is correct)

%total_err = sqrt((wt(2) .* t_dterr.^2 + wt(1) .* t_otimeerr.^2) ./ sum(wt));
total_err = (wt(2) .* dtn.* t_dterr + wt(1) .* n .* t_otimeerr) ./ (wt(2) .* dtn + wt(1) .* n);
[mintt, mintt_i] = min(total_err(:));

fprintf ('sr: combined minimum rms, total: %f.\n', mintt);

% error at chosen point
fprintf ('sr: best total point, onset error: %f, difference error: %f\n', t_otimeerr(mintt_i), t_dterr(mintt_i));

%% Plot result
az = 0; el = 90;

figure(fig(2)); clf('reset');
% subplot(1,2,1);
%   slice (oterr, 1:5,1:5,1:5);
S = [rgx(:) rgy(:) rgz(:) total_err(:)];

Sn = S(isnan(S(:,4)),:);
S  = S(~isnan(S(:,4)),:);
hold on;
h = nan(4,1);


scatter3(S(:,1), S(:,2), S(:,3), 100, S(:,4), 'filled', 'MarkerEdgeColor', 'k');
%scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'MarkerEdgeColor', 'k');
%scatter3(Sn(:,1), Sn(:,2), Sn(:,3), 100, Sn(:,4), 'rx');

% plot stations
cc = ['g', 'y', 'k'];

for i=1:nstations
  h(i) = scatter(stations(i,1), stations(i,2),  sprintf('%co', cc(stations(i,3)+1)), 'filled', 'SizeData', 100);
  plot3(stations(i,1), stations(i,2), -1, 'kx');
end




xlabel ('x [m] UPS');
ylabel ('y [m] UPS');
zlabel ('depth [m]');
colormap (hot); colorbar;
set(gca, 'ZDir', 'reverse');
view (az, el);

title (sprintf('%s: Combined RMS: Total error', figtitle));

[x, y, z] = ind2sub(size(total_err), mintt_i);

% residuals of phases
fprintf ('sr: residuals of onset times:\n');
for ks=1:nstations
  for ph=1:nphases
    s   = stations(ks,3);
    res = ((origtime(ks,ph) - t_otime(mintt_i)) - tttime(x,y,z,ks,ph));
    fprintf ('sr: station: %d, phase: %d, onset time residual: %f\n', s, ph, res);
  end
end

% calculated travel times
fprintf ('sr: calculated travel times:\n');
for ks=1:nstations
  for ph=1:nphases
    s   = stations(ks,3);
    res = tttime(x,y,z,ks,ph);
    fprintf ('sr: station: %d, phase: %d, travel time: %f\n', s, ph, res);
  end
end

%% Best match
% grid point: rms
% t_dterr_f = t_dterr(:);
% mindt_i   = find(t_dterr_f == mindt);
% s = [rb(mindt_i,:) t_dterr_f(mindt_i)];

% t_otimeerr_f = t_otimeerr(:);
% minot_i   = find(t_otimeerr_f == minot);
% s = [rb(minot_i,:) t_otimeerr_f(minot_i) otime(minot_i)];
% mindt_i = minot_i;

t_otimeerr = total_err;

s = [rgx(mintt_i) rgy(mintt_i) rgz(mintt_i) total_err(mintt_i) (t_otime(mintt_i) / 24 / 60 / 60) mintt_i];
mindt_i = mintt_i;

h(5) = scatter3(rgx(mintt_i), rgy(mintt_i), rgz(mintt_i), 'co', 'filled', 'SizeData', 100);

% plot HYPOSAT quakes
h(4) = scatter3(quakes(:,1), quakes(:,2),  6000 .* ones(size(quakes,1),1), 100, 'ob', 'filled');

legend (h, 'GAK2', 'GAK3', 'GAK4', 'HYPO', 'Solution');
end
