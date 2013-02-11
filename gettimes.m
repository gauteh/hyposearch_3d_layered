function [so, p, s, m, mm, ap, as, am, amm, alla] = gettimes (R_tri, rnx, rny, rnz, stations, quake, traveltimes)
% gettimes.m: get travel times for quake
%
% quake is desired location


[rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);

%% Figure out closest point
[rxx, ryy, rzz] = meshgrid (1:rnx, 1:rny, 1:rnz);
x = interp3 (rgx, rgy, rgz, rxx, quake(1), quake(2), quake(3), 'nearest');
y = interp3 (rgx, rgy, rgz, ryy, quake(1), quake(2), quake(3), 'nearest');
z = interp3 (rgx, rgy, rgz, rzz, quake(1), quake(2), quake(3), 'nearest');

assert (all(~isnan([x y z])), 'chosen quake point is not in model');

ap = ones(size(rgx));
as = ones(size(rgx));
am = ones(size(rgx));
amm = ones(size(rgx));

for ks=1:size(stations,1)
  %% Phase P
  ray_p = traveltimes{ks,1};
  rp    = traveltimes{ks,2};
  dtp = rp(:,:,:,3);  % travel time
  ap = ap + dtp;
  %%

  %% Phase S
  ray_s = traveltimes{ks,3};
  rs    = traveltimes{ks,4};
  dts = rs(:,:,:,3);  % travel time
  as = as + dts;
  %%

  %% Phase M
  ray_m = traveltimes{ks,5};
  rm    = traveltimes{ks,6};
  dtm = rm(:,:,:,3);  % travel time
  am = am + dtm;
  %%

  %% Phase MM
  ray_mm = traveltimes{ks,7};
  rmm    = traveltimes{ks,8};
  dtmm = rmm(:,:,:,3);  % travel time
  amm = amm + dtmm;
  %%

  p(ks) = dtp(x, y, z);
  s(ks) = dts(x, y, z);
  m(ks) = dtm(x, y, z);
  mm(ks) = dtmm(x, y, z);
end

so = [rgx(x,y,z) rgy(x,y,z) rgz(x,y,z) 0.0 0.0 sub2ind(size(rgx), x, y, z)];
fprintf ('gt: used quake hypocenter: %g, %g, %g\n', so(1:3));

% calculated travel times
fprintf ('gt: calculated travel times:\n');
for ks=1:size(stations,1)
  st   = stations(ks,3);
  fprintf ('gt: station: %d, phase: 0, travel time: %f\n', st, p(ks));
  fprintf ('gt: station: %d, phase: 1, travel time: %f\n', st, s(ks));
  fprintf ('gt: station: %d, phase: 2, travel time: %f\n', st, m(ks));
  fprintf ('gt: station: %d, phase: 3, travel time: %f\n', st, mm(ks));
end

% calculate for which available ttimes
alla = ap + as + am + amm;
ap   = find(~isnan(ap(:)));
as   = find(~isnan(as(:)));
am   = find(~isnan(am(:)));
amm  = find(~isnan(amm(:)));
alla = find(~isnan(alla(:)));

fprintf ('gt: times available for P: %d\n', length(ap));
fprintf ('gt: times available for S: %d\n', length(as));
fprintf ('gt: times available for M: %d\n', length(am));
fprintf ('gt: times available for MM: %d\n', length(amm));
fprintf ('gt: times available for all: %d\n', length(alla));

end

