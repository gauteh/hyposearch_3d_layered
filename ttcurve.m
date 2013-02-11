%% Plot traveltime curves


% use depth of earthquake (s_ps)


[rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);


% Use first stations
s = stations(1,:);

% Calculate from s to quake
%e = s_ps(1:3);
e = [20000 20000 5000];
s(3) = e(3);

r_ray_p = traveltimes{1,2};
r_ray_s = traveltimes{1,4};
r_ray_m = traveltimes{1,6};
r_ray_mm = traveltimes{1,8};

% figure out grid cells along s to e
rgxp = 1:rnx;
rgyp = 1:rny;
rgzp = 1:rnz;

[rgxp, rgyp, rgzp] = meshgrid (rgxp, rgyp, rgzp);


ix = linspace(0,1,1000);
x  = [0 1];

iyx = interp1(x, [s(1) e(1)], ix, 'linear');
iyy = interp1(x, [s(2) e(2)], ix, 'linear');
iyz = interp1(x, [s(3) e(3)], ix, 'linear');

iy = [iyx' iyy' iyz'];

x = interp3(rgx, rgy, rgz, rgxp, iyx, iyy, iyz, 'nearest');
y = interp3(rgx, rgy, rgz, rgyp, iyx, iyy, iyz, 'nearest');
z = interp3(rgx, rgy, rgz, rgzp, iyx, iyy, iyz, 'nearest');

arg = [x' y' z'];

% remove duplicates
arg = unique (arg, 'rows');

x = [];
l = [];
d = [];
tp = [];
ts = [];
tm = [];
tmm = [];
for k=arg'
  x = cat(1, x, [rgx(k(1), k(2), k(3)) rgy(k(1), k(2), k(3)) rgz(k(1), k(2), k(3))]);
  l = cat(1, l, r_ray_p (k(1), k(2), k(3), 2));
  tp = cat(1, tp, r_ray_p (k(1), k(2), k(3), 3));
  ts = cat(1, ts, r_ray_s (k(1), k(2), k(3), 3));
  tm = cat(1, tm, r_ray_m (k(1), k(2), k(3), 3));
  tmm = cat(1, tmm, r_ray_mm (k(1), k(2), k(3), 3));
end

d = vectnorm(x - repmat(s,size(x,1),1));
d = d ./ 1000;


figure(1); clf('reset');
plot (d, tp, 'b'); hold on;
plot (d, ts, 'r');
plot (d, tm, 'y');
plot (d, tmm, 'g');

legend ('P', 'S', 'M', 'MM');
title (sprintf('Traveltimes from station GAK2 (along depth: %g m)', e(3)));
xlabel ('Distance [km]');
ylabel ('Time [s]');
grid on;

%% Plot difference of P - M and P - MM as function of depth

% Calculate from s to quake
e = s_ps(1:3);
e = [2000 2000 7000];
e(3) = 7000;
s = e;
s(3) = 0;

ix = linspace(0,1,1000);
x  = [0 1];

iyx = interp1(x, [s(1) e(1)], ix, 'linear');
iyy = interp1(x, [s(2) e(2)], ix, 'linear');
iyz = interp1(x, [s(3) e(3)], ix, 'linear');

iy = [iyx' iyy' iyz'];

x = interp3(rgx, rgy, rgz, rgxp, iyx, iyy, iyz, 'nearest');
y = interp3(rgx, rgy, rgz, rgyp, iyx, iyy, iyz, 'nearest');
z = interp3(rgx, rgy, rgz, rgzp, iyx, iyy, iyz, 'nearest');

arg = [x' y' z'];

% remove duplicates
arg = unique (arg, 'rows');
n = isnan(arg(:,1)) | isnan(arg(:,2)) | isnan(arg(:,3));
arg = arg(~n,:);

x = [];
l = [];
d = [];
tp = [];
ts = [];
tm = [];
tmm = [];
for k=arg'
  x = cat(1, x, [rgx(k(1), k(2), k(3)) rgy(k(1), k(2), k(3)) rgz(k(1), k(2), k(3))]);
  l = cat(1, l, r_ray_p (k(1), k(2), k(3), 2));
  tp = cat(1, tp, r_ray_p (k(1), k(2), k(3), 3));
  ts = cat(1, ts, r_ray_s (k(1), k(2), k(3), 3));
  tm = cat(1, tm, r_ray_m (k(1), k(2), k(3), 3));
  tmm = cat(1, tmm, r_ray_mm (k(1), k(2), k(3), 3));
end

d = vectnorm(x - repmat(s,size(x,1),1));
%d = d ./ 1000;


figure(2); clf('reset');
subplot (1,2,1);
plot (d, tp, 'b'); hold on;
plot (d, ts, 'r');
plot (d, tm, 'y');
plot (d, tmm, 'g');

legend ('P', 'S', 'M', 'MM');
title (sprintf('Traveltimes from station GAK2 (distance from station: %g m)', norm(stations(1,1:2) - e(1:2))));
xlabel ('Detph [m]');
ylabel ('Time [s]');
grid on;

xlim ([5500 7000]);
subplot (1,2,2);
plot (d, ts - tp, 'k'); hold on;
plot (d, tm-tp, 'b');
plot (d, tmm-tp, 'r');
plot (d, tm - ts, 'g');
plot (d, tmm - ts, 'y');
xlim ([5500 7000]);

legend ('S - P', 'M - P', 'MM - P', 'M - S', 'MM - S');
title (sprintf('Traveltimes differences from station GAK2 (distance from station: %g m)', norm(stations(1,1:2) - e(1:2))));
xlabel ('Detph [m]');
ylabel ('Time [s]');
grid on;





