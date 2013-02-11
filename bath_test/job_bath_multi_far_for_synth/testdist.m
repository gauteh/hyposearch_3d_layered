% test distance on coordinates

% cartesian
qc = load('quakes.coor');
qc = qc(1:2);

sc = load('stations.coor');
sc = sc(:,1:2);

% geographic
qd = load('quakes.d');
qd = qd(1:2);

sd = load('stations.d');
sd = sd(:,1:2);

% cartesian distance
qc = repmat(qc, 3, 1);
distc = vectnorm(qc - sc);

fprintf ('Cartesian distance: %g m\n', distc);

% geographic distance
qd = repmat(qd, 3, 1);
distg = distance(sd(:,[2 1]), qd(:,[2 1]));

fprintf ('Geographic distance: %g deg\n', distg);

% geographic distance to cartesian distance
distgc = deg2km (distg) .* 1000;
fprintf ('Geographic to cartesian distance: %g m\n', distgc);

% error
err = distc - distgc;

fprintf ('Absolute error: %g m\n', err);