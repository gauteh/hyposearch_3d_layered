%% Run hyposearch job from job file in current directory and save result

% Start logging
diary ('hs_log.txt');
diary on;

runt = datestr(now);
fprintf ('hs: run: %s\n', runt);

[~, hostname] = system ('hostname');
hostname = strtrim (hostname);

if (strcmp (hostname, 'billy.geo.uib.no'))
  desired_workers = 8;
elseif (strcmp (hostname, 'Odd'))
  desired_workers = 8;
else
  desired_workers = 3; % no point in having more than the number of stations * phases
end

fprintf ('hs: host: %s\n', hostname);
fprintf ('hs: workers: %i\n', desired_workers);

% Set up parallell workers
cworkers = matlabpool ('size');
if (cworkers == 0)
  matlabpool ('open', desired_workers);
end

% Set up job
s_t0 = tic;
[bounds, interface_tris, stations, quakes, phases, velp, vels, R_tri, usephases, event, jobname, rnx, rny, rnz] = setuphsjob (false);
s_t  = toc(s_t0);

% Run Hyposearch to solve grid
h_t0 = tic;
traveltimes = hyposearch (bounds, interface_tris, stations, velp, vels, R_tri, rnx, rny, rnz);
h_t  = toc(h_t0);

% Search grid and detemine RMS values for various combinations of phase
% uses
sr_t0 = tic;

runsearchgrid;

sr_t  = toc(sr_t0);
t_t = toc(s_t0);

% Summarize
fprintf ('Hyposearch job done:\n\n');

fprintf ('Event ..: %s\n', event);
fprintf ('Job ....: %s\n', jobname);
fprintf ('Time ...: %s\n\n', runt);

fprintf ('----------------------\n');
fprintf ('Setup time ..........: %f secs\n', s_t);
fprintf ('Hyposearch time .....: %f secs\n', h_t);
fprintf ('Search grid time ....: %f secs\n', sr_t);
fprintf ('----------------------\n');
fprintf ('Total time ..........: %f secs\n', t_t);
fprintf ('----------------------\n\n');

% Save result and write output file
fprintf ('Result and solved rays saved in ..: hs_result.mat\n');
save ('hs_result.mat', '-mat', '-v7.3');

% Save figures
fprintf ('Figures stored in ................: hs_figures.fig\n');
figs = 5:10;
hgsave (figs, 'hs_figures.fig');

fprintf ('Location and error stored in .....: hs_out\n');
fid = fopen ('hs_out', 'w');
fprintf (fid, '#1 %s\n', event);    % write event
fprintf (fid, '#2 %s\n', jobname);  % write job name
fprintf (fid, '#3 %s\n', runt);     % write date and time

% write a line for each solution of the different phase usage
% first number details what phases have been used in solution
if (~any(isnan(s_ps)))
  for k=1:size(s_ps,1)
    fprintf (fid, '%f  %f %f  %f 0', s_ps(k,1:4));
    fprintf (fid, ' %s # P and S\n', datestr(s_ps(k,5), 'yyyy-mm-dd HH:MM:SS.FFF'));
  end
end
if (~any(isnan(s_psm)))
  for k=1:size(s_psm,1)
    fprintf (fid, '%f  %f %f  %f 1', s_psm(k,1:4));
    fprintf (fid, ' %s # P, S and M\n', datestr(s_ps(k,5), 'yyyy-mm-dd HH:MM:SS.FFF'));
  end
end
if (~any(isnan(s_psmm)))
  for k=1:size(s_psmm,1)
    fprintf (fid, '%f  %f %f  %f 2', s_psmm(k,1:4));
    fprintf (fid, ' %s # P, S, M and MM\n', datestr(s_ps(k,5), 'yyyy-mm-dd HH:MM:SS.FFF'));
  end
end

fclose (fid);

fprintf ('Output has been logged to ........: hs_log.txt\n');
diary off;
