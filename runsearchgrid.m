% Use only P and S
u = [ 1 1 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0];

fprintf ('hs: calculating RMS and best fit for P and S phases only..\n');
% mindt_i, s, t_dterr, dterrs, t_oterr, oterrs, otime, t_otime, t_otimeerr
[i_ps, s_ps, dterr_ps, dterrs_ps, oterr_ps, tttime_ps, oterrs_ps, t_otime_ps, t_otimeerr_ps] = searchgrid (R_tri, rnx, rny, rnz, stations, quakes, phases, traveltimes, u, 5:6, 'P and S');


% Use P, S and M
u = [ 1 1 1 0;
      1 1 1 0;
      0 0 0 0;
      0 0 0 0];

fprintf ('hs: calculating RMS and best fit for P, S and M phases only..\n');
[i_psm, s_psm, dterr_psm, dterrs_psm, oterr_psm, tttime_psm, oterrs_psm, t_otime_psm, t_otimeerr_psm]  = searchgrid (R_tri, rnx, rny, rnz, stations, quakes, phases, traveltimes, u, 7:8, 'P, S and M');

% Use P, S, M and MM
u = [ 1 1 1 1;
      1 1 1 1;
      1 1 1 1;
      1 1 1 1];
    
fprintf ('hs: calculating RMS and best fit for all phases..\n');
[i_psmm, s_psmm, dterr_psmm, dterrs_psmm, oterr_psmm, tttime_psmm, oterrs_psmm, t_otime_psmm, t_otimeerr_psmm]  = searchgrid (R_tri, rnx, rny, rnz, stations, quakes, phases, traveltimes, u, 9:10, 'P, S, M and MM');
