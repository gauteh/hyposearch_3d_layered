function plotslice (bathymetry, stations, traveltimes, quakes, R_tri, rnx, rny, rnz, fig, s, err, ph)
  % Plot result of eq localization
  %
  % ph = 1, P & S phase
  %    = 2, P, S & M phase
  %    = 3, P, S, M & MM phase
  
  figure (fig); clf('reset');
  
  %trisurf (bathymetry, 'FaceAlpha', 0.5);
  title ('Bathymetry triangulation (depth positive, constrained region)');
  set(gca, 'ZDir', 'reverse');
  colorbar;
  
  hold on
  
  % Plot region of interest
  %trisurf (R_tri, 'FaceAlpha', 0.2);

  [rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);
  
  % Plot stations
  nstations = size(stations,1);
  
  cc = ['g', 'y', 'k'];
  h = [];
  for i=1:nstations
    h(i) = scatter(stations(i,1), stations(i,2),  sprintf('%co', cc(stations(i,3)+1)), 'filled');
   
    ray_p = traveltimes{i,1}; % rays
    rp    = traveltimes{i,2}; % grid point <-> ray match
    ray_s = traveltimes{i,3};
    rs    = traveltimes{i,4};  
    ray_m = traveltimes{i,5};
    rm    = traveltimes{i,6};
    ray_mm = traveltimes{i,7};
    rmm    = traveltimes{i,8};
    
    %% plot phases
    if (ph == 1)
      title ('Error and location, P & S');
      %% P & S: P
      %ind  = rbg(s_ps(6),:);
      %ix = rbg(s_ps(6),1); iy = rbg(s_ps(6),2); iz = rbg(s_ps(6), 3);
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'r');

      % S:
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b');
      
    elseif (ph == 2)
      title ('Error and location, P, S & M');
      %% P, S & M
      %ix = rbg(s_psm(6),1); iy = rbg(s_psm(6),2); iz = rbg(s_psm(6), 3);
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'r');

      % S:   
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b');

      % M
      rind = rm(ix, iy, iz, 1);
      ray  = ray_m(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'g');

    elseif (ph == 3)
      title ('Error and location, P, S, M & MM');
      %% P, S, M & MM
      %ix = rbg(s_psmm(6),1); iy = rbg(s_psmm(6),2); iz = rbg(s_psmm(6), 3);
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'r');

      % S:   
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b');

      % M
      rind = rm(ix, iy, iz, 1);
      ray  = ray_m(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'g');    

      % MM
      rind = rmm(ix, iy, iz, 1);
      ray  = ray_mm(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'y');
    end
  end
  

  % Plot locations
  h(i+1) = scatter3 (s(:,1), s(:,2), s(:,3),  'co', 'filled');

  
  xlabel('x [m] UPS');
  ylabel('y [m] UPS');
  
  % Plot error
  slice (rgx, rgy, rgz, err, [0 s(1)], [0 s(2)], s(3));
  colormap (hot); colorbar;
  
  grid on;
  
  % Plot hyposat locations
  cc = ['b', 'k'];
  j = i;
  for i=1:size(quakes,1)
    h(j + 3 + i) = scatter3(quakes(i,1), quakes(i,2), 5000, sprintf('%co', cc(quakes(i,3)+1)), 'filled');
  end
  
  legend (h, 'GAK2', 'GAK3', 'GAK4', 'Solution', 'HYPO1', 'HYPO2');
end