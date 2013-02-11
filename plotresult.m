function plotresult (interface_tris, stations, s_ps, s_psm, s_psmm, quakes, traveltimes, R_tri, rnx, rny, rnz, fig)
  % Plot result of eq localization
  
  figure (fig); clf('reset');
  
  ninterfaces = size(interface_tris,1);
  for i = 1:3
    trisurf (interface_tris{i}, 'FaceAlpha', 0.8);
    hold on;
  end

  title ('Bathymetry triangulation (depth positive, constrained region)');
  set(gca, 'ZDir', 'reverse');
  colorbar;
  
  shading interp;
  
  hold on
  
  % Plot region of interest
  %trisurf (R_tri, 'FaceAlpha', 0.2);

  [rgx, rgy, rgz, ep_x, ep_y, ep_z, rdx, rdy, rdz] = makegrid (R_tri, rnx, rny, rnz);
  
  % convert to kms
  
  % Plot stations
  nstations = size(stations,1);
  
  cc = ['g', 'y', 'k'];
  h = [];
  for i=1:nstations
    h(i) = scatter(stations(i,1), stations(i,2),  sprintf('%co', cc(stations(i,3)+1)), 'filled');
    
    %% plot phases
    %% P & S: P
    ray_p = traveltimes{i,1}; % rays
    rp    = traveltimes{i,2}; % grid point <-> ray match
    ray_s = traveltimes{i,3};
    rs    = traveltimes{i,4};   
    
    if (~isnan(s_ps(4)))


      %ind  = rbg(s_ps(6),:);
      %ix = rbg(s_ps(6),1); iy = rbg(s_ps(6),2); iz = rbg(s_ps(6), 3);
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s_ps(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'rx', 'LineWidth', 2);

      % S:
 
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b', 'LineWidth', 2);
    end
    
    ray_m = traveltimes{i,5};
    rm    = traveltimes{i,6};
    
    %% P, S & M
    %ix = rbg(s_psm(6),1); iy = rbg(s_psm(6),2); iz = rbg(s_psm(6), 3);
    if (~isnan(s_psm(4)))
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s_psm(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'r', 'LineWidth', 2);

      % S:   
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b', 'LineWidth', 2);

      % M

      rind = rm(ix, iy, iz, 1);
      ray  = ray_m(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'g', 'LineWidth', 2);
    end
    
    %% P, S, M & MM
    %ix = rbg(s_psmm(6),1); iy = rbg(s_psmm(6),2); iz = rbg(s_psmm(6), 3);
    if ~isnan(s_psmm(4))
      [ix, iy, iz] = ind2sub(size(rp(:,:,:,1)), s_psmm(6));
      rind = rp(ix, iy, iz, 1);
      ray  = ray_p(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'r', 'LineWidth', 2);

      % S:   
      rind = rs(ix, iy, iz, 1);
      ray  = ray_s(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'b', 'LineWidth', 2);

      % M
      rind = rm(ix, iy, iz, 1);
      ray  = ray_m(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'g', 'LineWidth', 2);    

      % MM
      ray_mm = traveltimes{i,7};
      rmm    = traveltimes{i,8};
      rind = rmm(ix, iy, iz, 1);
      ray  = ray_mm(rind,:,:);

      x = ray(:,4,:);
      y = ray(:,5,:);
      z = ray(:,6,:);
      plot3 (x(:), y(:), z(:), 'y', 'LineWidth', 2);      
    end
  end

  

  % Plot locations
  cc = ['r', 'b', 'm'];
  if (~isnan(s_ps(4)))
    h(i+1) = scatter3 (s_ps(:,1), s_ps(:,2), s_ps(:,3),  sprintf('%co', cc(1)), 'filled');
  end
  if (~isnan(s_psm(4)))
    h(i+2) = scatter3 (s_psm(:,1), s_psm(:,2), s_psm(:,3),  sprintf('%co', cc(2)), 'filled');
  end
  if (~isnan(s_psmm(4)))
    h(i+3) = scatter3 (s_psmm(:,1), s_psmm(:,2), s_psmm(:,3),  sprintf('%co', cc(3)), 'filled');
  end
  
  
  xlabel('x [m] UPS');
  ylabel('y [m] UPS');
  
  
  
  % Plot hyposat locations
  cc = ['k', 'k'];
  j = i;
  for i=1:size(quakes,1)
    h(j + 3 + i) = scatter3(quakes(i,1), quakes(i,2), 5000, sprintf('%co', cc(quakes(i,3)+1)), 'filled');
  end
  
  legend (h, 'GAK2', 'GAK3', 'GAK4', 'PS', 'PSM', 'PSMM', 'HYPO1', 'HYPO2');
end