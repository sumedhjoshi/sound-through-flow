% data = sem_2D_inhomogenous( n, mx, my, x, y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 );
function data = sem_2D_inhomogenous( n, mx, my, x, y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 )

   % Set some grid parameters.
   Lx  = max(x(:)) - min(x(:));
   Ly  = max(y(:)) - min(y(:));

   % Set some time stepping parameters.
   dt   = t(2) - t(1);
   tmax = t(end);

   % Set the grid size.
   r = n * n * mx * my;

   % Set some grid parameters.
   xlim = [min(x(:)) max(x(:))];
   ylim = [min(y(:)) max(y(:))];

   % Build the GLL grid in one element.
   [eta junk junk] = lglnodes( n - 1 );
   eta = -eta;

   % Get the grid element width.
   hx = Lx / mx;
   hy = Ly / my;

   % Set the bulk modulus.
   c    = 343;
   rho0 = 1.02;
   beta = rho0 * c^2;

   % Set some time stepping parameters.
   dt   = t(2) - t(1);
   tmin = t(1);
   tmax = t(end);

   % Compute the spectral differentiation matrix of order 1.
   d1 = zeros( n, n );
   for ii = 1: n
      d1(:, ii) = make_lagrange( eta', eta', ii, 1 );
   end

      % Compute the spectral differentiation matrix in x.
      D1x = ( 2 / hx ) * kron( eye(mx,mx), d1 );

      % Update the interfacial boundary conditions to be identical.
      for ii = 1:mx - 1
         left  = ( ii - 1 ) * n + n;
         right = left + 1;
         ndx = n * ii;
         interface = ( D1x( left:right, : ) + flipud( D1x( left:right, :) ) );
         interface( 1:2, left:right ) =  interface( 1:2, left:right) + 0 * ( 2 / hx )^2 * [ 1.0, 1.0; 1.0, 1.0];
         D1x( left:right, : ) =  interface/2;
      end

      % Update the global matrix for Dirichlet boundary conditions.
      D1x(1,1)     =  2/hx * 0.5;
      D1x(end,end) = -2/hx * 0.5;

      % Compute the spectral differentiation matrix in y.
      D1y = ( 2 / hy ) * kron( eye(mx,mx), d1 );

      % Update the interfacial boundary conditions to be identical.
      for ii = 1:my - 1
         left  = ( ii - 1 ) * n + n;
         right = left + 1;
         ndx = n * ii;
         interface = ( D1y( left:right, : ) + flipud( D1y( left:right, :) ) );
         interface( 1:2, left:right ) =  interface( 1:2, left:right) + 0 * ( 2 / hy )^2 * [ 1.0, 1.0; 1.0, 1.0];
         D1y( left:right, : ) =  interface/2;
      end

      % Update the global matrix for Dirichlet boundary conditions.
      D1y(1,1)     =  2/hy * 0.5;
      D1y(end,end) = -2/hy * 0.5;

   % Build the 2D x and y differentiation matrices.
   D1x = kron( sparse( D1x ), speye( n * mx ));
   D1y = kron( speye( n * my ), sparse(D1y) );

   % Set some initial conditions.
   u0x = ux0;
   u0y = uy0;

   % Set the background field.
   Vx = zeros( r, 1 );
   Vy = zeros( r, 1 );
   Vx(:) = vx0(:);
   Vy(:) = vy0(:);

   % Start time stepping.
   Ux = zeros(r,length(t));
   Uy = zeros(r,length(t));
   S  = zeros(r,length(t));
   Ux(:,1) = u0x(:);
   Uy(:,1) = u0y(:);
   S(:,1)  = s0(:);
   for ii = 2:length(t)

      p(ii,length(t),1);

      % Build the condensation.
      S(:,ii)  =  S(:,ii-1) - dt * ( D1x * Ux(:,ii-1) + D1y * Uy(:,ii-1) ) ...
                            - dt * Vx(:,1) .* ( D1x * S(:,ii-1) ) ...
                            - dt * Vy(:,1) .* ( D1y * S(:,ii-1) );

      % Build the velocity in x.
      Ux(:,ii) = Ux(:,ii-1) - dt * beta * ( D1x * S(:,ii-1) ) / rho0  ...
                            - dt * Vx(:, 1) .* (D1x * Ux(:,ii-1)) ...
                            - dt * Vy(:, 1) .* (D1y * Ux(:,ii-1)) ...
                            - dt * Ux(:, ii-1) .* (D1x * Vx(:,1)) ...
                            - dt * Uy(:, ii-1) .* (D1y * Vx(:,1)) ...
                            - dt * S(:,ii-1) .* Vx(:, 1) .* (D1x * Vx(:,1)) ...
                            - dt * S(:,ii-1) .* Vy(:, 1) .* (D1y * Vx(:,1));

      % Build the velocity in y.
      Uy(:,ii) = Uy(:,ii-1) - dt * beta * ( D1y * S(:,ii-1) ) / rho0   ...
                            - dt * Vx(:, 1) .* ( D1x * Uy(:,ii-1)) ...
                            - dt * Vy(:, 1) .* ( D1y * Uy(:,ii-1)) ...
                            - dt * Ux(:, ii-1) .* (D1x * Vy(:,1)) ...
                            - dt * Uy(:, ii-1) .* (D1y * Vy(:,1)) ...
                            - dt * S(:,ii-1) .* Vx(:, 1) .* (D1x * Vy(:,1)) ...
                            - dt * S(:,ii-1) .* Vy(:, 1) .* (D1y * Vy(:,1));

   end

   % Pack the data for output.
   rx = n * mx;
   ry = n * my;
   data.S  = reshape( S,  [ ry, rx, length(t) ] );
   data.Ux = reshape( Ux, [ ry, rx, length(t) ] );
   data.Uy = reshape( Uy, [ ry, rx, length(t) ] );
   data.t  = t;
   data.x  = x;
   data.y  = y;
   data.physics.c0   = c0;
   data.physics.Vx   = reshape( Vx, ry, rx );
   data.physics.Vy   = reshape( Vy, ry, rx );
   data.physics.rho0 = rho0;

end




