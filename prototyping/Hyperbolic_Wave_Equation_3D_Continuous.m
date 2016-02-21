%
% Solve the 1D homogenous wave equation in split density/velocity form
% as a coupled set of hyperbolic conservation laws and a spectral
% element discretization.

% This equation is solved in an L^3 box.

   % Set some constants.
   rho  = 1.02;
   c    = 343;
   beta = rho * c^2;

   % Set some grid parameters.
   L  = 800;
   n  = 8;
   mx = 8;

   % Build the GLL grid in one element.
   [x junk junk] = lglnodes( n - 1 );
   eta = -x;
   x   = -x;

   % Compute the global grid.
   hx = L/mx;
   x  = hx * ( x + 1 ) / 2;
   x  = repmat( x, mx, 1 );
   mm = hx * [ 0: mx - 1];
   m2 = repmat(mm,n,1);
   m2 = m2(:);
   x = x + m2;

   % Set the time discretization.
   tmax = 0.25;
   dx   = diff( x );
   dx   = min( dx( dx > 0 ) );
   dt   = (1/20) * dx / c;
   Nt   = ceil( tmax / dt );
   t    = linspace( 0, tmax, Nt );

   % Compute the spectral differentiation matrix of order 1.
   d1 = zeros( n, n );
   for ii = 1: n
      d1(:, ii) = make_lagrange( eta', eta', ii, 1 );
   end

   % Tile the spectral differentiation matrix to build a 1D difference matrix on a multi-domain grid.
   D1 = ( 2 / hx ) * kron( eye(mx,mx), d1 );

   % Update the interfacial boundary conditions to be identical.
   for ii = 1:mx - 1
      left  = ( ii - 1 ) * n + n;
      right = left + 1;
      interface = D1( left:right, : ) + flipud( D1( left:right, :) );
      D1( left:right, : ) = 0.5* interface;
   end

   % Update the global matrix for Dirichlet boundary conditions.
   D1(1,1)     =  2/hx * 0.5;
   D1(end,end) = -2/hx * 0.5;

   % Consistency check.
   u = sin( pi*x/100 );
   Du_a = pi/100*cos( pi*x/100 );
   err = norm( Du_a - D1*u );

   % Assemble the 3D 1st derivative multidomain matrix.
   mz = mx;
   my = mx;
   nx = n * mx;
   ny = n * my;
   nz = n * mz;
   % Dx = Iz \kron Iy \kron D1;
   % Dy = Iz \kron D1 \kron Ix;
   % Dz = D1 \kron Iy \kron Ix;
   Dx = kron( speye( nz * ny ), sparse(D1) );
   Dy = kron( kron( speye( nz ), sparse(D1) ), speye(nx) );
   Dz = kron( sparse(D1), speye( nx*ny ) );

   % Test 3D consistency.
   z = x;
   y = x;
   [X Y Z] = meshgrid( x, y, z );
   X = X(:);
   Y = Y(:);
   Z = Z(:);
   U = sin( pi * X / 100 ) .* sin( pi * Y / 100 ) .* sin( pi * Z / 100 );
   X = reshape( X, [nx, ny, nz] );
   Y = reshape( Y, [nx, ny, nz] );
   Z = reshape( Z, [nx, ny, nz] );
   LU = Dx * U + Dy * U + Dz * U;
   LU_a = cos( pi * X / 100 ) .* sin( pi * Y / 100 ) .* sin( pi * Z / 100 ) * pi / 100 + ...
          sin( pi * X / 100 ) .* cos( pi * Y / 100 ) .* sin( pi * Z / 100 ) * pi / 100 + ...
          sin( pi * X / 100 ) .* sin( pi * Y / 100 ) .* cos( pi * Z / 100 ) * pi / 100;
   LU = reshape( LU , [nx, ny, nz] );
   LU_a = reshape( LU_a , [nx, ny, nz] );

   % Set the initial condition.
   sigma = 25;
   x0 = L/2;
   y0 = L/2;
   z0 = L/2;
   s0 = exp( -( (X - x0).^2 + (Y - y0).^2 + (Z - z0).^2 ) /(2 * sigma^2 ) );
   ux0 = 0 * s0;
   uy0 = 0 * s0;
   uz0 = 0 * s0;
   keyboard;

   % Compute the dt based on frequency of the initial pulse.
   lambda    = 2 * sigma;
   frequency = c / lambda;
   period    = 1 / frequency;

   % Build the operator matrix for eigenvalue identification.
   %r = mx * n;
   %A = eye( 2*r, 2*r );
   %A( 1:r, r+1:2*r) = -dt * D1;
   %A( r+1:2*r, 1:r) = -dt * beta * D1 / rho;

   % Set up time-lagging field matrices.
   Ux0 = ux0(:);
   Ux1 = Ux0;
   Ux2 = Ux0;
   Ux3 = Ux0;
   Uy0 = uy0(:);
   Uy1 = Uy0;
   Uy2 = Uy0;
   Uy3 = Uy0;
   Uz0 = uz0(:);
   Uz1 = Uz0;
   Uz2 = Uz0;
   Uz3 = Uz0;
   S0  = s0(:);
   S1  = S0;
   S2  = S0;
   S3  = S0;

   % Time-propagate.
   for ii = 2:length(t)
      p( ii, Nt, 1);

      % Depending on the time step, use a particular ABk method.
      switch ii
         %case 2
         %   S(:,ii)  = S(:,ii-1)  - dt * ( Dx * Ux(:,ii-1) + Dy * Uy(:,ii-1) + Dz * Uz(:,ii-1) );
         %   Ux(:,ii) = Ux(:,ii-1) - dt * ( beta / rho ) * Dx * S(:,ii-1);
         %   Uy(:,ii) = Uy(:,ii-1) - dt * ( beta / rho ) * Dy * S(:,ii-1);
         %   Uz(:,ii) = Uz(:,ii-1) - dt * ( beta / rho ) * Dz * S(:,ii-1);
         %case 3
         %   S(:,ii) = S(:,ii-1) - dt * Dx * ( 1.5 * Ux(:,ii-1) - 0.5 * Ux(:,ii-2) ) ...
         %                       - dt * Dy * ( 1.5 * Uy(:,ii-1) - 0.5 * Uy(:,ii-2) ) ...
         %                       - dt * Dz * ( 1.5 * Uz(:,ii-1) - 0.5 * Uz(:,ii-2) );
         %   Ux(:,ii) = Ux(:,ii-1) - dt * ( beta / rho ) * Dx * ( 1.5 * S(:,ii-1) - 0.5 * S(:,ii-2) );
         %   Uy(:,ii) = Uy(:,ii-1) - dt * ( beta / rho ) * Dy * ( 1.5 * S(:,ii-1) - 0.5 * S(:,ii-2) );
         %   Uz(:,ii) = Uz(:,ii-1) - dt * ( beta / rho ) * Dz * ( 1.5 * S(:,ii-1) - 0.5 * S(:,ii-2) );
         %otherwise
         %   S(:,ii) = S(:,ii-1) - dt * Dx * ( (23/12) * Ux(:,ii-1) - 4/3 * Ux(:,ii-2)  + 5/12 * Ux(:,ii-3)) ...
         %                       - dt * Dy * ( (23/12) * Uy(:,ii-1) - 4/3 * Uy(:,ii-2)  + 5/12 * Uy(:,ii-3)) ...
         %                       - dt * Dz * ( (23/12) * Uz(:,ii-1) - 4/3 * Uz(:,ii-2)  + 5/12 * Uz(:,ii-3));
         %   Ux(:,ii) = Ux(:,ii-1) - dt * ( beta / rho ) * Dx * ( 23/12 * S(:,ii-1) - 4/3 * S(:,ii-2) + 5/12*S(:,ii-3) );
         %   Uy(:,ii) = Uy(:,ii-1) - dt * ( beta / rho ) * Dy * ( 23/12 * S(:,ii-1) - 4/3 * S(:,ii-2) + 5/12*S(:,ii-3) );
         %   Uz(:,ii) = Uz(:,ii-1) - dt * ( beta / rho ) * Dz * ( 23/12 * S(:,ii-1) - 4/3 * S(:,ii-2) + 5/12*S(:,ii-3) );
         %end
         case 2
            S0  = S1  - dt * ( Dx * Ux1 + Dy * Uy1 + Dz * Uz1 );
            Ux0 = Ux1 - dt * ( beta / rho ) * Dx * S1;
            Uy0 = Uy1 - dt * ( beta / rho ) * Dy * S1;
            Uz0 = Uz1 - dt * ( beta / rho ) * Dz * S1;
         case 3
            S0 = S1 - dt * Dx * ( 1.5 * Ux1 - 0.5 * Ux2 ) ...
                    - dt * Dy * ( 1.5 * Uy1 - 0.5 * Uy2 ) ...
                    - dt * Dz * ( 1.5 * Uz1 - 0.5 * Uz2 );
            Ux0 = Ux1 - dt * ( beta / rho ) * Dx * ( 1.5 * S1 - 0.5 * S2 );
            Uy0 = Uy1 - dt * ( beta / rho ) * Dy * ( 1.5 * S1 - 0.5 * S2 );
            Uz0 = Uz1 - dt * ( beta / rho ) * Dz * ( 1.5 * S1 - 0.5 * S2 );
         otherwise
            S0 = S1 - dt * Dx * ( (23/12) * Ux1 - 4/3 * Ux2  + 5/12 * Ux3 ) ...
                    - dt * Dy * ( (23/12) * Uy1 - 4/3 * Uy2  + 5/12 * Uy3 ) ...
                    - dt * Dz * ( (23/12) * Uz1 - 4/3 * Uz2  + 5/12 * Uz3 );
            Ux0 = Ux1 - dt * ( beta / rho ) * Dx * ( 23/12 * S1 - 4/3 * S2 + 5/12*S3 );
            Uy0 = Uy1 - dt * ( beta / rho ) * Dy * ( 23/12 * S1 - 4/3 * S2 + 5/12*S3 );
            Uz0 = Uz1 - dt * ( beta / rho ) * Dz * ( 23/12 * S1 - 4/3 * S2 + 5/12*S3 );
         end

         fprintf( [ '|s0| : ' num2str( norm( S0 ) ) '\n' ] );

         % Time advance.
         Ux3 = Ux2;
         Ux2 = Ux1;
         Ux1 = Ux0;
         Uy3 = Uy2;
         Uy2 = Uy1;
         Uy1 = Uy0;
         Uz3 = Uz2;
         Uz2 = Uz1;
         Uz1 = Uz0;
         S3  = S2;
         S2  = S1;
         S1  = S0;

   end

   % Reshape the data for future plotting.
   S0  = reshape( S0,  [nx, ny, nz] );
   Ux0 = reshape( Ux0, [nx, ny, nz] );
   Uy0 = reshape( Uy0, [nx, ny, nz] );
   Uz0 = reshape( Uz0, [nx, ny, nz] );


