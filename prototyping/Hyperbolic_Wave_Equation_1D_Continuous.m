%
% Solve the 1D homogenous wave equation in split density/velocity form
% as a coupled set of hyperbolic conservation laws and a spectral
% element discretization.

   % Set some constants.
   rho  = 1.02;
   c    = 343;
   beta = rho * c^2;

   % Set some grid parameters.
   L  = 1000;
   n  = 16;
   mx = 32;

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
   tmax = 1.0;
   dx   = diff( x );
   dx   = min( dx( dx > 0 ) );
   dt   = (1/10) * dx / c;
   Nt   = ceil( tmax / dt );
   t    = linspace( 0, tmax, Nt );

   % Compute the spectral differentiation matrix of order 1.
   d1 = zeros( n, n );
   for ii = 1: n
      d1(:, ii) = make_lagrange( eta', eta', ii, 1 );
   end

   % Tile the spectral differentiation matrix to build a multi-domain grid.
   D1 = ( 2 / hx ) * kron( eye(mx,mx), d1 );

   % Update the interfacial boundary conditions to be identical.
   for ii = 1:mx - 1
      left  = ( ii - 1 ) * n + n;
      right = left + 1;
      ndx = n * ii;
      interface = ( D1( left:right, : ) + flipud( D1( left:right, :) ) );
      interface( 1:2, left:right ) =  interface( 1:2, left:right) + 0 * ( 2 / hx )^2 * [ 1.0, 1.0; 1.0, 1.0];
      %interface = D1( left:right, : ) + flipud( D1( left:right, :) );
      %interface = [ 1.0, -1.0; 1.0, -1.0 ]; % Weak C0 continuity.
      D1( left:right, : ) =  interface/2;
   end

   % Update the global matrix for Dirichlet boundary conditions.
   D1(1,1)     =  2/hx * 0.5;
   D1(end,end) = -2/hx * 0.5;

   % Consistency check.
   u = sin( pi*x/100 );
   Du_a = pi/100*cos( pi*x/100 );
   err = norm( Du_a - D1*u );

%   % Assemble the 3D 1st derivative multidomain matrix.
%   mz = mx;
%   my = mx;
%   nx = n * mx;
%   ny = n * my;
%   nz = n * mz;
%   % Dx = Iz \kron Iy \kron D1;
%   % Dy = Iz \kron D1 \kron Ix;
%   % Dz = D1 \kron Iy \kron Ix;
%   Dx = kron( speye( nz * ny ), sparse(D1) );
%   Dy = kron( kron( speye( nz ), sparse(D1) ), speye(nx) );
%   Dz = kron( sparse(D1), speye( nx*ny ) );
%
%   % Test 3D consistency.
%   z = x;
%   y = x;
%   [X Y Z] = meshgrid( x, y, z );
%   X = X(:);
%   Y = Y(:);
%   Z = Z(:);
%   U = sin( pi * X / 100 ) .* sin( pi * Y / 100 ) .* sin( pi * Z / 100 );
%   X = reshape( X, [nx, ny, nz] );
%   Y = reshape( Y, [nx, ny, nz] );
%   Z = reshape( Z, [nx, ny, nz] );
%   LU = Dx * U + Dy * U + Dz * U;
%   LU_a = cos( pi * X / 100 ) .* sin( pi * Y / 100 ) .* sin( pi * Z / 100 ) * pi / 100 + ...
%          sin( pi * X / 100 ) .* cos( pi * Y / 100 ) .* sin( pi * Z / 100 ) * pi / 100 + ...
%          sin( pi * X / 100 ) .* sin( pi * Y / 100 ) .* cos( pi * Z / 100 ) * pi / 100;
%   LU = reshape( LU , [nx, ny, nz] );
%   LU_a = reshape( LU_a , [nx, ny, nz] );

   % Set the initial condition.
   sigma = 25;
   x0 = L/2;
   s0 = exp( -(x - x0).^2/(2 * sigma^2 ) );
   u0 = 0 * s0;

   % Compute the dt based on frequency of the initial pulse.
   lambda    = 2 * sigma;
   frequency = c / lambda;
   period    = 1 / frequency;

   % Build the operator matrix for eigenvalue identification.
   r = mx * n;
   A = eye( 2*r, 2*r );
   A( 1:r, r+1:2*r) = -dt * D1;
   A( r+1:2*r, 1:r) = -dt * beta * D1 / rho;

   % Time-propagate.
   U = zeros( n * mx, length(t) );
   S = U;
   U(:,1) = u0;
   S(:,1) = s0;
   for ii = 2:length(t)
      p( ii, Nt, 1);

      % Depending on the time step, use a particular ABk method.
      switch ii
         case 2
            S(:,ii) = S(:,ii-1) - dt * D1 * U(:,ii-1);
            U(:,ii) = U(:,ii-1) - dt * beta * D1 * S(:,ii-1) / rho;
         case 3
            S(:,ii) = S(:,ii-1) - dt * D1 * ( 1.5 * U(:,ii-1) - 0.5 * U(:,ii-2) );
            U(:,ii) = U(:,ii-1) - dt * ( beta / rho ) * D1 * ( 1.5 * S(:,ii-1) - 0.5 * S(:,ii-2) );
         otherwise
            S(:,ii) = S(:,ii-1) - dt * D1 * ( (23/12) * U(:,ii-1) - 4/3 * U(:,ii-2)  + 5/12 * U(:,ii-3));
            U(:,ii) = U(:,ii-1) - dt * ( beta / rho ) * D1 * ( 23/12 * S(:,ii-1) - 4/3 * S(:,ii-2) + 5/12*S(:,ii-3) );
         end

   end
