%
% data = sem_1D_inhomogenous( n, mx, x, t, s0, u0, rho0, v0, c0 );
function data = fdtd_1D_inhomogenous( n, mx, x, t, s0, u0, rho0, v0, c0 )

   % Set some constants.
   rho  = 1.02;
   c    = 343;
   beta = rho * c^2;

   % Set some grid parameters.
   L  = x(end) - x(1);

   % Set the sound speed and background flow.
   c = c0;
   V = v0;
   beta = rho0 * c0.^2;

   % Get the grid element width.
   hx = L / mx;

   % Build the GLL grid in one element.
   [eta junk junk] = lglnodes( n - 1 );
   eta = -eta;

   % Set the time discretization.
   tmax = 1.0;
   dx   = diff( x );
   dx   = min( dx( dx > 0 ) );
   Nt   = length(t);
   dt   = t(2) - t(1);

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
      D1( left:right, : ) =  interface/2;
   end

   % Update the global matrix for Dirichlet boundary conditions.
   D1(1,1)     =  2/hx * 0.5;
   D1(end,end) = -2/hx * 0.5;

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

   % Pack the data for output.
   data.S = S;
   data.U = U;
   data.t = t;
   data.x = x;
   data.physics.c0   = c0;
   data.physics.v0   = v0;
   data.physics.rho0 = rho0;

end
