% data = fdtd_2D_inhomogenous( x, y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 );
function data = fdtd_2D_inhomogenous( x, y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 )

   % Get some grid constants.
   n  = sqrt( length( x ) );
   hx = ( max( x ) - min( x ) ) / n;
   hy = ( max( y ) - min( y ) ) / n;

   % Set some time stepping parameters.
   dt   = t(2) - t(1);
   tmax = t(end);

   % Set some grid parameters.
   xlim = [min(x) max(x)];
   ylim = [min(y) max(y)];

   % Set the bulk modulus.
   c    = 343;
   rho0 = 1.02;
   beta = rho0 * c^2;

   % Set the domain lengths.
   Lx = hx * n;
   Ly = hy * n;

   % Set some time stepping parameters.
   dt   = t(2) - t(1);
   tmin = t(1);
   tmax = t(end);

   % Set some initial conditions.
   u0x = ux0;
   u0y = uy0;

    % The hyperbolic conservation laws are:
    %
    %       1. dsdt = - div( u )
    %       2. dudt = - beta * grad(s) / rho0
    %
    % Translated to finite difference approximations, they are:
    %
    %       1. s(k+1) = s(k) - dt * div( u(k) )
    %       2. u(k+1) = u(k) - dt * beta * grad( s(k) ) / rho0
    %
    % Where beta is the fluid bulk modulus.

    % Build the first order central difference.
    D1  = diag(0*ones(n,1)) + diag(1*ones(n-1,1),1) - diag(1*ones(n-1,1),-1);
    D1x = sparse(kron(sparse(eye(n)),D1)*(1/(2*hx)));
    D1y = sparse(kron(D1,sparse(eye(n)))*(1/(2*hy)));

   % Set the background field.
   Vx = vy0;
   Vy = vx0;

   % Start time stepping.
   Ux = zeros(n^2,length(t));
   Uy = zeros(n^2,length(t));
   S  = zeros(n^2,length(t));
   Ux(:,1) = u0x;
   Uy(:,1) = u0y;
   S(:,1)  = s0;
   for ii = 2:length(t)

      p(ii,length(t),1);

      % Build the condensation.
      S(:,ii)  =  S(:,ii-1) - dt * ( D1x * Ux(:,ii-1) + D1y * Uy(:,ii-1) );

      % Build the velocity in x.
      Ux(:,ii) = Ux(:,ii-1) - dt * beta * D1x * S(:,ii-1) / rho0  ...
                            - dt * Vx(:, 1) .* (D1x * Ux(:,ii-1)) ...
                            - dt * Vy(:, 1) .* (D1y * Ux(:,ii-1)) ...
                            - dt * ( (D1x * S(:,ii-1)).*Vx(:,1) + (D1x * S(:,ii-1)).*Vx(:,1) ).*Vx(:,1);

      % Build the velocity in y.
      Uy(:,ii) = Uy(:,ii-1) - dt * beta * D1y * S(:,ii-1) / rho0   ...
                            - dt * Vx(:, 1) .* ( D1x * Uy(:,ii-1)) ...
                            - dt * Vy(:, 1) .* ( D1y * Uy(:,ii-1)) ...
                            - dt * ( (D1x * S(:,ii-1)).*Vx(:,1) + (D1y * S(:,ii-1)).*Vy(:,1) ).*Vy(:,1);

   end

   % Pack the data for output.
   data.S  = S;
   data.Ux = Ux;
   data.Uy = Uy;
   data.t  = t;
   data.x  = x;
   data.y  = y;
   data.physics.c0   = c0;
   data.physics.Vx   = Vx;
   data.physics.Vy   = Vy;
   data.physics.rho0 = rho0;

end




