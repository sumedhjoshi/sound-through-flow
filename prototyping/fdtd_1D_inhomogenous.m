% data = fdtd_1D_inhomogenous( x, t, s0, u0, rho0, v0, c0 );
function data = fdtd_1D_inhomogenous( x, t, s0, u0, rho0, v0, c0 )

   % Get some grid constants.
   h = x(2) - x(1);
   n = length(x);

   % Set some time stepping parameters.
   dt   = t(2) - t(1);
   tmax = t(end);

   % Set the sound speed and background flow.
   c = c0;
   V = v0;
   beta = rho0 * c0.^2;

   % Build the differentiation matrices.

    % First order central difference.
    D1 = diag(0*ones(n,1)) + diag(1*ones(n-1,1),1) - diag(1*ones(n-1,1),-1);
    D1 = (1/(2*h))*D1;
    sD1  = sparse(D1);

    % Build Dirichlet-Dirichlet and Neumann-Neumann versions of these
    % operators.
    sD1d = sD1;
    sD1n = sD1;
    sD1n(1,1) = 1/h;
    sD1n(1,2) = -1/h;
    sD1n(end,end) = 1/h;
    sD1n(end,end-1) = -1/h;

   % Build the time-advancment operator.
   A = sparse( [ speye(n,n) - dt * diag(v0) * sD1d, -dt * sD1n; ...
                -dt * diag(c0).^2 * sD1d, speye(n,n) - dt * diag(v0) .* sD1n ] );

   % Assemble the data and solve.
   U = zeros(n,length(t));
   S = zeros(n,length(t));
   S = zeros(n,length(t));
   U(:,1) = u0;
   S(:,1) = s0;
   for ii = 2:length(t)
     iiphi = A * [ S(:,ii-1); U(:,ii-1) ];
     S(:,ii) = iiphi(1:n);
     U(:,ii) = iiphi(n+1:end);
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
