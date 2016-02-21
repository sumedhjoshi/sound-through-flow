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
   if length(v0) == 1
      V = v0 * ones( length(x), 1);
   else
      V = v0;
   end
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
   %A = sparse( [ speye(n,n) - dt * diag(v0) * sD1d, -dt * sD1n; ...
   %             -dt * diag(c0).^2 * sD1d, speye(n,n) - dt * diag(v0) .* sD1n ] );

   % Assemble the data and solve.
%   U = zeros(n,length(t));
%   S = zeros(n,length(t));
%   S = zeros(n,length(t));
%   U(:,1) = u0;
%   S(:,1) = s0;
%   for ii = 2:length(t)
%     iiphi = A * [ S(:,ii-1); U(:,ii-1) ];
%     S(:,ii) = iiphi(1:n);
%     U(:,ii) = iiphi(n+1:end);
%   end
   U = zeros(n,length(t));
   S = zeros(n,length(t));
   U(:,1) = u0;
   S(:,1) = s0;
   M0 = 0 * u0;
   M1 = 0 * u0;
   M2 = 0 * u0;
   C0 = 0 * u0;
   C1 = 0 * u0;
   C2 = 0 * u0;
   for ii = 2:length(t)

      p( ii, length(t), 1);

      % Compute the continuity and momentum updates.
      C0 = ( V .* ( sD1d * S(:,ii-1) ) + sD1n * U(:,ii-1) );
      M0 = ( V .* ( sD1n * U(:,ii-1) ) + U(:,ii-1) .* ( sD1n * V ) + S(:,ii-1) .* V .* ( sD1n * V ) +  beta * sD1d * S(:,ii-1) / rho0 );

      % Depending on the time step, use a particular ABk method.
      switch ii
         case 2
            S(:,ii) = S(:,ii-1) - dt * C0;
            U(:,ii) = U(:,ii-1) - dt * M0;
            %S(:,ii) = S(:,ii-1) - dt * ( V .* ( D1 * S(:,ii-1) ) + D1 * U(:,ii-1) );
            %U(:,ii) = U(:,ii-1) - dt * ( V .* ( D1 * U(:,ii-1) ) +  +  beta * D1 * S(:,ii-1) / rho );
         case 3
            S(:,ii) = S(:,ii-1) - dt * ( 1.5 * C0 - 0.5 * C1 );
            U(:,ii) = U(:,ii-1) - dt * ( 1.5 * M0 - 0.5 * M1 );
            %S(:,ii) = S(:,ii-1) - dt * D1 * ( 1.5 * U(:,ii-1) - 0.5 * U(:,ii-2) );
            %U(:,ii) = U(:,ii-1) - dt * ( beta / rho ) * D1 * ( 1.5 * S(:,ii-1) - 0.5 * S(:,ii-2) );
         otherwise
            S(:,ii) = S(:,ii-1) - dt * ( 23/12 * C0 - 4/3 * C1 + 5/12 * C2 );
            U(:,ii) = U(:,ii-1) - dt * ( 23/12 * M0 - 4/3 * M1 + 5/12 * M2 );
            %S(:,ii) = S(:,ii-1) - dt * D1 * ( (23/12) * U(:,ii-1) - 4/3 * U(:,ii-2)  + 5/12 * U(:,ii-3));
            %U(:,ii) = U(:,ii-1) - dt * ( beta / rho ) * D1 * ( 23/12 * S(:,ii-1) - 4/3 * S(:,ii-2) + 5/12*S(:,ii-3) );
      end

      % Time advance the lagging terms.
      M2 = M1;
      C2 = C1;
      M1 = M0;
      C1 = C0;
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
