%%
%
% 2D finite difference time domain simulation of the inhomogenous wave
% equation computed from the hyperbolic conservation laws derived
% previously, for the case of a Rankine vortex tornado. - 3/18/2013

% Loop over some parameters.
GammaVec = [-12500, -10000, -5000, -2500, 0, 2500, 5000, 10000, 12500];
GammaVec = [-12000, -12000, -12000];
GammaVec = [-12000];
GammaVec = [-10000];
SigmaVec = [32];

% Movie snapshot every n timesteps.
nsteps_snapshot = 5000;

for jj = 1:length(GammaVec)

    fprintf([' Running tornado model for Gamma value ', num2str(GammaVec(jj)) '\n \n \n']);
    keep GammaVec jj nsteps_snapshot SigmaVec

    % Set some grid parameters.
    n    = 501;
    xlim = [-1600 1600];
    ylim = [-1600 1600];

    % Set some new parameters, just for the data visualization.
    n    = 801;
    xlim = [-1600 1600];
    ylim = [-1600 1600];

    % Set some new parameters, just for the data visualization.
    n    = 501;
    xlim = [-1600 1600];
    ylim = [-1600 1600];

    % Set the bulk modulus.
    c    = 343;
    rho0 = 1.02;
    beta = rho0 * c^2;

    % Build the grid.
    x = linspace(xlim(1),xlim(2),n)';
    y = linspace(ylim(1),ylim(2),n)';
    [X Y] = meshgrid(y,x);
    X = X(:);
    Y = Y(:);
    hx = x(2) - x(1);
    hy = y(2) - y(1);
    Lx = hx * n;
    Ly = hy * n;

    % Setup some arrays to hold the snapshots.
    Ssnap = zeros(n*n,16);
    Ssnapdiff = [];
    snapcount = 1;

    % Set some time stepping parameters.
    dt   = hx/c/20;
    tmin = 0.0;
    tmax = (Lx / 2) / c;
    t    = linspace(tmin,tmax,round((tmax-tmin)/dt));
    nt   = length(t);

    % Set some constants.
    sigma = SigmaVec(jj);

    % Set some initial conditions.
    u0x = zeros(n*n,1);
    u0y = zeros(n*n,1);
    x0 = -0;
    y0 = 0;
    s0 = exp( - ((X - x0).^2 + (Y - y0).^2 ) / sigma^2);

    % Set some constants for the Rankine vortex.
    tornado_center = max(abs(X)) - 200;
    tornado_center = 800;
    two_way_travel_time = 2*tornado_center/c;
    Gamma = GammaVec(jj);
    Rv    = 25;
    nu    = 10e-6;
    nu    = 0;
    tornado_velocity = 0.0; % Meters per second ==> 70 miles per hour.
    x0 = tornado_center;
    theta = atan(Y./(X-x0));

        % Sometimes if a grid point coincides with the center of the tornado,
        % theta will be undefined.  Fix this.
        nan_ndx = find(isnan(theta));
        theta(nan_ndx) = 0;

    R     = sqrt((X-x0).^2 + Y.^2);
    U1    = Gamma.*R/(2*pi*Rv^2);
    U2    = Gamma./(2*pi*R);
    U(R<Rv) = U1(R<Rv);
    U(R>=Rv) = U2(R>=Rv);

    T = sign(X - x0);
    T(T == 0) = 1;
    Vx = sin(theta).*U'.*T; %
    Vy = cos(theta).*U'.*T; % Tornado velocity.
    Vx0 = Vx;
    Vy0 = Vy;
    Vx1 = Vx;
    Vy1 = Vy;

    % Call the finite difference model.
    data = fdtd_2D_inhomogenous( X, Y, t, s0, u0x, u0y, rho0, Vx, Vy, c );

end
