

% Set some constants.
n = 256;
L = 250;
x = linspace( -L, L, n );
y = linspace( -L, L, n );
[X Y] = meshgrid( x, y );
X = X(:);
Y = Y(:);

% Set the initial condition.
sigma = 10;
s0    = exp( 1 - ( X.^2 + Y.^2 ) / sigma^2 );
ux0   = 0 * X;
uy0   = 0 * Y;

% Set some physical constants.
c0   = 343;
rho0 = 1.0;
vx0  = 0;
vy0  = 0;

% Set the time discretization.
dt = ( x(2) - x(1) ) / c0 / 20;
tf = 0.250;
t  = linspace( 0, tf, round( tf / dt ) );

% Solve the problem.
data = fdtd_2D_inhomogenous( X, Y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 );

% Plot the solution with arrival time lines.
imagesc( x, y, reshape( data.S(:,end), n, n ) );
