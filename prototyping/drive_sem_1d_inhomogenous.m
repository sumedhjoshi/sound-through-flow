

% Set some constants.
n = 16;
m = 16;
L = 500;
x = linspace( -L, L, n );

% Build a 1D sem cartesian mesh.
[x junk junk] = sem_build_cartesian_mesh( n, m, m, m, [ -L, L ], [ -L, L ], [ -L, L ] );
x = x(:,1,1);

% Set the initial condition.
sigma = 10;
s0    = exp( - x.^2 / sigma^2 );
u0    = 0 * x;

% Set some physical constants.
c0   = 343;
rho0 = 1.0;
v0   = 0;

% Set the time discretization.
dt = ( x(2) - x(1) ) / c0 / 20;
tf = 0.5;
t  = linspace( 0, tf, round( tf / dt ) );

% Solve the problem with finite differences.
data = sem_1D_inhomogenous( n, m, x, t, s0, u0, rho0, v0, c0 );

% Plot the solution with arrival time lines.
plot( data.x, data.S(:,end) );
vline( t(end) * ( mean( v0 ) + c0 ), 'r' );
vline( - t(end) * ( c0 - mean( v0 ) ), 'r' );
