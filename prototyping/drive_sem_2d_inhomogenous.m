

% Set some constants.
n = 10;
m = 32;
L = 500;
x = linspace( -L, L, n );

% Build a 1D sem cartesian mesh.
mx = m;
my = m ;
[x y z] = sem_build_cartesian_mesh( n, mx, my, 1, [ -L, L ], [ -L, L ], [ -L, L ] );
X = squeeze( x(1,:,:) );
Y = squeeze( y(1,:,:) );
clear z;
x = X(1,:);
y = Y(:,1);

% Set the initial condition.
sigma = 10;
s0    = exp( - ( X.^2 + Y.^2 ) / sigma^2 );
ux0   = 0 * X;
uy0   = 0 * Y;

% Set some physical constants.
c0   = 343;
rho0 = 1.0;
vx0  = 100;
vy0  = 0;

% Set the time discretization.
dt = ( Y(2) - Y(1) ) / c0 / 25;
tf = 0.250;
t  = linspace( 0, tf, round( tf / dt ) );

% Solve the problem.
data = sem_2D_inhomogenous( n, m, m, X, Y, t, s0, ux0, uy0, rho0, vx0, vy0, c0 );

% Plot the solution with arrival time lines.
ndx = round( length(t) * [ 0.25, 0.5, 0.75, 1 ] );
for ii = 1:4
   subplot( 1, 4, ii );
   contourf( data.x, data.y, data.S(:,:,ndx(ii) ) );
   xlabel( 'x - m' );
   ylabel( 'y - m' );
   title( [ 't = ', num2str( t(ndx(ii)) ), ' seconds ']);
   if ii == 1
      ca = caxis;
   else
      caxis( ca );
   end
end
