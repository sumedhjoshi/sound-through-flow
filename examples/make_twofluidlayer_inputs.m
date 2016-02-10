%
% Make input files for the two-fluid-layer tests.
%
% Sumedh Joshi
% 22 June 2015.

% Assumes that sigma, the Gaussian width, is set to 25.0.

% Add some paths.
addpath( genpath( '/home/smj96/Dropbox/repos/spectral-element-method-acoustics' ) );
addpath( genpath( '/Users/joshi/Dropbox/repos/spectral-element-method-acoustics' ) );

% Set the run-name.
runname = 'twofluidlayer';

% Set some physics constants.
rho1 = 1.02;
c1   = 343.0;
rho2 = rho1;
c2   = 2 * c1;

% Set some parameters.

   % Set the Gaussian width.
   sigma = 25.0;

   % Calculate some wavelength.
   lambda1 = sigma;
   lambda2 = sigma * c2 / c1;

   % Set the height over the interface.
   h = 15.0 * sigma;

   % Set the final time in the simulation as the two-way travel time.
   tfinal = 2 * h / c1;

   % Set the domain sizes.
   cx = 2 * c1 * tfinal;
   L1 = cx + h;
   L2 = c2 / c1 * 0.5 * L1;
   Lz = L1 + L2;
   Ly = 2 * c1 * tfinal;
   Lx = 2 * c1 * tfinal;

   % Set the discretization constants.
   lambda_min = min( lambda1, lambda2 );
   n  = 10;
   mx = ceil( 3 * Lx / ( n * lambda_min ) );
   my = ceil( 3 * Ly / ( n * lambda_min ) );
   mz = 2^nextpow2( ceil( 3 * Lz / ( n * lambda_min ) ) );

% Build the grid.
[x y z] = sem_build_cartesian_mesh( n, mx, my, mz, [-Lx/2, Lx/2], [-Ly/2, Ly/2], [-cx, Lz - cx] );

% Build the initial conditions.
s = exp( -( x.^2 + y.^2 + z.^2 ) / ( sigma.^2 ) );

% Build a smoothly varying background bulk modulus.
rho  = 1.02;
z1 = squeeze(z(1,1,:));
thickness = 20;
beta1 = ( rho * c2^2 - rho * c1^2 ) * 0.5 * ( tanh( ( z1 - 250.0 ) / thickness ) + 1 ) + rho * c1^2;

% Set the background parameters.
beta = zeros( n * mx, n * my, n * mz );
for ii = 1: n * mx
  for jj = 1: n * my
     beta( ii, jj, : ) = beta1;
  end
end

% Write the initial conditions file.
sem_write_initfile( n, mx, my, mz, ...
                    x, y, z, ...
                    s, 0, 0, 0, ...
                    rho, beta, ...
                    0, 0, 0, ...
                    [ pwd '/' runname '_init.h5' ] );

% Now build the inputs.
input.fname_runname = runname;
input.fname_init    = [runname '_init.h5'];
input.n  = n;
input.mx = mx;
input.my = my;
input.mz = mz;
input.Lx = [ -Lx/2, Lx/2 ];
input.Ly = [ -Ly/2, Ly/2 ];
input.Lz = [ -cx, Lz - cx ];
input.mu = 0.0;
input.dt = 1e3;
input.t_final = tfinal;
input.report_every_n_steps = 1;
input.timesteps_between_writes = 100000; % Don't write until the end.

% Write the input file.
sem_write_inputfile( [ runname '_in'], input );
