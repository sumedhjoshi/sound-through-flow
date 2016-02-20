%
% Make input files for the Gaussian-source in a homogenous fluid test.
%
% Sumedh Joshi
% 23 June 2015.

% Add some paths.
addpath( genpath( '/home/smj96/Dropbox/repos/sound-through-flow' ) );
addpath( genpath( '/Users/joshi/Dropbox/repos/sound-through-flow' ) );

% Set the run-name.
runname = 'crossflow';

% Set some physics constants.
rho = 1.02;
c   = 343.0;

% Set some parameters.

   % Set the Gaussian width.
   sigma = 10.0;

   % Set the domain size.
   Lx = [-250, 250];
   Ly = [-250, 250];
   Lz = [-250, 250];

   % Set the final time in the simulation.
   tfinal = 0.25;

   % Set the discretization constants.
   n  = 10;
   mx = 8;
   my = 8;
   mz = 8;

% Build the grid.
[x y z] = sem_build_cartesian_mesh( n, mx, my, mz, Lx, Ly, Lz );

% Build the initial conditions.
s = exp( -( x.^2 + y.^2 + z.^2 ) / ( sigma.^2 ) );

% Set the background parameters.
rho  = 1.02;
beta = rho * c^2;

% Set the background velocities.
vx = 300.0;
vy = 0.0;
vz = 0.0;

% Write the initial conditions file.
sem_write_initfile( n, mx, my, mz, ...
                    x, y, z, ...
                    s, 0, 0, 0, ...
                    rho, beta, ...
                    vx, vy, vz, ...
                    [ pwd '/' runname '_init.h5' ] );

% Set the time-step based on the CFL condition.
x1d = unique( x(:,1,1) );
dx = min( diff( x1d ) );
dt = dx / c;

% Now build the inputs.
input.fname_runname = runname;
input.fname_init    = [runname '_init.h5'];
input.n  = n;
input.mx = mx;
input.my = my;
input.mz = mz;
input.Lx = Lx;
input.Ly = Ly;
input.Lz = Lz;
input.mu = 0.0;
input.dt = dt / 25.0;
input.t_final = tfinal;
input.report_every_n_steps = 1;
input.timesteps_between_writes = 100; % Don't write until the end.

% Write the input file.
sem_write_inputfile( [ runname '_in'], input );
