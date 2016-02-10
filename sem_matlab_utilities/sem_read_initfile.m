function [data] = sem_read_initfile( init_file_name )
% data = sem_read_initfile( init_file_name )
%
% Read an initial conditions file and returns its contents.
%
% Takes 1 argument:
%
%   init_file_name - String indicating the initial conditions file to read
%                    from disk.
%
% Returns 1 value:
%
%   data        - Struct with fields specified as below.
%     .grid     - Struct with fields specified as below.
%       n           - Number of GLL points per direction, per subdomain.
%       mx          - Number of subdomains in the x-direction.
%       my          - Number of subdomains in the y-direction.
%       mz          - Number of subdomains in the z-direction.
%       x           - Matrix, of size mx*n by my*n by mz*n, containing the
%                     x-coordinates of the mesh associated with the field.
%       y           - Matrix, of size mx*n by my*n by mz*n, containing the
%                     z-coordinates of the mesh associated with the field.
%       z           - Matrix, of size mx*n by my*n by mz*n, containing the
%                     z-coordinates of the mesh associated with the field.
%     .ic       - Struct with fields specified as below.
%        s0          - Matrix, of dim mx * n by my * n by mz * n containing the
%                      initial density at each grid point.
%        ux0         - Matrix, of dim mx * n by my * n by mz * n containing the
%                      initial x-velocity at each grid point.
%        uy0         - Matrix, of dim mx * n by my * n by mz * n containing the
%                      initial z-velocity at each grid point.
%        uz0         - Matrix, of dim mx * n by my * n by mz * n containing the
%                      initial z-velocity at each grid point.
%     .environment - Struct containing the flow.  All arrays are matrices of size
%                      mx*n by my*n by mz*n.
%
%        rho         - background density.
%        beta        - fluid bulk modulus.
%        vx          - background velocity in x.
%        vy          - background velocity in y.
%        vz          - background velocity in z.
%
% 23 Jun 2013
% Sumedh Joshi

    % Read the grid information and the grid.
    data.grid.n  = h5read( init_file_name, '/grid/n' );
    data.grid.mx = h5read( init_file_name, '/grid/mx' );
    data.grid.my = h5read( init_file_name, '/grid/my' );
    data.grid.mz = h5read( init_file_name, '/grid/mz' );
    data.grid.x  = h5read( init_file_name, '/grid/x' );
    data.grid.y  = h5read( init_file_name, '/grid/y' );
    data.grid.z  = h5read( init_file_name, '/grid/z' );

    % Read the initial conditions.
    data.ic.ux = h5read( init_file_name, '/ic/ux' );
    data.ic.uy = h5read( init_file_name, '/ic/uy' );
    data.ic.uz = h5read( init_file_name, '/ic/uz' );
    data.ic.s  = h5read( init_file_name, '/ic/s' );

    % Read the environment.
    data.environment.rho  = h5read( init_file_name, '/environment/rho' );
    data.environment.beta = h5read( init_file_name, '/environment/beta' );
    data.environment.vx   = h5read( init_file_name, '/environment/vx' );
    data.environment.vy   = h5read( init_file_name, '/environment/vy' );
    data.environment.vz   = h5read( init_file_name, '/environment/vz' );

end
