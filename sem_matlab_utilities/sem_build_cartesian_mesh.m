function [x, y, z] = sem_build_cartesian_mesh( n, mx, my, mz, xlims, ylims, zlims )
% [x, y, z] = sem_build_cartesian_mesh( n, mx, my, mz, xlims, ylims, zlims )
%
% Generate a cartesian (undeformed) mesh in xi-first indexing for reading by
% the SMPM solver code.
%
% Takes 5 arguments:
%
%   n              - Number of GLL points per direction, per subdomain.
%   mx             - Number of subdomains in the x-direction.
%   mz             - Number of subdomains in the z-direction.
%   xlims          - Row vector, of length 1 x 2, specifying bounds for the
%                    x-dimension.
%   ylims          - Row vector, of length 1 x 2, specifying bounds for the
%                    y-dimension.
%   zlims          - Row vector, of length 1 x 2, specifying bounds for the
%                    z-dimension.
%
% Returns 2 values:
%
%   (x,y,z) - 3D arrays representing the grid.
%
% 23 June 2013
% Sumedh Joshi

    % Get the GLL points.
    xi  = lglnodes( n - 1 );
    xi  = sort( xi );
    eta = xi;
    chi = xi;

    % Build the 1D multidomain mesh in z.
    Lz     = (zlims(2) - zlims(1)) / mz;
    xi     = (xi + 1) * Lz / 2 - (zlims(2) - zlims(1)) / 2;
    XI     = repmat( xi, mz, 1 );
    offset = repmat( [0:mz-1]' * Lz, 1, n )';
    offset = offset(:);
    XI     = XI + offset;

    % Build the 1D multidomain mesh in y.
    Ly     = (ylims(2) - ylims(1)) / my;
    chi    = (chi + 1) * Ly / 2 - (ylims(2) - ylims(1)) / 2;
    CHI    = repmat( chi, my, 1 );
    offset = repmat( [0:my-1]' * Ly, 1, n )';
    offset = offset(:);
    CHI    = CHI + offset;

    % Build the 1D multidomain mesh in x.
    Lx     = (xlims(2) - xlims(1)) / mx;
    eta    = (eta + 1) * Lx / 2 - (xlims(2) - xlims(1)) / 2;
    ETA    = repmat( eta, mx, 1 );
    offset = repmat( [0:mx-1]' * Lx, 1, n )';
    offset = offset(:);
    ETA    = ETA + offset;

    % Set the range.
    ETA = ETA - min(ETA) + xlims(1);
    CHI = CHI - min(CHI) + ylims(1);
    XI  = XI  - min(XI)  + zlims(1);

    % Do the meshgrid.
    [x y z] = ndgrid( ETA, CHI, XI );

end
