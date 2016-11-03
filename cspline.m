function [xx,yy] = cspline(x,y,npts,varargin)
%
% [xx,yy] = cspline(x,y,npts)
% [xx,yy] = cspline(x,y,npts,stype,<type>)
% [xx,yy] = cspline(x,y,npts,pscale,<scale>)
% [xx,yy] = cspline(x,y,npts,<Property>,<Value>,...)
%
% Reconstruct curve of generic set of points and any shape, on a substrate
% of "npts" points. DO NOT EXTRAPOLATE.
%
% Input arguments:
% - x        : x-data array of size 1-by-n
% - y        : y-data array of same size of x or size n-by-k
% - npts     : number of points to consider in the (x,y) domain.
% - varargin :   Use 'option',<val> for optional input arguments (see ProduceCorrectVarargin).
%   Accepted 'option' strings are:
%   + stype  : {'pchip'} | 'spline'   Spline interpolator type.
%   + xscale, yscale: {'linear'} or {'lin'} | 'log'  Use linear or
%                     logarithmic spacing for the npts points to
%                     interpolate.
%
% see also PCHIP, SPLINE, PRODUCECORRECTVARARGIN.
% 
% v1.1
% Revised and added help.
% Maurizio De Pitta', The University of Chicago, Chicago, April 28th, 2016.
%
% v1.0
% Maurizio De Pitta', Tel Aviv University, Tel Aviv, Israel, January 28th, 2012.
% 
% https://sites.google.com/site/mauriziodepitta/home
% maurizio.depitta@gmail.com

%--------------------------------------------------------------------------
% Defaults
%--------------------------------------------------------------------------
opts.xscale = 'lin';
opts.yscale = 'lin';
opts.stype = 'pchip';

%--------------------------------------------------------------------------
% User-defined values
%--------------------------------------------------------------------------
if ~isempty(varargin)
    varargin = ProduceCorrectVarargin(varargin);
    for i = 1:length(varargin)/2
        if isfield(opts,varargin{2*i-1})
            opts.(genvarname(varargin{2*i-1})) = varargin{2*i};
        end
    end
end

%--------------------------------------------------------------------------
% Initialize auxiliary substrate
%--------------------------------------------------------------------------
% Auxiliary variable
t = 1:length(x);

% Define interpolant substrate (X-scale)
switch lower(opts.xscale)
    case {'linear','lin'}
        ttx = linspace(t(1),t(end),npts);
    case 'log'
        ttx = logspace(log10(t(1)),log10(t(end)),npts);
end

% Define interpolant substrate (Y-scale)
switch lower(opts.yscale)
    case {'linear','lin'}
        tty = linspace(t(1),t(end),npts);
    case 'log'
        tty = logspace(log10(t(1)),log10(t(end)),npts);
end

%--------------------------------------------------------------------------
% Effective spline interpolation
%--------------------------------------------------------------------------
switch opts.stype
    case 'pchip'
        % pchip length(x)==size(y,1) for correct interpolation, i.e. each trial is
        % given in columns
        xx = pchip(t,x,ttx);
        yy = pchip(t,y',tty);
    case 'spline'
        % spline length(x)==size(y,2) for correct interpolation, i.e. each trial is
        % given in rows
        xx = spline(t,x,ttx);
        % flip y data
        yy = spline(t,y',tty)';
end

% Restore format
xx = xx';
yy = yy';