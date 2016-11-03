function fh = niceplot(fh,varargin)
%
% fh = niceplot(fh,varargin)
%
% Consider a generic figure with data within (and interpolate data) and
% replots them to have equispaced markers or points.
%
% Input arguments:
% - fh       : figure handle where data have ALREADY been plotted but need 
%              to be 'nicely' re-plotted.
% - varargin :   Use 'option',<val> for optional input arguments (see ProduceCorrectVarargin).
%   Accepted 'option' strings are:
%   + minfitpts : {2} | Integer   Minimal set of points to perform cspline interpolation
%   + dense     : {2} | Integer   Considers 'dense' times the number of 
%                 points originally computed by AUTO in each branch to
%                 perform cspline interpolation in the (x,y) domains.
%   + step      : 2e-2| Double    Space between points along an interpolated branch.
%   + special   : {[]}| Valued Array  Specify which branches from the
%                 plotted bifurcation diagram are processed. The array must
%                 be in the form of pairs like [graphobj1, step1,... graphobjn,stepn]
%                 Default: branches are treated in the same fashion.
%
% Returns:
% - fh        : figure handle of modified figure.
%
% see also CSPLINE, PRODUCECORRECTVARARGIN.
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
opts.minfitpts = 2;     % minimal set of points to perform spline
opts.dense = 2;         % multiplicative factor: considers twice the number of points within each branch for spline fitting
opts.step = 2e-2;       % [n.u.] % space between points
opts.special = [];      % No branch-specific special settings (i.e. all data sets are treated in the same fashion)
ispec = 2;              % Initial index of special branch settings

%--------------------------------------------------------------------------
% User-defined options
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
% Retrieve info from figure
%--------------------------------------------------------------------------
% Get Children of figure
% Assumes that your data are of line type and marked also by markers
obj = findobj(fh,'Type','line','-not','Marker','none');
ax = findobj(fh,'Type','axes');

% Retrieve Axis limits
XLim = get(ax,'XLim');
YLim = get(ax,'YLim');

% Retrieve axis scaling
IsXLog = strcmp(get(ax,'XScale'),'log');
IsYLog = strcmp(get(ax,'YScale'),'log');

% First retrieve min-max ranges or plotted data
[xrange,yrange] = GetDataRanges(fh);

% Before processing data, sort opts.special in ascending order if specified
if ~isempty(opts.special)
    % Reshape special for suitable manipulation
    aux = reshape(opts.special,2,length(opts.special)/2)';
    % Order graphic objects accoridng to index in the first column
    aux = sortrows(aux,1);
    % Reshape back to opts.special
    opts.special = reshape(aux',1,length(opts.special));
end

%--------------------------------------------------------------------------
% Effective resizing
%--------------------------------------------------------------------------
for i = 1:length(obj)
    % First retrieve parameter values depending on specific branches that
    % need specific tuning
    stepdist = opts.step;
    % ispec starts from 1 so this condition is verified automatically only
    % is opts.special is set
    if ispec<=size(opts.special,2)
        if i==opts.special(ispec-1)
            stepdist = opts.special(ispec);
            ispec = ispec + 2;
        end
    end
    
    % Get data points
    xdata = get(obj(i),'XData')';
    ydata = get(obj(i),'YData')';
    
    % First check axis scaling and convert accordingly
    [xdata,ydata] = loglinConversion(xdata,ydata,IsXLog,IsYLog,1);
    
    % Normalize data
    [xdata,ydata] = RescaleData(xdata,ydata,xrange,yrange,1);

    % Spline fitting od data on regular substrate (at least two points)
    % Needed for nice computation of distances
    % WARNING: it will probably give an error if data points are all equal
    % (should be improved in the future)
    [xdata,ydata] = cspline(xdata,ydata,max([opts.dense*floor(size(xdata,1)),opts.minfitpts]));
    
    % Compute line incremental length and associated 'tag' index
    slen = [0,cumsum(sqrt(sum(diff([xdata,ydata],1,1).^2 ,2)),1)'];
    isl = 1:length(slen);
    % Points to mark    
    pts = 0:stepdist:slen(end);

    % Locate this points within the line
    [auxp,idx] = sort([slen,pts]);
    % Find marker points within sorted merged line (given that their 'tag' index
    % is greater than isl(end)
    ipts = find(idx>isl(end));
    auxisl = find(idx<=isl(end));

    % Define set of point triplets to use for computation of coordinates:
    % C is the marker point, i.e. the point whose coordinates we want to
    % compute; A is the point in the line that immediately precedes C, and
    % B the one that immediately follows C; D is the distance on the line
    % between A and C.
    % First allocate space for indexes and compute them
    a = zeros(length(ipts),1);
    b = zeros(length(ipts),1);
    for j = 1:isl(end)-1
        auxind = auxisl(j)+1:auxisl(j+1)-1;
%         auxind = auxisl(j):auxisl(j+1)-1;
        a(auxind-j) = idx(auxisl(j));
        b(auxind-j) = idx(auxisl(j+1));
    end
    
    % When you have two points the above does not work (might be a potential source of error, check why)
    if a==0
        a = 1;
        b = 2;
    end
    
    A = [xdata(a),ydata(a)];           
    B = [xdata(b),ydata(b)];
    D = [pts-slen(a)]';

    % Compute coordinates of points to be marked
    C = A + repmat(D./sqrt(sum((B-A).^2,2)),1,2).*(B-A);
    x = C(:,1);
    y = C(:,2);
    
    % Rescale back
    [xfin,yfin] = RescaleData(x,y,xrange,yrange,-1);
    
    % Reconvert to original axis scale
    [xfin,yfin] = loglinConversion(xfin,yfin,IsXLog,IsYLog,0);
    
    % Replace data
    set(obj(i),'XData',xfin,'YData',yfin);
end

%--------------------------------------------------------------------------
% Auxiliary functions
%--------------------------------------------------------------------------
function [xc,yc] = loglinConversion(x0,y0,IsXLog,IsYLog,mode)
%
% Convert data from/to log scale to/from linear scale given axis scales
% Mode:
% - 1|log2lin: logarithmic to linear (provide exponents)
% - 0|lin2log: linear to logarithmic
%
% Maurizio De Pitta', Tel Aviv, January 30th, 2012

switch mode
    case {0,'lin2log'}
        if IsXLog
            xc = 10.^x0;
        else
            xc = x0;
        end
        if IsYLog
            yc = 10.^y0;
        else
            yc = y0;
        end
    case {1, 'log2lin'}
        if IsXLog
            aux = x0<=0;
            if sum(aux)~=0
                warning('Y Negative data dropped!');
            end
            xc = log10(x0);
        else
            xc = x0;
        end
        if IsYLog
            aux = y0<=0;
            if sum(aux)~=0
                warning('Y Negative data dropped!');
            end
            yc = log10(y0);
        else
            yc = y0;
        end
end

% Control check
if isempty(xc)
    error('Conversion gave No X data to plot. Try different scale perhaps?');
end
if isempty(yc)
    error('Conversion gave No Y data to plot. Try different scale perhaps?');
end


function [xrange,yrange] = GetDataRanges(fh)
%
% Retrieve data min-max ranges in the plotted figure (usually needed for
% correct normalization
%
% Maurizio De Pitta', Tel Aviv, january 29th, 2012.

objh = findobj(fh,'Type','line');
xall = [];
yall = [];
for i=1:length(objh);
    xall = [xall,get(objh(i),'XData')];
    yall = [yall,get(objh(i),'YData')];
end
xrange = [min(xall),max(xall)];
yrange = [min(yall),max(yall)];

%--------------------------------------------------------------------------
function [xfin,yfin] = RescaleData(x0,y0,xrange,yrange,mode)
%
% Rescale data points x0,y0 to xfin and yfin, this latter in the exact data
% range xrange, yrange
% Mode:
% - 1|normalize: normalize x0,y0 w/ respect to their ranges xrange, yrange
% - -1|restore:  restore x0,y0 to ranges xrange, yrange giving xfin and
% yfin
%
% Maurizio De Pitta', Tel Aviv, January 28th, 2012.

switch mode
    case {1,'normalize'}
        xfin = (x0-xrange(1))./diff(xrange);
        yfin = (y0-yrange(1))./diff(yrange);
    case {-1,'restore'}
        xfin = x0*diff(xrange)+xrange(1);
        yfin = y0*diff(yrange)+yrange(1);
end
