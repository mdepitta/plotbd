function h = plotbd(databd,varargin)
%
% h = plotbd(databd,varargin)
% h = plotbd(databd,'var',1)
% h = plotbd(databd,'var',1,'interp',1,...)
%
% Plots bifurcation data from XPPAUT/AUTO.
% 
% Input arguments:
% - databd   :  Structure of bifurcation data as produced by READBD.
% - varargin :  Use 'option',<val> for optional input arguments (see ProduceCorrectVarargin).
%   Accepted 'option' strings are:
% Basic options:
%   + var    : {1} | Integer  What variable to plot 
%              (variables are saved from the 6th column on in the
%              DAT files and they follow the order whereby they were
%              specified in the original ODE file.
%   + pthr   : {1} | Integer  Max number of isolated fixed points to filter out.
%   + othr   : {1} | Integer  Max number of isolated orbit points to filter out.
%   + period : {0} | 1        Oscillation period plotting. Default: does not plot period of oscillations.
%   + freq   :  0  | {1}      If period==1, show frequency of oscillatios instead of period.
%   + noe    : {0} | 1        Exclude equilibria (show only oscillation branches). Default: always plot equilibria.
%
% Plotting options:
%   + fh     : {[]}| Figure handle.  Default: create a new figure.
%                             If figure handle is given instead and bd needs 
%                             be superimposed, than fh must be a vector of 
%                             exactly the same number of elements of h
%                             (i.e. the return argument).
%   + interp : {0} | 1        CSPLINE interpolation of data to plot things nicely. Default: does not interpolate.
%                             Assume that you master AUTO or AUTO has been nice with you. The changes for the last 
%                             are often equivalent to me believing in God...
%   + stype  : {'pchip'} | 'spline'   Spline interpolator type.
%   + nice = : {0} | 1        Post-process plotted bifurcation diagram doing some magic by NICEPLOT.
%   + special: {[]}| Array    (only with nice==1) Specify which branches from the
%              plotted bifurcation diagram are post-processed. The array must
%              be in the form of pairs like [graphobj1, step1,... graphobjn,stepn]
%              Default: branches are treated in the same fashion.
%   + dense  : {2} | Integer   (only with nice==1) Considers 'dense' times 
%              the number of points originally computed by AUTO in each branch to
%              perform cspline interpolation in the (x,y) domains.
%   + step   : 2e-2| Double    (only with nice==1) Space between points along 
%              an interpolated branch.
%   + scol   : {1} | Integer   Sampling point number to plot stable point data.
%   + ucol   : {1} | Integer   Sampling point number to plot unstable point data.
%   + normalize : {0} | 1      Normalize data w/ respect to their maxima.
% 
% Axis options:
%   + xscale, yscale: {'linear'} or {'lin'} | 'log'  Use linear or
%                     logarithmic scales for x,y axes.
%   + lims   : {'auto'} | 'manual' Automatic adjustment of axis limits.
%   + xlim, ylim : {[]} | [inf,sup] Axis limits (only for lims=='manual').
% 
% Graphic options:
%   + box    : 0 | {1}  Defualt: Draw Box around figure.
%   + afs    : {12} | Value      Axis Font Size for axis ticks.
%   + alw    : {1.5}| Value      Axis Line Width (in pt).
%   + MarkerSize      : {5} | Value   Marker Size.
%   + MarkerLineWidth : {1} | Value   Marker Line Width (in pt).
%   + factor : {0.1} | Value  Factor to adjust exat limits by to avoid data 
%                             touching the axes.
%   + FontName   : {'Arial'} | Font Name   Axis Label Font.
%   + FontWeight : {'bold'}  | Style       Axis Label Font Weight.
%   + FontSize   : {14} | Value            Axis Label Font Size.
%
% Line Colors ([min max]):
%   + sec : {'k'} | Color    Stable equilibrium color.
%   + uec : {'r'} | Color    Unstable equilibrium color.
%   + soc : {'k'} | Color    Stable Orbit color.
%   + uoc : {'k'} | Color    Unstable Orbit color.
% Line Style:
%   + ses : {'-'} | line style    Stable equilibrium line style.
%   + ues : {'--'}| line style    Unstable equilibrium line style.
%   + sos : {'none'}| line style  Stable Orbit line style.
%   + uos : {'none'}| line style  Unstable Orbit line style.
% LineWidth:
%   + selw : {2.1}| line width  Stable equilibrium line width.
%   + uelw : {2.1}| line width  Unstable equilibrium line width.
%   + solw : {2.1}| line width  Stable Orbit line width.
%   + uolw : {2.1}| line width  Unstable Orbit line width.
%
% Markers:
%   + sem : {'none'} | Marker    Stable equilibrium marker.
%   + uem : {'none'} | Marker    Unstable equilibrium marker.
%   + som : {'o'} | Marker       Stable Orbit marker.
%   + uom : {'o'} | Marker       Unstable Orbit marker.
% MarkerFaceColor:
%   + semfc : {'w'} | Color      Stable equilibrium marker facecolor. 
%   + uemfc : {'w'} | Color      Unstable equilibrium marker facecolor.
%   + somfc : {'k'} | Color      Stable Orbit marker facecolor.
%   + uomfc : {'none'} | Color   Unstable Orbit marker facecolor.
% MarkerEdgeColor:
%   + semec : {'none'} | Color   Stable equilibrium marker edgecolor. 
%   + uemec : {'none'} | Color   Unstable equilibrium marker edgecolor.
%   + somec : {'k'} | Color      Stable Orbit marker edgecolor.
%   + uomec : {'k'} | Color      Unstable Orbit marker edgecolor.
%
%
% Returns:
% - h        :  figure handle where bifurcation diagram is plotted. 
%
% see also READBD, POLISHBD, NICEPLOT, CSPLINE, PRODUCECORRECTVARARGIN.
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
% Main bifurcation parameter column (XPPAUT version=>8.0 it is 4; XPPAUT version<7.0 it used to be 3
opts.par1 = 4;

% Period column (XPPAUT version=>7.0 it is 6; XPPAUT version<7.0 it used to be 5
opts.per = 6;

% Variables that you want to plot (can also be a vector)
opts.var = 1;

% Filtering threshold for isolated points
opts.pthr = 1;
opts.othr = 1;

% Period plotting
opts.period = 0;    % Default: does not plot period. 0/1 {do not plot period} | plot period;
opts.freq = 1;      % If opts.period==1; plots frequency instead of period      

% Exclusion of equilibria
opts.noe = 0;       % Default: plots always equilibria

% Plotting options
opts.interp = 0;      % spline interpolation of data (nice it works!)
opts.stype = 'pchip'; % spline interpolation type
opts.nice = 1;        % refine visualized data by nice-plotting ;)
opts.special = [];    % (niceplot option): treat all branches in the same fashion. When specified is an array of pairs like [graphobj1, step1,... graphobjn,stepn]
opts.dense = 2;       % (niceplot option): density factor, counts points times this factor for spline fitting
opts.step = 2e-2;     % (niceplot option): distance step to plot markers
opts.scol = 1;        % point step to pick up data
opts.ucol = 1;        % point step to pick up data
opts.normalize = 0;   % normalize data w/ respect to their maxima (across all compartments)

% Figure handles
opts.fh = [];         % Default: create new figures. If instead plotting must be superimposed, than handles must be a vector of exactly the same number of elements of h

% Axis and 
opts.xscale = 'lin';  % linear X scale for spline interpolation;
opts.yscale = 'lin';  % linear Y scale for spline interpolation;
opts.lims = 'auto';   % automatic adjustment of axis (alternatively: 'manual')
opts.xlim = [];
opts.ylim = [];

% Graphic settings
opts.box = 1;       % Box around figure on (otherwise: 0: off)
opts.afs = 12;      % Axis Font Size
opts.alw = 1.5;     % Axis Line Width

opts.MarkerSize = 5;
opts.MarkerLineWidth = 1;
opts.factor = 1;

opts.FontName = 'Arial';
opts.FontWeight = 'bold';
opts.FontSize = 14;

% Line Colors ([min max])
opts.sec = 'k';
opts.uec = 'r';
opts.soc = 'k';
opts.uoc = 'k';

% Markers
opts.sem = 'none';
opts.uem = 'none';
opts.som = 'o';
opts.uom = 'o';

% MarkerFaceColor
opts.semfc = 'w';
opts.uemfc = 'w';
opts.somfc = 'k';
opts.uomfc = 'none';

% MarkerEdgeColor
opts.semec = 'none';
opts.uemec = 'none';
opts.somec = 'k';
opts.uomec = 'k';

% Line Style
opts.ses = '-';
opts.ues = '--';
opts.sos = 'none';
opts.uos = 'none';

% LineWidth
% opts.LineWidth = 2.1;
% % Could also assign case by case LineWidth, but perhaps it is not worth it
opts.selw = 2.1;
opts.uelw = 2.1;
opts.solw = 2.1;
opts.uolw = 2.1;

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

% Adapt-to-MATLAB interface some parameters
switch opts.box
    case 0
        opts.box = 'off';
    case 1
        opts.box = 'on';
end

if strcmp(opts.xscale,'lin')
    opts.xscale = 'linear';
end
if strcmp(opts.yscale,'lin')
    opts.yscale = 'linear';
end

%--------------------------------------------------------------------------
% Data format specifics (it is useful to save them also within the opts struct)
%--------------------------------------------------------------------------
opts.par1 = databd.par1;
opts.per = databd.per;

%--------------------------------------------------------------------------
% Pre-processing of data
%--------------------------------------------------------------------------
iset = 1;   % index of data2plot.pts
for i = 1:length(databd.pts)
    % Temporary aux variable
    dataset = databd.pts{i};
    % If dataset is only less than 'thr' number of points then drops it
    switch databd.type{i}
        case {'se','ue'}
            thr = opts.pthr;
        case {'so','uo'}
            thr = opts.othr;
    end
    if (size(dataset,1)<=thr)||((opts.noe==1)&&(strcmp(databd.type{i},'se')||strcmp(databd.type{i},'ue')))
        continue;
    end
    if opts.interp==1
        %  Do something
        return;
    end
    % Tag to relate data2plot.pts to the original databd.pts set
    tagi(iset) = i;
    % Data to plot and type
    data2plot.pts{iset} = dataset;
    data2plot.type{iset} = databd.type{i};
    % update index
    iset = iset+1;
end

%--------------------------------------------------------------------------
% Re-order data2plot sets
%--------------------------------------------------------------------------
% In typical plotting settings, you want stable orbits to show off better
% than other orbits. And in very complicated bifurcation diagrams, this is
% achieved plotting them last
for i = 1:length(data2plot.pts)
    auxi(i) = PointTypeInverse(data2plot.type{i});
end
% Makes the stable orbit point the last in the index
i3 = auxi==3;
i4 = auxi==4;
auxi(i3) = 4;
auxi(i4) = 3;
[auxi,idx] = sort(auxi);
% Permute accordingly data2plot and tagi
tagi = tagi(idx);
data2plot.pts = data2plot.pts(idx);
data2plot.type = data2plot.type(idx);

%--------------------------------------------------------------------------
% Normalize data if option is selected
%--------------------------------------------------------------------------
if opts.normalize
    % NOTE::DOES NOT NORMALIZE EIGENVALUES!!!
    % Take all data to plot and build auxiliary matrix of all values to
    % compute min/max from
    aux = vertcat(data2plot.pts{:});
    % Considers only variable columns
    aux = aux(:,databd.per+1:databd.per+2*databd.eqs);
    % Compute maxima and minima of variables
    vmin = min(aux);
    vmax = max(aux);
    % Rearrange values internally to treat orbits correctly in the
    % normalization
    vmin(1:databd.eqs) = vmin(databd.eqs+1:end);
    vmax(databd.eqs+1:end) = vmax(1:databd.eqs);
    % Effective normalization
    for i = 1:length(data2plot.pts)
        % Effective data normalization
%         data2plot.pts{i}(:,6:5+databd.eqs) = (data2plot.pts{i}(:,6:5+databd.eqs)-repmat(vmin,size(data2plot.pts{i},1),1))./repmat(vmax-vmin,size(data2plot.pts{i},1),1));
        % Normalize data only by their maxima
        data2plot.pts{i}(:,databd.per+1:databd.per+2*databd.eqs) = data2plot.pts{i}(:,databd.per+1:databd.per+2*databd.eqs)./repmat(vmax,size(data2plot.pts{i},1),1);
    end
end

%--------------------------------------------------------------------------
% Effective plotting
%--------------------------------------------------------------------------
% Number of graphic objects to be plotted
[nobj,objval] = ComputeNumObjects(data2plot.type,opts.period);
% Plot variables or period
switch opts.period
    case 0
        % Standard bifurcation diagram
        for j = 1:length(opts.var)
            if isempty(opts.fh)
                h(j) = figure;
                %             set(gca,'units','centimeters');
            else
                h(j) = figure(opts.fh(j));
            end
            for i = 1:length(data2plot.pts)
                PlotBDPoints(data2plot.pts{i},data2plot.type{i},databd.par1,j,opts,databd.eqs,tagi(i),nobj-sum(objval(1:i))+1);
            end
            % If period plotting is selectedd needs to plot only one figure
        end
    case 1
        % Period plotting
        if isempty(opts.fh)
            h(1) = figure;
            %         set(gca,'units','centimeters');
        else
            h(1) = figure(opts.fh(1));
        end
        for i = 1:length(data2plot.pts)
            switch data2plot.type{i}
                case {'se','ue'}
                    % do nothing
                case {'so','uo'}
                    PlotBDPoints(data2plot.pts{i},data2plot.type{i},databd.par1,1,opts,databd.eqs,tagi(i),nobj-sum(objval(1:i))+1);
            end
        end
end

%--------------------------------------------------------------------------
% Resize plotting and some basic restyling (box, alw, afs, xscale, yscale)
%--------------------------------------------------------------------------
% Overrun opts.lim option in case it is plotting on previously existitng
% figure (this is assuming that in the existing figure a bifurcation
% diagram with already nice axis settings exists and you do not want to
% loose these latter)
if ~isempty(opts.fh)
    opts.xlim = get(gca,'XLim');
    opts.ylim = get(gca,'YLim');
    opts.lim = 'manual';
end

% Define axis limits
switch opts.lims
    case 'manual'
        % Manual adjustment of axis limits
        XLim = opts.xlim;
        YLim = opts.ylim;
    case 'auto'
        % Automatically adjust axis values accordingly to min and max of
        % data
        switch opts.period
            case 0
                auxdata = vertcat(data2plot.pts{:});
                XLim = [min(auxdata(:,databd.par1)),max(auxdata(:,databd.par1))];
                YLim = [min(min(auxdata(:,databd.eqs+databd.per+opts.var))),max(max(auxdata(:,databd.per+opts.var)))];
            case 1
                % Select only data that contains orbit points
                auxdata = data2plot.pts((auxi==3|auxi==4));  % NOTE: This might not work if more than 2 branches for orbits are present
                auxdata = vertcat(auxdata{:});
                XLim = [min(auxdata(:,databd.par1)),max(auxdata(:,databd.par1))];
                YLim = [min(auxdata(:,databd.per)),max(auxdata(:,databd.per))];
        end
        % adjust exact limits by a factor .1 to make data not touch the axis
        % NOTE: I think that .1 factor could fit most cases, worse case add
        % option or use 'manual'
        XLim = [XLim(:,1)-.1*diff(XLim),XLim(:,2)+.1*diff(XLim)];
        YLim = [YLim(:,1)-.1*diff(YLim),YLim(:,2)+.1*diff(YLim)];
end

% Effective resizing
switch opts.period
    case 0
        % Count handles
        for i = 1:length(h)
            % Assumes that figures have only one Axis child
            ax = get(h(i),'Children');
            set(ax,'XLim',XLim,'YLim',YLim,'XScale',opts.xscale,'YScale',opts.yscale,...
                   'Box',opts.box,'LineWidth',opts.alw,'FontSize',opts.afs);
        end
    case 1
        ax = get(h(1),'Children');
        set(ax,'XLim',XLim,'YLim',YLim,'XScale',opts.xscale,'YScale',opts.yscale,...
               'Box',opts.box,'LineWidth',opts.alw,'FontSize',opts.afs);
end

%--------------------------------------------------------------------------
% Refine plotting (opts.interp==1)
%--------------------------------------------------------------------------
if opts.nice==1
%     switch opts.period
%         case 0
%             % Count handles
            for i = 1:length(h)
                % Assumes that figures have only one Axis child
                niceplot(h(i),'MarkerSize',opts.MarkerSize,'factor',opts.factor,'special',opts.special,'dense',opts.dense,'step',opts.step);
            end
%         case 1
%             niceplot(h(1),'MarkerSize',opts.MarkerSize,'factor',opts.factor,'special',opts.special,'dense',opts.dense,'step',opts.step);
%     end
end

%--------------------------------------------------------------------------
% Auxiliart functions
%--------------------------------------------------------------------------
function [nobj,objval] = ComputeNumObjects(list_obj_type,period_flag)
%
% Utility that computes a priori the number of graphic objects to be
% included within the final figure, based on list_obj_type (that is 
% data2plot.type in our case).
%
% Maurizio De Pitta', Tel Aviv, January 30th, 2012.

% First classify (only orbits have a contribution in terms of graphic
% objects that depends on period_flag)
orb_idx = strcmp([list_obj_type(:)],'so')|strcmp([list_obj_type(:)],'uo');

% Assign contribution to stack index of each branch depending on period
% plotting option
objval = ones(1,length(list_obj_type));
if period_flag==0
    % Each orbit branch contributes with two graphic objects: max and
    % min in a regular bifurcation diagram otherwise just with one onject
    % in period plot
    objval(orb_idx) = 2;
end

% Compute number of graphic objects to be plotted
nobj = sum(objval);

%--------------------------------------------------------------------------
function auxi = PointTypeInverse(ptype)
%
% Assign index to type as specified by stype in databd
% It performs the reverse operation of PointType in readbd
%
% Maurizio De Pitta', Tel Aviv, January 28th, 2012

switch ptype
    case 'se'
        auxi = 1;
    case 'ue'
        auxi = 2;
    case 'so'
        auxi = 3;
    case 'uo'
        auxi = 4;
end
        
%--------------------------------------------------------------------------
function PlotBDPoints(pset,type,par1,ivar,opts,neqs,idx,idx_stack)
% Plot bifurcation data according to type
% 
% Inputs:
% - pset: data set to plot
% - type: type of data set
% - par1: index of column of main bifurcation parameter
% - ivar: index of variable to plot
% - opts: option structure as generated by plotbd
% - eqs:  equations of the model
% - idx:  branch index (namely what branch of data2plot.pts you are plotting)
% % 'idx' is in the case you need to identify set of data that need manual
% % editing: it will be saved in the 'Tag' string of each plotted data set
%
% v1.1
% Added par1 as input. 
% Maurizio De Pitta', Chicago, Nov 2, 2016.
%  
% v1.0
% Maurizio De Pitta', Tel Aviv, January 28th, 2012.

% Define plotting parameters according to type
switch type
    case 'se'
        LineStyle = opts.ses;
        LineWidth = opts.selw;
        Color = opts.sec;
        Marker = opts.sem;
        MarkerSize = opts.MarkerSize;
        MarkerEdgeColor = opts.semec;
        MarkerFaceColor = opts.semfc;
    case 'ue'
        LineStyle = opts.ues;
        LineWidth = opts.uelw;
        Color = opts.uec;
        Marker = opts.uem;
        MarkerSize = opts.MarkerSize;
        MarkerEdgeColor = opts.uemec;
        MarkerFaceColor = opts.uemfc;
    case 'so'
        LineStyle = opts.sos;
        LineWidth = opts.solw;
        Color = opts.soc;
        Marker = opts.som;
        MarkerSize = opts.MarkerSize;
        MarkerEdgeColor = opts.somec;
        MarkerFaceColor = opts.somfc;
    case 'uo'
        LineStyle = opts.uos;
        LineWidth = opts.uolw;
        Color = opts.uoc;
        Marker = opts.uom;
        MarkerSize = opts.MarkerSize;
        MarkerEdgeColor = opts.uomec;
        MarkerFaceColor = opts.uomfc;
end

switch opts.period
    case 0
        % Plot variables
        icol = opts.per+opts.var(ivar);
    case 1
        % Plots period
        icol = opts.per;
        % if opts.freq on then plot frequency data
        if opts.freq==1
            pset(:,icol) = 1./pset(:,icol);
%             % Some error checking in case some frequency data goes to zero
%             % (usually impossible)
%             aux = find(pset(:,icol)==0);
%             % if you have data that go to zero with frquency means that
%             you have period that goes to INF
%             if ~isempty(aux)
%                 warning('Frequency computation!!::DIV0::niceplot option might fail if set!');
%             end
        end
end

% Plotting
hold on;
switch type
    case {'se','ue'}
        plot(pset(:,par1),pset(:,icol),...
            'LineStyle',LineStyle,...
            'LineWidth',LineWidth,...
            'Color',Color,...
            'Marker',Marker,...
            'MarkerSize',MarkerSize,...
            'MarkerEdgeColor',MarkerEdgeColor,...
            'MarkerFaceColor',MarkerFaceColor,...
            'Tag',num2str([idx,idx_stack]));
    case {'so','uo'}
         hold on;
         % Maxima
         plot(pset(:,par1),pset(:,icol),...
            'LineStyle',LineStyle,...
            'LineWidth',LineWidth,...
            'Color',Color,...
            'Marker',Marker,...
            'MarkerSize',MarkerSize,...
            'MarkerEdgeColor',MarkerEdgeColor,...
            'MarkerFaceColor',MarkerFaceColor,...
            'Tag',num2str([idx,idx_stack+1]));
        % Minima (only if not period plotting)
        if ~opts.period
            plot(pset(:,par1),pset(:,neqs+icol),...
                'LineStyle',LineStyle,...
                'LineWidth',LineWidth,...
                'Color',Color,...
                'Marker',Marker,...
                'MarkerSize',MarkerSize,...
                'MarkerEdgeColor',MarkerEdgeColor,...
                'MarkerFaceColor',MarkerFaceColor,...
                'Tag',num2str([idx,idx_stack]));
        end
end
hold off;        