function data = polishbd(filename,varargin)
%
% data = polishbd(filename,varargin)
%
% Read allinfo.dat-like filenames from <[filename_opts.first]> to
% <[filename_opts.last]> plotting time by time the different branches and
% allowing user to select which branches to keep and which ones to drop
% Merge all data into a single 'data' table and save it accordingly as
% 'filename'.
%
% This utility is very useful when you want to have control on nasty
% XPPAUT/AUTO data, and polish your final plot.
%
% Simple case: one data file akin to 'allinfo_1.dat'
%
% In the continuation, XPPAUTO often fails requiring user to refine
% parameters by trial-and-error before actually being able to continue
% branches successfully. While we only need data of those successfully
% continued branches, AUTO saves to file all data, that is, also those
% branches that are started but for some error are not continued and/or are
% interrupted by the user. Using polishbd on this data file allows to drop
% these broken branches providing only the informative data.
%
%
% General case: several data sets from the same bifurcation diagram,
% numbered in files akin to '<problem_name>_<val dataset>.dat', e.g.
% 'allinfo_1.dat','allinfo_2.dat', etc.
%
% When the problem is complex, memory may become an issue. Often, after
% many branches computed, XPPAUT crashes making you loose all the work. It
% is however possible to build the bifurcation diagram in steps, saving
% time by time different branches to different files. This is possible by the
% "grab" option in XPPAUT whereby we can start from any point along
% previously computed branches, to continue new branches, while discharging
%  (and thus freeing memory) previously computed branches. If we save each
%  branch set to a different file, then POLISHBD will read all branch sets
%  (or selected ones, see options) and provide the whole data set to
%  ultimately plot the complete bifurcation diagram.
%
% Input arguments:
% - filename : String with file name. NO extension has to be specified. It
%              is assumed that the file has '.dat' extension, and MUST BE in
%              the format produced by XPPAUTO / AUTO by File > All Info.
% - varargin :   Use 'option',<val> for optional input arguments (see ProduceCorrectVarargin).
%   Accepted 'option' strings are:
%   + par1   : {4} | 3   Main bifurcation par column in DAT file (>=v8.0 is 4; <v8.0 is 3)
%   + per    : {6} | 5   Period Column in DAT file (>=v7.0 is 6; <v7.0 is 5) 
%   + plot   : 0 | {1}   Interactive plot of what branches are selected.
%   + var    : {1} | Integer  What variable to show in the branch selection
%              process (variables are saved from the 6th column on in the
%              DAT files and they follow the order whereby they were
%              specified in the original ODE file.
%   + first  : {1} | Integer First branch data set to consider.
%   + last   : {1} | Integer Last branch data set to consider. It MUST be last >= first.  
%   + drop   : {[]} | Cell Array  Specify branch data sets not to consider a priori
%              in the final bifurcation data. There must be as many sets as
%              the number of file read (i.e. no branch to drop, set = []).
%              Currently works only in plot mode = 0.
%   + mergefirst : 0 | {1} (Orbits only) Join first point of the set with 
%                  last of previous set. Set it to 0 to turn it off.
%   + mergelast  : 0 | {1} (Orbits only) Join last point of the set with 
%                  first of following set. Set it to 0 to turn it off.
%   + save   : 0 | {1}  Save data to file <filename>.dat.
%
% Returns:
% - data  : MATLAB array with polished and saved data.
%
% NOTE: Known bug: the saved data file cannot be further polished by POLISHBD. 
% Additional polishing needs to be done again starting from the original 
% AUTO-produced dat files. As a rule of thumbs: DO NOT delete the original
% DAT files.
%
% see also READBD, PRODUCECORRECTVARARGIN.
%
% v2.0
% Solved issue of extra column in DAT file. Allows user to specify in case 
% Period column in DAT file for backcompatibility with XPPAUT <=v7.0 and 
% main bifurcation parameter, with backcompatibility with XPPAUT <8.0.
% Maurizio De Pitta', The University of Chicago, Chicago, October 12th, 2016. 
%
% v1.1
% Revised and added help.
% Maurizio De Pitta', The University of Chicago, Chicago, April 28th, 2016.
%
% v1.0
% Maurizio De Pitta', Tel Aviv University, Tel Aviv, Israel, January 29th 2012.
% 
% https://sites.google.com/site/mauriziodepitta/home
% maurizio.depitta@gmail.com

close all;

%--------------------------------------------------------------------------
% Defaults
%--------------------------------------------------------------------------
% Main bifurcation parameter column (XPPAUT version=>8.0 it is 4; XPPAUT version<7.0 it used to be 3
opts.par1 = 4;

% Period column (XPPAUT version=>7.0 it is 6; XPPAUT version<7.0 it used to be 5
opts.per = 6;

% Plotting / Manipulation
opts.first = 1;
opts.last = 1;
opts.plot = 1;
opts.drop = [];

% Read options
opts.mergefirst = 1;
opts.mergelast = 1;

% Variable to plot
opts.var = 1;

% Saving options
opts.save = 1;

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
% Read Files
%--------------------------------------------------------------------------
% Generate file names
NamesOfFiles = cell(1,opts.last-opts.first+1);
for i = 1:length(NamesOfFiles)
    NamesOfFiles{i} = [filename,'_',num2str(opts.first+i-1),'.dat'];
end

% Prepare opts.drop
% If specified must be {[b00,...,b0k],...,[b10,...,b1n]} with as many
% 'drops' as the branches (known a priori) otherwise set them to empty
if iscell(opts.drop)
    if length(opts.drop)==length(NamesOfFiles)
        error(['Option ''drop'' must be cell of ',num2str(length(NamesOfFiles)),' elements']);
    end
else
    opts.drop = cell(1,length(NamesOfFiles));
    for i = 1:length(NamesOfFiles)
        opts.drop{i} = [];
    end
end
    

% Read files and discard those that you select
if opts.plot==1
    h = figure;
    set(h,'WindowStyle','docked');
end
for i = 1:length(NamesOfFiles)
    fprintf(['\nReading ',NamesOfFiles{i}]);
    databd = readbd(NamesOfFiles{i},'mergefirst',opts.mergefirst,'mergelast',opts.mergelast,'par1',opts.par1,'per',opts.per);
    % Create a sample final data w/ correct column numbers (needed for
    % further concatenation)
    if i == 1
        data = zeros(0,size(databd.pts{i},2));
    end
    if opts.plot==1
        % First branch to drop
        k = 1; 
        for j = 1:size(databd.pts,2)
            % Plotting
            newh = BasicPlotting(h,databd.pts{j},databd.type{j},opts.par1,opts.per+opts.var,databd.eqs);
            if size(databd.pts{j},1)==1
                fprintf('\nIsolated point');
            end
            
            % Simple GUI to handle exceptions
            keep = 'xl';
            while or(strcmp(keep,'xl'),strcmp(keep,'yl'))
                keep = input('\nActions:: y(return)/n: Keep this set; xl: XLim; yl: YLim; e: Edit: ','s');
                switch keep
                    case 'n'
                        opts.drop{i}(k) = j;
                        k = k+1;
                        % Drop last plotted data from plot
                        delete(newh)
                    case 'xl'
                        XLim = input('\nEnter XLim: ');
                        set(gca,'XLim',XLim);
                    case 'yl'
                        YLim = input('\nEnter YLim: ');
                        set(gca,'YLim',YLim);
                    case 'e'
                        axlim = axis;
                        databd.pts{j} = EditData(databd.pts{j},opts.par1,opts.var,axlim);
                        delete(newh);
                        BasicPlotting(h,databd.pts{j},databd.type{j},opts.par1,opts.per+opts.var,databd.eqs);
                end
            end
            
            % Recolorize and prepare for next data set
            Recolorize(h);
        end
    end
    % Keep only data specified by user
    if ~isempty(opts.drop)
        databd.pts(opts.drop{i}) = [];
    end;
    data = vertcat(data,vertcat(databd.pts{:}));
end

if opts.save==1
    save([filename,'.dat'],'data','-ascii','-double','-tabs');
end

%--------------------------------------------------------------------------
% Auxiliary functions
%--------------------------------------------------------------------------
function pts = EditData(pts,par1,icol,axlim)
%
% Basic editing of data
%
% v1.1
% Added par1 as input. 
% Maurizio De Pitta', Chicago, Nov 2, 2016.
%  
% v1.0
% Maurizio De Pitta', Tel Aviv, January 29th, 2012.

% Prompt user for editing action 
sel = input(['\n1: Resize data range',...
             '\n2: Flip data',... 
             '\nAction: ']);
switch sel
    case 1
        % Resize range of data
        xrange = input('\nXRange=');
        pts(pts(:,par1)<xrange(1)|pts(:,par1)>xrange(2),:) = [];
    case 2
        pts = flipud(pts);
    otherwise
        % Do nothing
        pts = pts;
end

%--------------------------------------------------------------------------
function Recolorize(h)
%
% Plotting routine that consider a given figure, handled by handler 'h' and
% recolorize data in grey (as if the data have been already treated)
%
% Maurizio De Pitta', Tel Aviv, January 29th, 2012.

Color = 0.75.*[1,1,1];

obj = findobj(h,'Type','line');
for i = 1:length(obj)
    set(obj(i),'Color',Color','MarkerFaceColor',Color','MarkerEdgeColor',Color);
end

%--------------------------------------------------------------------------
function newh = BasicPlotting(h,pset,type,par1,icol,neqs)
%
% Basic plotting of branches (default view as in XPPAUT)
%
% v1.1
% Added par1 as input. 
% Maurizio De Pitta', Chicago, Nov 2, 2016.
%  
% v1.0
% Maurizio De Pitta', Tel Aviv, January 29th, 2012.

LineWidth = 2.5;
MarkerSize = 8;

% First retrieve old handles of data
oldh = findobj(h,'Type','line');

hold on;
switch type
    case 'se'
        plot(pset(:,par1),pset(:,icol),...
            'LineStyle','-',...
            'LineWidth',LineWidth,...
            'Color','r',...
            'Marker','none');
    case 'ue'
        plot(pset(:,par1),pset(:,icol),...
            'LineStyle','-',...
            'LineWidth',LineWidth,...
            'Color','k',...
            'Marker','none');
    case 'so'
         % Maxima
         plot(pset(:,par1),pset(:,icol),pset(:,par1),pset(:,neqs+icol),...
            'LineStyle','none',...
            'LineWidth',LineWidth,...
            'Color','g',...
            'Marker','.',...
            'MarkerSize',MarkerSize,...
            'MarkerEdgeColor','g',...
            'MarkerFaceColor','g');
    case 'uo'
            plot(pset(:,par1),pset(:,icol),pset(:,par1),pset(:,neqs+icol),...
                'LineStyle','none',...
                'LineWidth',LineWidth,...
                'Color','b',...
                'Marker','o',...
                'MarkerSize',MarkerSize,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','none');
end
hold off;  

% Retrieve all handles
newh = findobj(h,'Type','line');

% Get new handles
newh = setdiff(newh,oldh);