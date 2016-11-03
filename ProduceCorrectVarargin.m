function TrueVarargin = ProduceCorrectVarargin(VararginCell)
%
%TrueVarargin = ProduceCorrectVarargin(VararginCell)
%
%Utility routine needed to correct varargin when passed in functions called
%by other functions 
%NOTE: if varargin is already in the right form it turns it back in the output
%(the condition is checked within the function).
%
%Input arguments:
% - VararginCell: must be a vector of cells (in the worst case):
%   VararginCell = [1x1 cell] such that:
%   VararginCell{1} = varargin = {'String1',Val1,...,'StringN',ValN}
%
%Output arguments:
% - TrueVarargin: is varargin = {'String1',Val1,...,'StringN',ValN}
% suitable for use in all functions of the GUI
%
%Ex.: 
% 
% varargin = ProduceCorrectVarargin(varargin)
%
% Directly produces varargin in the correct form
%
% NOTE: if TrueVarargin is not provided in the stated form then it is
% likely that varargin is not provided in the required form.
%
% v4.0
% Solved back-compatibility issues with code in MATLAB previous to R2013a.
% Maurizio De Pitta', The University of Chicago, Chicago, April 28th, 2016.
%
% v3.0 
% Update handling of TrueVarargin: direct indexing / this way can
% handle any type of oject in a cell, such as function handles or '[]'
% empty vectors for example
% Maurizio De Pitta', Tel Aviv University, Tel Aviv, Israel, January 12th, 2012.
%
% v2.0 
% Remove Multiple identical varargins, leaving only the last-defined one
% Maurizio De Pitta', Tel Aviv University, Tel Aviv, Israel, September 26th, 2010.
%
% v1.0
% Maurizio De Pitta', Tel Aviv University, Tel Aviv, Israel, June 23rd, 2007.

% OLD
TrueVarargin = {};
ind = 0;
if isempty(VararginCell)==0
    if iscell(VararginCell)&&iscell(VararginCell{1})
%         TrueVarargin = [];
        for i = 1:length(VararginCell)
            for j = 1:length(VararginCell{i})
                % OLD IMPLEMENTATION
                TrueVarargin{ind+j} = VararginCell{i}{j};
                % CORRECTED
%                 if ~isempty(VararginCell{i}{j})
%                     TrueVarargin{ind+j} = VararginCell{i}{j};
%                 end
                % OBSOLETE
%                 if isa(VararginCell{i}{j},'function_handle')
%                     TrueVarargin{length(TrueVarargin)+1} = VararginCell{i}{j};
%                 else
%                     TrueVarargin = [TrueVarargin VararginCell{i}{j}];
%                 end
            end
            % Save last index
            % OLD
            if ~isempty(j)
                % Case when VararginCell{i} is not empty
                ind = j;
            end
        end
    else
        TrueVarargin = VararginCell;
    end
else
    TrueVarargin = {};
end

% Remove multiple fields leaving the last-defined one
% NOTE: If TrueVarargin is not 'string1,<val1>,'string1',<val2>... then you
% get and error here (modify the root routine of varargin: you are not
% likely passing varargin correctly
if ~isempty(TrueVarargin)
    [Tags,index] = unique(TrueVarargin(1:2:end),'last');
    index = sort([2*index-1,2*index]);
    TrueVarargin = TrueVarargin(index);
    % Somewhen MATLAB started return unique indexes as row vectors. The
    % following line is added for back compatibility with code previous to
    % R2013a
    TrueVarargin = reshape(TrueVarargin',[],numel(TrueVarargin));
end
    
