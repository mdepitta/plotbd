% plotbd_example.m
%
% Sample Script to build a bifurcation diagram whose data was
% "sequentially" saved to files with labels '<filename>_i.dat' where
% i=1,2,...
%
% Maurizio De Pitta', The University of Chicago, Chicago, April 28th, 2016.
% 
% https://sites.google.com/site/mauriziodepitta/home
% maurizio.depitta@gmail.com

%--------------------------------------------------------------------------
% Plot bifurcation diagram w/out polishing data
%--------------------------------------------------------------------------
% First read data as they have been produced by XPPAUT and plot them
% without polishing them.
databd = readbd('allinfo_edited_1.dat');
plotbd(databd);

%--------------------------------------------------------------------------
% Plot bifurcation diagram w/ polishing data first
%--------------------------------------------------------------------------
% First load all your bifurcation data, allows to interactively polish
% them (see command line options) and save them to '<filename>.dat' (that
% is no number!)
data = polishbd('allinfo_edited'); % Need to specify only file name (and in case which data sets)

% Now the file 'allinfo_edited.dat' is produced with all the polished data
% to build the bifurcation diagram from. This file needs to be read and its
% info saved in the databd structure:
databd = readbd('allinfo_edited.dat');

% Finally databd structure is passed to the plotting routine
plotbd(databd);