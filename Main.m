% *************************************************************************
% 
% This script is used in ASEN 3111 Exp. Lab #1: Wake Aerodynamics (Gp. #10)
% 
% This script:
%       - Reads in experimental data
%       - Averages each point capture
% 		- Calculates drag for each experiment
% 		- Calculates velocity deficit and half wake for each experiment
% 		- Creates velocity deficit plots
% 		- Plots the velocity deficit profile
% 		- Plots wake half-width vs. x-position
% 		- Sorts data
%
% Authors:      Lucas Calvert
%               Keith Covington
%               Daniel Mastick
%               Ginger Beerman
% 
% Created:      09/08/2017
% Modified:     11/03/2017
%
% *************************************************************************


% Housekeeping
clearvars -except allData % clear everything except raw data
% NOTE: This is just for convenience so you don't have to watch it load 
% 		data every time you want to test something.
% 		If you add or delete data files in the 'Data' folder, run 'clear'
% 		in the command window before you run 'Main'.
close all; clc;
addpath('Data')


%% Define filenames

% ex. filename: ExperimentalLab1_Section1_V15_x28_Airfoil_Case1
filenames = cell(120, 1);							% 120 total files to read (jeeeez)
labSections = repelem(1:3, 1, 40);					% 40 files for each lab section
velocities = repmat([15 15 25 25], 1, 30);			% 30 groups test 2 @ 15 m/s and 2 @ 30 m/s
locations = repmat(repelem( ...
	[90 180 240 13 18 23 28 33 38 43],1,4), 1,3);	% 10 x-locations for 4 tests for 3 sections
objectTypes = repmat(repelem({'Cylinder', ...
	'Cylinder', 'Cylinder', 'Airfoil', ...
	'Airfoil', 'Airfoil', 'Airfoil', ...
	'Airfoil', 'Airfoil', 'Airfoil'},1,4), 1,3);	% 10 object types for 4 tests for 3 sections
testCases = repmat([1 2 1 2], 1, 30);				% 4 tests conducted by 30 total groups

% Synthesize names for all 120 tests
for file = 1:length(filenames)
    filenames{file} = [ ...
		'ExperimentalLab1_Section'	num2str(labSections(file)) ... 
        '_V'						num2str(velocities(file)) ... 
        '_x'						num2str(locations(file)) ...
        '_'							objectTypes{file} ...
        '_Case'						num2str(testCases(file)) ...
		'.csv'];
end
%disp(filenames);
%ls Data


%% Read in data
if ~exist('allData')		% only read in data if it doesn't exist in workspace

	% If data for specific files are bad, replace filenames to read with replacement file with good data
	filenames_bad = {	'ExperimentalLab1_Section3_V25_x28_Airfoil_Case2.csv'
						'ExperimentalLab1_Section3_V15_x28_Airfoil_Case2.csv'
						'ExperimentalLab1_Section3_V25_x18_Airfoil_Case1.csv'
						'ExperimentalLab1_Section2_V15_x33_Airfoil_Case2.csv'
						'ExperimentalLab1_Section2_V25_x23_Airfoil_Case2.csv'
						'ExperimentalLab1_Section2_V15_x23_Airfoil_Case1.csv'
						'ExperimentalLab1_Section2_V15_x240_Cylinder_Case1.csv'
					};	% files we don't want to use
	filenames_good = {'','','','','','',''};	% files we'll use instead
	for bFiles = 1:length(filenames_bad)
		badFileInd = find(strcmp(filenames_bad{bFiles}, filenames));	% find indices of bad files in array of filenames
		filenames{badFileInd} = filenames_good{bFiles};	% replace name of bad file with name of good file
	end

	% Read the data
	firstFile = true;
	wbar = waitbar(0,'Reading data...');
	for file = 1:length(filenames)
		% If file exits, read it
		if ~exist(filenames{file}, 'file') == 0
			disp(['reading  ' filenames{file} ' ...']);

			try
				[A,N,eld_y,scal_press,q_inf,t_atm,M,aoa,p_atm,rho_atm,aux_inf,...
					eld_x,v_inf, v] = load_data( filenames{file} ); % read file
					
					% create allData if first iteration
					if firstFile
						%allData = NaN(length(A), nargout('load_data'), length(filenames));
						dataFileSize = size([A,N,eld_y,scal_press,q_inf,t_atm, ...
											 M,aoa,p_atm,rho_atm,aux_inf,eld_x,v_inf, v]);
						allData = NaN([dataFileSize, length(filenames)]); % preallocate
						firstFile = false;
					end
					
					% add data to allData matrix
					allData(:,:,file) = [A,N,eld_y,scal_press,q_inf,t_atm,M, ...
										aoa,p_atm,rho_atm,aux_inf,eld_x,v_inf,v];
			catch
				disp(['File:  ' filenames{file} '  could not be read...skipping file.'])
			end
		else
			disp(['File:  ' filenames{file} '  is missing...skipping file.'])
		end
		waitbar(file/length(filenames));
	end
	close(wbar);
end

% Re-synthesize names for all 120 tests
for file = 1:length(filenames)
    filenames{file} = [ ...
		'ExperimentalLab1_Section'	num2str(labSections(file)) ... 
        '_V'						num2str(velocities(file)) ... 
        '_x'						num2str(locations(file)) ...
        '_'							objectTypes{file} ...
        '_Case'						num2str(testCases(file)) ...
		'.csv'];
end
dataSize = size(allData);	% get size of matrix containing all raw data


%% Average data for each 20-point capture
dividedDataCell = mat2cell(allData, repelem(500, 20), dataSize(2), dataSize(3));
avgDataCell = cellfun(@(C) nanmean(C,1), dividedDataCell, 'UniformOutput',false);
avgData = cell2mat(avgDataCell);

% Show that data is averaged by demonstrating averaging of eld_y probe position
%{
figure
plot(allData(:,3,1), 'o')
figure
plot(avgData(:,3,1), 'o')
%}


%% Do things with averaged, unsorted data

% Extract data from matrix...3rd dimension is test #
A_avg_all			= avgData(:,1,:);
N_avg_all			= avgData(:,2,:);
eld_y_avg_all		= avgData(:,3,:);
scal_press_avg_all	= avgData(:,4:18,:);
q_inf_avg_all		= avgData(:,19,:);
t_atm_avg_all		= avgData(:,20,:);
M_avg_all			= avgData(:,21,:);
aoa_avg_all			= avgData(:,22,:);
p_atm_avg_all		= avgData(:,23,:);
rho_atm_avg_all		= avgData(:,24,:);
aux_inf_avg_all		= avgData(:,25,:);
eld_x_avg_all		= avgData(:,26,:);
v_inf_avg_all		= avgData(:,27,:);
v_avg_all			= avgData(:,28,:);


% Calculate drag
for i = 1:length(A_avg_all(1,1,:))
	dPrime(i) = calcD(v_inf_avg_all(:,:,i), 	...
					  rho_atm_avg_all(:,:,i), 	...
					  v_avg_all(:,:,i), 		...
					  eld_y_avg_all(:,:,i));
end

[a15avgCd,a25avgCd,c15avgCd,c25avgCd] = compareDrag(dPrime);


%% Plot Velocity Deficit

% Preallocate
vel_def 	= NaN(length(A_avg_all(1,1,:)), 20);
wake_half 	= NaN(length(A_avg_all(1,1,:)), 1);
vel_def_max	= NaN(length(A_avg_all(1,1,:)), 1);

% Loop through all experiments
for i = 1:length(A_avg_all(1,1,:))
	try
		[ vel_def(i,:), wake_half(i) ,vel_def_max(i) ] ...
			= vel_deficit(	...
							eld_x_avg_all(:,:,i), 	...
							eld_y_avg_all(:,:,i), 	...
							v_inf_avg_all(:,:,i), 	...
							v_avg_all(:,:,i), 		...
							'PlotsOff');
							%'SavePlots', i, filenames);
	catch e
		%fprintf(1,'There was an error creating velocity deficit plot for file:  %s \nThe identifier was:\n%s\n', ...
		%		filenames{i},e.identifier);
		%fprintf(1,'The message was:\n%s\n',e.message);
		%fprintf(1,'Skipping file...\n\n');
	end
end

% Find indices pertaining to specific sets of conditions
files_vel_15 = ~cellfun(@isempty, strfind(filenames, '15'));
files_vel_25 = ~cellfun(@isempty, strfind(filenames, '25'));
files_cyl = ~cellfun(@isempty, strfind(filenames, 'Cylinder'));
files_af = ~cellfun(@isempty, strfind(filenames, 'Airfoil'));




%% Plot Dimensionless Velocity Deficit

% Preallocate
vel_def_dimless 	= NaN(length(A_avg_all(1,1,:)), 20);
wake_half_dimless 	= NaN(length(A_avg_all(1,1,:)), 1);
vel_def_max_dimless	= NaN(length(A_avg_all(1,1,:)), 1);

% Loop through all experiments
for i = 1:length(A_avg_all(1,1,:))
	try
		[ vel_def_dimless(i,:), wake_half_dimless(i) ,vel_def_max_dimless(i) ] ...
			= vel_deficit_dimensionless(	...
							eld_x_avg_all(:,:,i), 	...
							eld_y_avg_all(:,:,i), 	...
							v_inf_avg_all(:,:,i), 	...
							v_avg_all(:,:,i),		...
							'PlotsOff');
							%'SavePlots', i, filenames);
	catch e
		%fprintf(1,'There was an error creating velocity deficit plot for file:  %s \nThe identifier was:\n%s\n', ...
		%		filenames{i},e.identifier);
		%fprintf(1,'The message was:\n%s\n',e.message);
		%fprintf(1,'Skipping file...\n\n');
	end
end



eld_x_avg_all = squeeze(eld_x_avg_all (:,:,:));
eld_y_avg_all = squeeze(eld_y_avg_all (:,:,:));

% Plot surface of velocity deficit
x_surf = eld_x_avg_all';
y_surf = eld_y_avg_all';
z_surf = -vel_def;

x_surf = x_surf(files_vel_15 & files_af, :);
y_surf = y_surf(files_vel_15 & files_af, :);
z_surf = z_surf(files_vel_15 & files_af, :);

notNanInd = ~isnan(x_surf) & ~isnan(y_surf) & ~isnan(z_surf);
x_surf = x_surf(notNanInd);
y_surf = y_surf(notNanInd);
z_surf = z_surf(notNanInd);

f = fit( [x_surf(:), y_surf(:)], z_surf(:), 'cubicinterp' );
figure
grid on
plot(f, [x_surf, y_surf], z_surf)
set(gca,'Ydir','reverse')
title('Flow Velocity Normalized to Freestream Behind Airfoil at 15 m/s')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('Velocity Deficit (m/s)')
colormap(hot)




eld_x_avg_all = mean(eld_x_avg_all, 1);
eld_x_avg_all = eld_x_avg_all';
groupings = 2;


figure
hold on
sz = 30;
lineWidth = 1.5;

% For cylinder at 15 m/s
edgeColor = [0 .5 .5];
faceColor = [0 .7 .7];
x = eld_x_avg_all(files_vel_15 & files_cyl);
y = wake_half(files_vel_15 & files_cyl);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

% For cylinder at 25 m/s
edgeColor = [0 .2 .2];
faceColor = [0 .3 .3];
x = eld_x_avg_all(files_vel_25 & files_cyl);
y = wake_half(files_vel_25 & files_cyl);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

title('Wake Half-Width for Cylinder')
xlabel('Distance behind TE (mm)')
ylabel('Wake Half-Width (mm)')
legend('15 m/s Data', '15 m/s Fit', '25 m/s Data', '25 m/s Fit')



figure
hold on
% For airfoil at 15 m/s
edgeColor = [0 .5 .5];
faceColor = [0 .7 .7];
x = eld_x_avg_all(files_vel_15 & files_af);
y = wake_half(files_vel_15 & files_af);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

% For airfoil at 25 m/s
edgeColor = [0 .2 .2];
faceColor = [0 .3 .3];
x = eld_x_avg_all(files_vel_25 & files_af);
y = wake_half(files_vel_25 & files_af);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)



title('Wake Half-Width for Airfoil')
xlabel('Distance behind TE (mm)')
ylabel('Wake Half-Width (mm)')
legend('15 m/s Data', '15 m/s Fit', '25 m/s Data', '25 m/s Fit')





%% Plot wake max velocity deficit vs x-position

figure
hold on
sz = 30;
lineWidth = 1.5;

% For cylinder at 15 m/s
edgeColor = [.5 0 .5];
faceColor = [.9 0 .9];
x = eld_x_avg_all(files_vel_15 & files_cyl);
y = vel_def_max(files_vel_15 & files_cyl);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

% For cylinder at 25 m/s
edgeColor = [.5 0 0];
faceColor = [.7 .2 0];
x = eld_x_avg_all(files_vel_25 & files_cyl);
y = vel_def_max(files_vel_25 & files_cyl);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

title('Max Velocity Deficit for Cylinder')
xlabel('Distance behind TE (mm)')
ylabel('Velocity Deficit (m/s)')
legend('15 m/s Data', '15 m/s Fit', '25 m/s Data', '25 m/s Fit')


figure
hold on
% For airfoil at 15 m/s
edgeColor = [.5 0 .5];
faceColor = [.9 0 .9];
x = eld_x_avg_all(files_vel_15 & files_af);
y = vel_def_max(files_vel_15 & files_af);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

% For airfoil at 25 m/s
edgeColor = [.5 0 0];
faceColor = [.7 .2 0];
x = eld_x_avg_all(files_vel_25 & files_af);
y = vel_def_max(files_vel_25 & files_af);
x = arrayfun(@(i) nanmean(x(i:i+groupings-1)), ...
				1:groupings:length(x)-groupings+1)'; % code from the innerwebs that works
x = nanmean(reshape(x, [], 3),2);
y = arrayfun(@(i) nanmean(y(i:i+groupings-1)), ...
				1:groupings:length(y)-groupings+1)'; % code from the innerwebs that works
y = nanmean(reshape(y, [], 3),2);
nanInd = isnan(x) | isnan(y);
x(nanInd) = [];
y(nanInd) = [];
scatter(x, y, sz,'MarkerEdgeColor',edgeColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)
coefs = polyfit(x,y, 1);
x_fit = linspace(min(x), max(x));
y_fit = polyval(coefs, x_fit);
plot(x_fit, y_fit, 'MarkerEdgeColor',edgeColor,...
              'Color',faceColor,...
              'MarkerFaceColor',faceColor,...
              'LineWidth',lineWidth)

title('Max Velocity Deficit for Airfoil')
xlabel('Distance behind TE (mm)')
ylabel('Velocity Deficit (m/s)')
legend('15 m/s Data', '15 m/s Fit', '25 m/s Data', '25 m/s Fit')






























%% Sort data

% This data is to be categorized by experimental conditions (e.g. Cylinder, 
% V_inf = 25 m/s, x = 240 mm) and then sorted with respect to eld_y probe 
% positioning. This means that all files from experiments corresponding to 
% a specific experimental condition will be organized by their eld_y probe 
% positions. This is done in order to develop a more comprehensive picture 
% of the wake behind the test objects at various distances from their TE.

% Group files by experimental conditions

% Get indices for each experiemental condition
IND_cyl_v15_x90 = find(contains(filenames,'Cylinder') & contains(filenames,'V15') & contains(filenames,'x90'));
IND_cyl_v15_x180 = find(contains(filenames,'Cylinder') & contains(filenames,'V15') & contains(filenames,'x180'));
IND_cyl_v15_x240 = find(contains(filenames,'Cylinder') & contains(filenames,'V15') & contains(filenames,'x240'));
IND_cyl_v25_x90 = find(contains(filenames,'Cylinder') & contains(filenames,'V25') & contains(filenames,'x90'));
IND_cyl_v25_x180 = find(contains(filenames,'Cylinder') & contains(filenames,'V25') & contains(filenames,'x180'));
IND_cyl_v25_x240 = find(contains(filenames,'Cylinder') & contains(filenames,'V25') & contains(filenames,'x240'));
%--
IND_af_v15_x13 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x13'));
IND_af_v15_x18 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x18'));
IND_af_v15_x23 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x23'));
IND_af_v15_x28 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x28'));
IND_af_v15_x33 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x33'));
IND_af_v15_x38 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x38'));
IND_af_v15_x43 = find(contains(filenames,'Airfoil') & contains(filenames,'V15') & contains(filenames,'x43'));
IND_af_v25_x13 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x13'));
IND_af_v25_x18 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x18'));
IND_af_v25_x23 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x23'));
IND_af_v25_x28 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x28'));
IND_af_v25_x33 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x33'));
IND_af_v25_x38 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x38'));
IND_af_v25_x43 = find(contains(filenames,'Airfoil') & contains(filenames,'V25') & contains(filenames,'x43'));

% Get data for each experiemental condition
DATA_cyl_v15_x90 = avgData(:,:,IND_cyl_v15_x90); 
DATA_cyl_v15_x180 = avgData(:,:,IND_cyl_v15_x180); 
DATA_cyl_v15_x240 = avgData(:,:,IND_cyl_v15_x240); 
DATA_cyl_v25_x90 = avgData(:,:,IND_cyl_v25_x90); 
DATA_cyl_v25_x180 = avgData(:,:,IND_cyl_v25_x180); 
DATA_cyl_v25_x240 = avgData(:,:,IND_cyl_v25_x240); 
%--
DATA_af_v15_x13 = avgData(:,:,IND_af_v15_x13); 
DATA_af_v15_x18 = avgData(:,:,IND_af_v15_x18); 
DATA_af_v15_x23 = avgData(:,:,IND_af_v15_x23); 
DATA_af_v15_x28 = avgData(:,:,IND_af_v15_x28); 
DATA_af_v15_x33 = avgData(:,:,IND_af_v15_x33); 
DATA_af_v15_x38 = avgData(:,:,IND_af_v15_x38); 
DATA_af_v15_x43 = avgData(:,:,IND_af_v15_x43); 
DATA_af_v25_x13 = avgData(:,:,IND_af_v25_x13); 
DATA_af_v25_x18 = avgData(:,:,IND_af_v25_x18); 
DATA_af_v25_x23 = avgData(:,:,IND_af_v25_x23); 
DATA_af_v25_x28 = avgData(:,:,IND_af_v25_x28); 
DATA_af_v25_x33 = avgData(:,:,IND_af_v25_x33); 
DATA_af_v25_x38 = avgData(:,:,IND_af_v25_x38); 
DATA_af_v25_x43 = avgData(:,:,IND_af_v25_x43); 

% Sort data for each experimental condition
DATA_SORTED_cyl_v15_x90 = sortData(DATA_cyl_v15_x90);
DATA_SORTED_cyl_v15_x180 = sortData(DATA_cyl_v15_x180);
DATA_SORTED_cyl_v15_x240 = sortData(DATA_cyl_v15_x240);
DATA_SORTED_cyl_v25_x90 = sortData(DATA_cyl_v25_x90);
DATA_SORTED_cyl_v25_x180 = sortData(DATA_cyl_v25_x180);
DATA_SORTED_cyl_v25_x240 = sortData(DATA_cyl_v25_x240);
%--
DATA_SORTED_af_v15_x13 = sortData(DATA_af_v15_x13);
DATA_SORTED_af_v15_x18 = sortData(DATA_af_v15_x18);
DATA_SORTED_af_v15_x23 = sortData(DATA_af_v15_x23);
DATA_SORTED_af_v15_x28 = sortData(DATA_af_v15_x28);
DATA_SORTED_af_v15_x33 = sortData(DATA_af_v15_x33);
DATA_SORTED_af_v15_x38 = sortData(DATA_af_v15_x38);
DATA_SORTED_af_v15_x43 = sortData(DATA_af_v15_x43);
DATA_SORTED_af_v25_x13 = sortData(DATA_af_v25_x13);
DATA_SORTED_af_v25_x18 = sortData(DATA_af_v25_x18);
DATA_SORTED_af_v25_x23 = sortData(DATA_af_v25_x23);
DATA_SORTED_af_v25_x28 = sortData(DATA_af_v25_x28);
DATA_SORTED_af_v25_x33 = sortData(DATA_af_v25_x33);
DATA_SORTED_af_v25_x38 = sortData(DATA_af_v25_x38);
DATA_SORTED_af_v25_x43 = sortData(DATA_af_v25_x43);

% Plotting eld_y for demonstation
%{
figure
plot(DATA_cyl_v15_x90(:,3),'o')
figure
plot(DATA_SORTED_cyl_v15_x90(:,3),'o')
figure
plot(DATA_SORTED_af_v25_x18(:,3), 'o')
%}



%{
eldY = DATA_SORTED_af_v15_x13(:,3);
Vinf = DATA_SORTED_af_v15_x13(:,end-1);
V = DATA_SORTED_af_v15_x13(:,end);

figure
scatter(eldY, Vinf - V)
%}
