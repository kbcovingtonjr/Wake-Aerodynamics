function [vel_def,wake_half,vel_def_max ] = vel_deficit(eld_x, eld_y,v_inf,v,varargin)
%VEL_DEFICIT  calculates the velocity deficit at each y-location in the 
%  wind tunnel as well as the corresponding half wake at half maximum. 
%  Returns these values and optionally plots and/or saves plots as a jpeg.
%  
%  Inputs:		eld_x	- x-position in wind tunnel test section
%  				eld_y	- y-position in wind tunnel test section
%  				v_inf	- freestream velocity
%  				v		- measured velocity at (x,y) position
%  				'PlotsOff' [optional]			- do not display
%  				'SavePlots',expNum [optional]	- save plot as jpeg with
%  												  experiment number label
%  Outputs:		c_l     - sectional coefficient of lift
%  
%  
%  Example calls:
%  	[ vel_def, wake_half ,vel_def_max ] = vel_deficit(	eld_xData,	...
%  														eld_yData,	...
%  														v_infData,	...
%  														v_Data);
%  	[ vel_def, wake_half ,vel_def_max ] = vel_deficit(	eld_xData,	...
%  														eld_yData,	...
%  														v_infData,	...
%  														v_Data,		...
%  														'PlotsOff');
%  	[ vel_def, wake_half ,vel_def_max ] = vel_deficit(	eld_xData,	...
%  														eld_yData,	...
%  														v_infData,	...
%  														v_Data,		...
%  														'SavePlots', expNum);
%  	[ vel_def, wake_half ,vel_def_max ] = vel_deficit(	eld_xData,	...
%  														eld_yData,	...
%  														v_infData,	...
%  														v_Data,		...
%  														'PlotsOff',	...
%  														'SavePlots', expNum);
%  
%  Created by:		Ginger Beerman
%  Modified by:		Keith Covington
%  Created on:		09/08/2017
%  Last modified:	11/02/2017
% *************************************************************************


% Defaults
plotsOn = true;
savePlots = false;

% Change defaults if input arguments say to do so
if ~isempty(varargin)
	if any( strcmp('PlotsOff', varargin) )
		plotsOn = false;
	end
	if any( strcmp('SavePlots', varargin) )
		savePlots = true;
		ind = find(strcmp('SavePlots', varargin));
		expNum = varargin{ind+1};
		filenames = varargin{ind+2};
	end
end


%% Velocity Deficit vs Vertical Location y
vel_def= v_inf - v;%find velocity deficit
vel_def = vel_def';
eld_spline=linspace(min(eld_y),max(eld_y),5000);
s=spline(eld_y,vel_def,eld_spline);

%% Wake Half-width vs horizontal location X 

vel_def_max=max(s);
half_vel_def=0.5*max(s);
ind1=find(s >= half_vel_def,1);
ind2=find(s >= half_vel_def,1,'last');

y1=eld_spline(ind1);
y2=eld_spline(ind2);
wake_half=abs(.5*(y1-y2));

if plotsOn || savePlots

	set(0, 'defaulttextInterpreter', 'latex') % plotting necessities

	if plotsOn
		fig = figure;
	else % if only saving but not showing figure
		fig = figure('visible', 'off');
	end

	scatter(vel_def, eld_y)%plot vel deficit
	hold on;
	filename = filenames{expNum};
	%section = filename( strfind(filename,'Section')+7 )
	velocity = filename( strfind(filename,'V')+1 : strfind(filename,'V')+2 );
	xPosition = filename( strfind(filename,'_x')+2 : strfind(filename,'_x')+4 );
	if xPosition(end) == '_'
		xPosition(end) = [];
	end
	if any(strfind(filename,'Airfoil'));
		objectType = 'Airfoil';
	elseif any(strfind(filename,'Cylinder'));
		objectType = 'Cylinder';
	else
		warning('Object type not found in filename.')
	end

	title([ 'Velocity Deficit For ' objectType ' at ' velocity ' m/s, ' ...
			xPosition ' mm behind TE']);
	ylabel('ELD vertical location (mm)');
	xlabel('Velocity Deficit (m/s)');
	plot(s,eld_spline);
	%plot(s(ind1),y1,'*');
	%plot(s(ind2),y2,'*');
	plot([s(ind1) s(ind2)], [y1 y2])
	legend('Data','Fit','Full Wake at Half Maximum');

	if savePlots
		print( [ 'VelDeficit_Exp' num2str(expNum) ], '-djpeg') % save figure as jpeg TODO: save to specific folder
		if ~plotsOn % if only saving but not showing figure
			close(fig) % close figure so we don't fill your computer's RAM and swap and eventually kill it
		end
	end
end


end

