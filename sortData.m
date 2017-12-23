function DATA_SORTED = sortData(unsorted_data_files)
%SORTDATA  Sorts data based on the eld_y probe position. This function 
%  		   takes in a 3D matrix, reshapes it, and sorts based on eld_y.
%  
%  Inputs:		unsorted_data_files	- 3D matrix where each experiment's 
% 									  data is separated in the 3rd dim.
%
%  Outputs:		DATA_SORTED			- 2D matrix of sorted data
%  
%  Created by:     Keith Covington
%  Created on:     10/02/2017
%  Last modified:  10/02/2017
% *************************************************************************


dataSz = size(unsorted_data_files);

DATA_UNSORT = [];
for i = 1:dataSz(3)
	DATA_UNSORT = cat(1, DATA_UNSORT, ...
						 unsorted_data_files(:,:,i));
end

[~,sortKey] = sort(DATA_UNSORT(:,3)); % create a sort order defined by the eld_y probe
DATA_SORTED= DATA_UNSORT(sortKey, :); % sort data based on sortKey



