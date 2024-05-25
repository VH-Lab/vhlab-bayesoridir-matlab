clear all; close all;clc;
%%
%To check which cells exhibit strong responses before motion exposure, use:

load stevesolddata.mat

good_indexes_before = [];
good_indexes_motionexposure = [];

for i=1:numel(cell_structures),
	for j=1:numel(cell_structures{i}),
		if strcmp(cell_structures{i}(j).type,'TP Ach OT vec varies p'),
			if cell_structures{i}(j).data<0.05,
				good_indexes_before(end+1) = i;
			end;
		end;
		if strcmp(cell_structures{i}(j).type,'TP ME Ach OT vec varies p'),
			if cell_structures{i}(j).data<0.05,
				good_indexes_motionexposure(end+1) = i;
			end;
		end;
	end;
end;
%%
% To extract direction tuning data for each cell before motion exposure, use:

for i=1:numel(good_indexes_before),
	for j=1:numel(cell_structures{good_indexes_before(i)}),
		if strcmp(cell_structures{good_indexes_before(i)}(j).type,'TP Ach OT Response curve'),
			data_before{i} = cell_structures{good_indexes_before(i)}(j).data,
		end;
	end;
end;
%% 
% To extract direction tuning data for each cell after motion exposure, use:

for i=1:numel(good_indexes_motionexposure),
	for j=1:numel(cell_structures{good_indexes_motionexposure(i)}),
		if strcmp(cell_structures{good_indexes_motionexposure(i)}(j).type,'TP ME Ach OT Response curve'),
			data_motionexposure{i} = cell_structures{good_indexes_motionexposure(i)}(j).data;
		end;
	end;
end;
%% 

i = 100;

for j=1:numel(cell_structures{good_indexes_before(i)}),
    if strcmp(cell_structures{good_indexes_before(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
        i,j,
        data = cell_structures{good_indexes_before(i)}(j).data,
    end;
end;

