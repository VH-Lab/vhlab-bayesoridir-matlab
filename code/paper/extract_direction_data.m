function [dirinfo] = extract_direction_data(dir_doc)
% EXTRACT_DIR_DATA - extracts direction tuning curves for each SF, TF
%
% DIRINFO = EXTRACT_DIR_DATA(DIR_DOC)
%
% Returns a set of DIRINFO structures with the following fields:
%   sf - the spatial frequency used
%   tf - the temporal frequency used
%   angle - the set of angles examined
%   mean_responses - the mean responses
%  


independent_variable_value = dir_doc.document_properties.stimulus_tuningcurve.independent_variable_value;

dirinfo.sf = [];
dirinfo.tf = [];
dirinfo.angle = [];
dirinfo.mean_responses = [];

dirinfo = dirinfo([]);

sfs = unique(independent_variable_value(:,1));
tfs = unique(independent_variable_value(:,2));

for s = 1:numel(sfs),
    for t = 1:numel(tfs),
        my_responses = find(independent_variable_value(:,1)==sfs(s) & independent_variable_value(:,2)==tfs(t));
        clear mydirinfo;
        mydirinfo.sf = sfs(s);
        mydirinfo.tf = tfs(t);
        angle_values = independent_variable_value(my_responses,3);
        [angles_sorted, angle_indexes] = sort(angle_values);
        mydirinfo.angle = angles_sorted;
        mydirinfo.mean_responses = dir_doc.document_properties.stimulus_tuningcurve.response_mean(my_responses(angle_indexes));
        dirinfo(end+1) = mydirinfo;
    end;
end;

