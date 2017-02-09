function [ data_median, data_outliers, data_number_of_outliers ] = calculate_median_and_outliers( data, cutoff, tile_ids, determination_method )
data_median = median(data);
data_outliers = [];
if strcmp(determination_method, 'fixed_cutoff')
    data_outliers = tile_ids(abs(data-1)>=cutoff); % Cutoff is the percent varia
elseif strcmp(determination_method, 'std')
    data_outliers = tile_ids(abs(data-mean(data))>cutoff*std(data)); % Cutoff is the number of standard deviations
end
data_number_of_outliers = numel(data_outliers);
end

