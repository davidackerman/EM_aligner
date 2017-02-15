function [ data_median, data_outlier_count, data_outlier_indices, data_outlier_tile_ids ] = calculate_median_and_outliers( data, cutoff, tile_ids, determination_method, only_greater_than )
if nargin < 5
    only_greater_than = false; % Used if only want outliers greater than mean, for eg. residuals where smaller than the mean is better
end
data_median = nanmedian(data);
data_outlier_tile_ids = [];
if strcmp(determination_method, 'fixed_cutoff')
    if only_greater_than
        data_outlier_indices = find(data-1>=cutoff); % Cutoff is the accepted fixed deviation magnitude
    else
        data_outlier_indices = find(abs(data-1)>=cutoff);
    end
    data_outlier_tile_ids = tile_ids(data_outlier_indices);
elseif strcmp(determination_method, 'std')
    if only_greater_than
        data_outlier_indices = find(data-nanmean(data)>=cutoff*nanstd(data)); % Cutoff is the number of standard deviations, but only care about positive
    else
        data_outlier_indices = find(abs(data-nanmean(data))>=cutoff*nanstd(data)); % Cutoff is the number of standard deviations, but only care about positive
    end
    data_outlier_tile_ids = tile_ids(data_outlier_indices); 
end
data_outlier_count = numel(data_outlier_tile_ids);
end

