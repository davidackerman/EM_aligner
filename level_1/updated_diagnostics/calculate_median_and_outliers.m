function [ data_median, data_outliers, data_number_of_outliers ] = calculate_median_and_outliers( data, number_of_standard_deviations, tile_ids )
data_median = median(data);
data_outliers = tile_ids(abs(data-mean(data))>number_of_standard_deviations*std(data));
data_number_of_outliers = numel(data_outliers);
end

