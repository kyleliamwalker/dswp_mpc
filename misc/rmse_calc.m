% Author: Dr. Kyle L. Walker
% Description: Calculates Root Mean Square Error (RMSE)

function [ rmse ] = rmse_calc( data, target )

sq_errors =( data - target' ).^2;

rmse = sqrt(mean(sq_errors));

end