

function [ rmse ] = rmse_calc( data, target )

sq_errors =( data - target' ).^2;

rmse = sqrt(mean(sq_errors));

end