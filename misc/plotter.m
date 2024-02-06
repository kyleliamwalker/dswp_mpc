% Author: Dr. Kyle L. Walker
% Description: Plots the robots positional evolution in time and updates
% all positional errors at the same time.

function [ robot ] = plotter( robot, count )

Y = [ robot.state.x, robot.state.z, robot.state.theta, ...
    robot.state.u, robot.state.w, robot.state.q, ...
    robot.state.udot, robot.state.wdot, robot.state.qdot ]; 

[ robot.robotPlots ] = updatePlotHistory( Y, robot.robotPlots, count, 1 );

[ robot ] = updateErrors( robot, count, robot.errors.pError );

end

