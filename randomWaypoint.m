
% RANDOMWAYPOINT Simulates mobility using the Random Waypoint model.
% 
% INPUT:
%   position: Current position of the node (structure with fields x, y, d1, phi, etc.)
%   rxHeight: Height of the receiver
%
% OUTPUT:
%   position: Updated position of the node after movement.
function [position]=randomWaypoint(position, rxHeight)
    c = myPackageConstant();
    
    % Random speed between speedMin and speedMax
    speedMin = 0;
    speedMax = 10;
    speed = speedMin + (speedMax-speedMin)*rand(1); % kph
    
    % Calculate the distance moved
    d = (speed*c.time/3600)*1000;  % in meters
    
    position.d1 = position.d1 + d;
    
    % If node goes beyond a boundary, reset its angle randomly
    if (position.d1 > 50)
        rng('shuffle')
        position.phi = 360*rand(1);  % Angle in degrees
        position.d1 = 0;
    end
    
    % Update x and y positions based on movement direction
    position.x = position.x + d*cosd(position.phi);
    position.y = position.y + d*sind(position.phi);
    
    % Update distance and 3D distance attributes
    position.distance = sqrt(position.x^2 + position.y^2);
    position.distance3d = sqrt(position.distance^2 + (c.txHeight-rxHeight)^2);
    
    % Compute the angle theta
    position.theta = acosd((position.distance^2 + position.x^2 - position.y^2) / (2*position.distance*position.x));

end
