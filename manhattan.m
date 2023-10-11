
% MANHATTAN Mobility model replicating the Manhattan movement in urban areas.
% Simulates grid-based city layout with right-angle turns, characteristic of 
% most urban areas.
%
% INPUTS:
% position - Structure with current position and related mobility metrics.
% rxHeight - Height of the receiver node/device.
%
% OUTPUTS:
% position - Updated structure with new positional information post-movement.
function [position] = manhattan(position, rxHeight)
  % Get constants from package
  c = myPackageConstant();
  
  % Define minimum and maximum speeds for urban movement (in kph)
  speedMin = 11;
  speedMax = 40;
  
  % Randomly select a speed within the defined range
  speed = speedMin + (speedMax-speedMin) * rand(1);
  
  % Calculate the distance covered based on the speed and the time interval
  d = (speed * c.time / 3600) * 1000;  % in meters
  
  % Update the internal distance counter
  position.d1 = position.d1 + d;
  
  % Check if the object has traveled more than a set distance (e.g., 50 meters)
  if (position.d1 > 50)
      rng('shuffle')
      
      % Assign one of the four right-angle directions (0, 90, 180, 270 degrees)
      n = randi([1,4],1);
      position.phi = 90 * n;        % angle in degree
      
      % Reset the internal distance counter
      position.d1 = 0;
  end
  
  % Update x and y coordinates based on the chosen direction and distance traveled
  position.x = position.x + d * cosd(position.phi);
  position.y = position.y + d * sind(position.phi);
  
  % Compute the 2D distance from origin
  position.distance = sqrt(position.x^2 + position.y^2);
  
  % Compute the 3D distance considering the height difference
  position.distance3d = sqrt(position.distance^2 + (c.txHeight - rxHeight)^2);
  
  % Compute the angle based on the updated positions
  position.theta = acosd((position.distance^2 + position.x^2 - position.y^2) / (2 * position.distance * position.x));

end
