function position = initialPosition(userParams)
% initialPosition: Calculate the initial random position of a user within defined constraints.
%
% Input:
%   userParams : Structure with user parameters that includes maximum and minimum distances, 
%                maximum angle, and receiver height.
%
% Output:
%   position   : Structure with the calculated initial position details of the user.

    % Retrieve constants from myPackageConstant function
    c = myPackageConstant();
    
    % Initialize the position structure
    position = struct([]);

    % Shuffle the random number generator for unpredictability
    rng('shuffle');

    % Calculate random theta within the maximum angle
    rTheta = userParams.mAngle * rand(1);

    % Calculate maximum and minimum x, y coordinates based on random theta
    xMax = cosd(rTheta) * userParams.maxDist;
    yMax = sind(rTheta) * userParams.maxDist;
    xMin = cosd(rTheta) * userParams.minDist;
    yMin = sind(rTheta) * userParams.minDist;

    % Choose a random sign for y-coordinate and theta
    sign = [-1,1];
    i = randi(2);
    s = sign(i);

    % Shuffle the random number generator again
    rng('shuffle');

    % Randomly select x within its min and max range
    position(1).x = xMin + rand(1) * (xMax - xMin);

    % Shuffle the random number generator and randomly select y within its min and max range, apply sign
    rng('shuffle');
    position(1).y = (yMin + rand(1) * (yMax - yMin)) * s;
    position(1).theta = rTheta * s; % Apply sign to theta

    % Calculate 2D and 3D distances
    position(1).distance = sqrt(position.x^2 + position.y^2);
    position(1).distance3d = sqrt(position.distance^2 + (c.txHeight - userParams.rxHeight)^2);

    % Generate a random azimuth angle (phi) in degrees
    rng('shuffle');
    position(1).phi = 360 * rand(1);

    % Initialize an additional parameter d1
    position(1).d1 = 0;

end
