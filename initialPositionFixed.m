function position = initialPositionFixed(userParams)
% initialPositionFixed: Calculates the initial position of a user based on given parameters.
%
% Input:
%   userParams : A structure containing user parameters such as maximum angle, maximum distance, etc.
%
% Output:
%   position   : A structure containing the calculated initial position of the user.

    % Retrieve constants from myPackageConstant function
    c = myPackageConstant();
    
    % Initialize the position structure
    position = struct([]);

    % Use a random number generator with shuffling
    rng('shuffle');

    % Randomly select a sign for the angle
    sign = [-1, 1];
    i = randi(2);
    s = sign(i);

    % Calculate initial position details
    position(1).theta = s * userParams.mAngle * rand(1); % Initial angle (in degrees)
    position(1).x = userParams.maxDist * cosd(position.theta); % x-coordinate
    position(1).y = userParams.maxDist * sind(position.theta); % y-coordinate
    position(1).distance = sqrt(position.x^2 + position.y^2); % 2D distance from the origin
    position(1).distance3d = sqrt(position.distance^2 + (c.txHeight - c.rxHeight)^2); % 3D distance considering heights
    
    rng('shuffle'); % Shuffle the random number generator again
    position(1).phi = 360 * rand(1); % Random azimuth angle (in degrees)
    position(1).d1 = 0; % Initialize another parameter (purpose should be documented)

end
