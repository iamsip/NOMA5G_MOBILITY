function c = myPackageConstant()
% myPackageConstant: A function that initializes and returns a structure containing various constants for the package.
%
% Output:
%   c : A structure containing various predefined constants.

    % Initialize the structure 'c' with various constants
    c = struct;

    c.time = 20;       % Time-related constant, perhaps duration (the unit should be defined, e.g., seconds)
    c.freq = 1.8;      % Frequency in GHz
    c.txHeight = 25;   % Transmitter height (the unit should be defined, e.g., meters)
    c.scenario = "UMa"; % Scenario type (RMa, UMa, UMi are possible values; specifics should be documented)
    c.rxHeight = 1.5;  % Receiver height (the unit should be defined, e.g., meters)

end
