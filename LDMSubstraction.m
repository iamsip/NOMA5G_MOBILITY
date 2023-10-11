function [rLDM] = LDMSubstraction(rxLDM, rxUpper, injectionLevel)
% LDMSubstraction: This function calculates the reduced LDM signal by 
% subtracting the upper received signal and normalizing with the injection level.
%
% Inputs:
%   rxLDM         : The received LDM (Lower Density Multiplexing) signal, expected to be a vector.
%   rxUpper       : The received upper layer signal, expected to be a scalar or vector of the same length as rxLDM.
%   injectionLevel: The level of signal injection, used for normalization, expected to be a scalar.
%
% Output:
%   rLDM          : The reduced LDM signal after subtraction and normalization.

% Reshape rxLDM to a column vector to ensure consistency in the subtraction operation.
% This step is necessary to avoid dimension mismatch errors if rxLDM is not originally in column format.
rxLDM = reshape(rxLDM, [length(rxLDM), 1]);

% Subtract the received upper signal from the received LDM signal.
% This step assumes that rxUpper is broadcasted (if scalar) or has the same length as rxLDM.
% The subtraction models the interference or signal removal process in the receiver.
rLDM = rxLDM - rxUpper;

% Normalize the subtracted signal by the injection level.
% This step scales the signal based on the predefined injection level, affecting the final signal strength or quality.
% The division operation assumes injectionLevel is a non-zero scalar.
rLDM = (rLDM ./ injectionLevel);
end
