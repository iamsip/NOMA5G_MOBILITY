function [ldmData] = LDMSuperpositionN(L1, L2, L3, L4, L5, g2, g3, g4, g5)
% LDMSuperpositionN: This function performs a weighted superposition of multiple signal layers.
% The weights are the gains applied to each signal layer (L2 to L5), and the superposition
% is a linear combination of all these layers including L1. The function is used for
% creating a composite signal from multiple layers in a NOMA (Non-Orthogonal Multiple Access) system.
%
% Inputs:
%   L1, L2, L3, L4, L5: These are the signal layers, where each Lx represents a distinct signal layer.
%                       Each Lx is expected to be a scalar or a vector.
%   g2, g3, g4, g5    : These are the gain factors applied to the signal layers L2 to L5, respectively.
%                       These gains are dimensionless and determine the contribution of each 
%                       layer to the final superposed signal. Each gx is expected to be a scalar.
%
% Output:
%   ldmData           : The resultant data after the superposition of the five layers. 
%                       It's a linear combination of the input layers, each weighted by its respective gain.

% Calculate the superposed signal by adding each layer L2 to L5, each scaled by its respective 
% product of gains, to L1. This forms a composite signal from multiple layers, with each layer 
% contributing to the final signal proportional to its gain factors.
ldmData = L5.*g2*g3*g4*g5 + L4.*g2*g3*g4 + L3.*g2*g3 + L2.*g2 + L1;
end
