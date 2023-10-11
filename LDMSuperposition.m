function [ldmData] = LDMSuperposition(upperLayerData, lowerLayerData, injectionLevel)
% LDMSuperposition: This function performs the superposition of two signal layers, typically used in a 
% NOMA (Non-Orthogonal Multiple Access) communication system. The lower layer data is scaled by an 
% 'injection level,' which is a power ratio factor before it's combined with the upper layer data. 
% The process forms a composite signal that carries the information from both layers.
%
% Inputs:
%   upperLayerData  : Signal data of the upper layer. It represents one of the two NOMA layers 
%                     (usually the layer intended for the user with a better channel condition).
%                     Expected to be a scalar or a vector.
%   lowerLayerData  : Signal data of the lower layer. It represents the other NOMA layer 
%                     (usually the layer intended for the user with a poorer channel condition).
%                     Expected to be a scalar or a vector.
%   injectionLevel  : The power ratio used to scale the lower layer data before superposition. 
%                     It's a dimensionless quantity and represents the power allocation 
%                     for the lower layer's signal. Expected to be a scalar.
%
% Output:
%   ldmData         : The resultant data after the superposition of the upper and lower layers. 
%                     This composite signal contains the information from both the upper and 
%                     lower layers, with the lower layer data scaled according to the injection level.

% Scale the lower layer data by the injection level. The injection level represents the power ratio 
% allocated to the lower layer's signal. This scaling is crucial for maintaining the NOMA principle 
% of superposing signals with different power levels.
lowerLayerData = lowerLayerData .* injectionLevel;

% Create the superposed signal by linearly combining the scaled lower layer data with the upper layer data. 
% The resultant 'ldmData' is a composite signal that carries both layers' information, ready for transmission 
% over a shared communication medium.
ldmData = lowerLayerData + upperLayerData;
end

