% LOSPROBABILITY - Determines if there's Line-of-Sight (LOS) based on a given probability
%
% This function takes in a probability value `pr` representing the likelihood
% of having a LOS connectivity. The function then returns a binary decision (1 or 0)
% indicating the presence (1) or absence (0) of LOS, based on the provided probability.
%
% Syntax:  los = losProbability(pr)
%
% Inputs:
%    pr - Probability of having a LOS connectivity (range: 0-1).
%
% Outputs:
%    los - Binary decision indicating LOS presence (1) or absence (0).
%
% Example: 
%    result = losProbability(0.7)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Md Shantanu Islam
% Work address
% email: your.email@example.com
% Website: https://www.example.com
% Created: Date using MATLAB R20XXx 
% Revisions: 

function [los] = losProbability(pr)
    % Total number of events for sampling
    noOfEvent = 10000;
    
    % Calculate the number of entries representing LOS
    losEntries = floor(pr * noOfEvent);
    
    % Initialize an array for N-LOS events
    nLosArrays = zeros(1, noOfEvent);
    
    % Randomly assign the LOS entries based on the provided probability
    entryPosition = randperm(noOfEvent);
    for i = 1:1:losEntries
        n = entryPosition(i);
        nLosArrays(n) = 1;
    end
    
    % Randomly make a decision based on the constructed probability array
    makeDecision = randi([1 noOfEvent], 1, 1);
    los = nLosArrays(makeDecision);
end
