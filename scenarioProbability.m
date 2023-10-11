% Function to calculate Line-of-Sight (LOS) probability in different scenarios
function [los] = scenarioProbability (d)
    % Get constants for the current scenario
    c = myPackageConstant();

    % Initialize probability variable
    pr = 0; % It's good practice to initialize variables to avoid potential issues in code execution

    % Calculate the LOS probability based on the scenario
    switch (c.scenario)
        case("RMa") % For Rural Macro-cell scenario
            % Probability decreases exponentially as distance increases beyond 10 meters
            pr = exp(-(d - 10) / 1000);

        case("UMa") % For Urban Macro-cell scenario
        case("UMi") % For Urban Micro-cell scenario
            % Both UMa and UMi scenarios share the same LOS probability calculation
            % Probability has a base value inversely proportional to distance, and decreases exponentially with a base dependent on distance
            pr = (18 / d) + exp(-d / 63) * (1 - (18 / d));
    end

    % Ensure probability doesn't exceed 1
    if pr > 1
        pr = 1;
    end

    % Calculate if LOS exists based on the probability
    los = losProbability(pr);
end
