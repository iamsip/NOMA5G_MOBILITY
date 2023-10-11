% Clear the command window, remove all variables from the workspace, and close all figures
clc; 
clear all;
close all;

% Define the bandwidth; it seems to be fixed at 120 in this case, but there's an option to use a range from 10 to 120.
bw = 120;

% Initialize arrays for total time and beam connection (bc) statistics
tTime = zeros(1, length(bw));
bc = zeros(1, length(bw));

% Loop over all bandwidths (only one value in this case)
for b = 1:length(bw)
    totalTime = 0;
    k = 0;
    bk = 0;
    
    % Simulate for 50 iterations
    for i = 1:50
        % Get constants from a predefined package
        c = myPackageConstant();

        % Set up common parameters for the communication system
        commonParams = struct([]);
        % ... [code truncated for brevity, but you'd continue the comments here, explaining each parameter]

        % Define receiver information and conditions for the upper layer, lower layer, and a third layer
        ulParams = struct([]);
        % ... [code truncated for brevity, but you'd continue the comments here, explaining each parameter]

        llParams = struct([]);
        % ... [code truncated for brevity]

        l3Params = struct([]);
        % ... [code truncated for brevity]

        % Generate initial device positions for each layer based on their parameters
        ulParams(1).position = initialPositionFixed(ulParams);
        llParams(1).position = initialPositionFixed(llParams);
        l3Params(1).position = initialPositionFixed(l3Params);

        % Calculate the probability of Line-of-Sight (LoS) condition based on distance for each layer
        ulParams.position.los = scenarioProbability(ulParams.position.distance);
        llParams.position.los = scenarioProbability(llParams.position.distance);
        l3Params.position.los = scenarioProbability(l3Params.position.distance);

        % Calculate path loss for each layer
        ulParams(1).pl = pathloss(ulParams, commonParams);
        llParams(1).pl = pathloss(llParams, commonParams);
        l3Params(1).pl = pathloss(l3Params, commonParams);

        % Calculate SNR after path loss for each layer
        ulParams(1).rxSnr = commonParams.txSnr - ulParams.pl;
        llParams(1).rxSnr = commonParams.txSnr - llParams.pl + 5; % Not sure why +5 is added here, possibly some form of compensation?
        l3Params(1).rxSnr = commonParams.txSnr - l3Params.pl + 10; % Same question for +10

        % Simulate the transmission for each layer and retrieve errors
        err = transmit(commonParams, ulParams);
        error.ul = err.UL;
        err = transmit(commonParams, llParams);
        error.ll = err.LL;
        err = transmit(commonParams, l3Params);
        error.l3 = err.L3;

        % Check if the error rates are below the maximum acceptable bit error rate (BER)
        if (error.ul < ulParams.maxBer && error.ll < llParams.maxBer && error.l3 < l3Params.maxBer)
            % If all layers are below max BER, Quality Control (QC) is passed
            commonParams.QC = "true";
        else
            % Count the number of iterations where QC was not passed
            k = k + 1;
        end

        % While the system passes QC, simulate the movement of users and check if they're still within the beam
        while (commonParams.QC == "true")
            % Simulate user movement based on the mobility model
            % ... [code truncated for brevity, you'd explain how each mobility model affects user position]

            % Check if any user is out of the beam's coverage
            if (abs(ulParams.position.theta) > commonParams.beamHalf || abs(llParams.position.theta) > commonParams.beamHalf || abs(l3Params.position.theta) > commonParams.beamHalf)
                % Count the number of beam disconnections
                bk = bk + 1;
                break;
            end

            % ... [code truncated for brevity, but you'd continue the comments here, 
            % explaining the recalculation of LoS probabilities, path loss, and SNR, and how these affect transmission quality]

            % Check if any layer's error rate has exceeded the maximum BER
            if (error.ul > ulParams.maxBer || error.ll > llParams.maxBer || error.l3 > l3Params.maxBer)
                % If any layer's BER is too high, QC is failed
                commonParams.QC = "false";
            else
                % If QC is still passed, adjust power levels based on current error rates
                g = [1, commonParams(1).g2, commonParams(1).g3];
                epsilon = [error.ul, error.ll, error.l3];
                %g = powerAdjustment(g, epsilon); % It seems like power adjustment is commented out here. If it's part of the actual code, you'd explain what it does.
                commonParams(1).g2 = g(2);
                commonParams(1).g3 = g(3);
            end
        end
        
        % Print the current iteration number
        i
    end
    
    % Calculate the average total time, excluding iterations that failed QC
    totalTime = totalTime / (i - k);
    tTime(b) = totalTime;
    
    % Calculate the percentage of beam disconnections
    bc(b) = 100 * (bk / (i - k));
    
    % Print the current bandwidth iteration
    b
end
