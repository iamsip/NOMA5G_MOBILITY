clc; % Clear Command Window
clear variables; % Clear variables from the workspace
tic % Start stopwatch timer

%% Scenario Information %%
c = myPackageConstant(); % Load constants from a custom package
f=1:0.5:7; % Frequency range initialization
f=1; % Overwriting frequency value for this simulation
bw=10:5:120; % Bandwidth range initialization
% Time in seconds
totalTime=zeros(1,length(bw)); % Initializing totalTime array with zeros

% Loop over different bandwidths
for i=1:length(bw)
    k=0; % Counter initialization
    for j=1:1 % Loop runs once (can be extended for multiple iterations)

        % Initializing common parameters for the scenario
        commonParams=struct([]);
        commonParams(1).slots=1;
        commonParams(1).nTxAnts=8; % Number of transmit antennas
        commonParams(1).channel=nrTDLChannel; % Specifying the channel model
        commonParams(1).delayProfile="TDL-C"; % Delay profile for the channel
        commonParams(1).g=-5; % Gain
        commonParams(1).txSnr=100; % Transmit SNR
        commonParams(1).modulationUL="16QAM"; % Modulation scheme for Upper Layer
        commonParams(1).modulationLL="64QAM"; % Modulation scheme for Lower Layer
        commonParams(1).QC="false"; % Quick convergence flag (QC)
        commonParams(1).beamHalf=bw(i)/2; % Half of the beamwidth
        commonParams(1).freq=f(i); % Frequency in GHz

        %% Receiver Information %%

        % Upper Layer (UL) receiver parameters
        ulParams=struct([]);
        ulParams(1).nRxAnts=2; % Number of receive antennas
        ulParams(1).maxDist=100; % Maximum distance from transmitter
        ulParams(1).minDist=10; % Minimum distance from transmitter
        ulParams(1).mAngle=30; % Main angle
        ulParams(1).maxBer=1e-4; % Maximum Bit Error Rate (BER)
        ulParams(1).mobilityModel="manhattan"; % Mobility model (wayPoint or manhattan)
        ulParams(1).rxHeight=1.5; % Receiver height

        % Lower Layer (LL) receiver parameters
        llParams=struct([]);
        llParams(1).nRxAnts=8; % Number of receive antennas
        llParams(1).maxDist=50; % Maximum distance from transmitter
        llParams(1).minDist=10; % Minimum distance from transmitter
        llParams(1).mAngle=15; % Main angle
        llParams(1).maxBer=1e-2; % Maximum Bit Error Rate (BER)
        llParams(1).mobilityModel="manhattan"; % Mobility model (wayPoint or manhattan)
        llParams(1).rxHeight=1.5; % Receiver height

        %% Generating initial device position %%

        % Initial position for UL and LL based on their respective parameters
        ulParams(1).position=initialPositionFixed(ulParams);
        llParams(1).position=initialPositionFixed(llParams);

        % Determine Line of Sight (LoS) probability based on distance
        ulParams.position.los=scenarioProbability(ulParams.position.distance);
        llParams.position.los=scenarioProbability(llParams.position.distance);

        % Calculate pathloss for both UL and LL
        ulParams(1).pl=pathloss(ulParams,commonParams);
        llParams(1).pl=pathloss(llParams,commonParams);

        % Calculate receive SNR for UL and LL after pathloss
        ulParams(1).rxSnr=commonParams.txSnr-ulParams.pl;
        llParams(1).rxSnr=commonParams.txSnr-llParams.pl;

        % Transmit and receive process for UL and LL, and check for errors
        err=transmit(commonParams,ulParams); % Transmitting UL
        error.ul=err.UL; % Error for UL
        err=transmit(commonParams,llParams); % Transmitting LL
        error.ll=err.LL; % Error for LL

        % Activate Quick Convergence (QC) if error is below threshold
        if (error.ul<ulParams.maxBer || error.ll<llParams.maxBer)
            commonParams.QC="true";
        end

        % Continue the process until Quick Convergence (QC) is achieved
        while (commonParams.QC=="true")
            %% Calling Mobility Function %%
            % Update position based on mobility model
            switch (ulParams.mobilityModel)
                case("wayPoint")
                    ulParams.position=randomWaypoint(ulParams.position, ulParams.rxHeight); % Random waypoint mobility model
                case("manhattan")
                    ulParams.position=manhattun(ulParams.position, ulParams.rxHeight); % Manhattan mobility model
            end

            switch (llParams.mobilityModel)
                case("wayPoint")
                    llParams.position=randomWaypoint(llParams.position, llParams.rxHeight); % Random waypoint mobility model
                case("manhattan")
                    llParams.position=manhattun(llParams.position, llParams.rxHeight); % Manhattan mobility model
            end

            % Check if device is out of beam, if so, break the loop
            if(abs(ulParams.position.theta)>commonParams.beamHalf || abs(llParams.position.theta)>commonParams.beamHalf)
                k=k+1; % Increment counter
                break; % Break the loop as device is out of beam
            end

            % Recalculate LoS probability and pathloss
            ulParams.position.los=scenarioProbability(ulParams.position.distance);
            llParams.position.los=scenarioProbability(llParams.position.distance);
            ulParams.pl=pathloss(ulParams,commonParams);
            llParams.pl=pathloss(llParams,commonParams);

            % Recalculate receive SNR for UL and LL after pathloss
            ulParams.rxSnr=commonParams.txSnr-ulParams.pl;
            llParams.rxSnr=commonParams.txSnr-llParams.pl;

            % Transmit and receive process for UL and LL, and check for errors
            err=transmit(commonParams,ulParams);
            error.ul=err.UL;
            err=transmit(commonParams,llParams);
            error.ll=err.LL;

            % Accumulate total time
            totalTime(i)=totalTime(i)+c.time;

            % Adjust power ratio based on the current error rates
            commonParams=setPowerRatio(error,commonParams,ulParams,llParams);

        end
    end
    k % Display counter value
    i % Display current iteration
end

% Average total time over iterations
totalTime=totalTime./j;

toc % Stop stopwatch timer and display elapsed time
