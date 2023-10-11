% Define the main function that calculates path loss
function [pl] = pathloss (receiver, common)
% Initialize path loss to 0
pl = 0;

% Retrieve constants for the scenario
c = myPackageConstant();

% Determine the scenario and calculate path loss based on the scenario
switch (c.scenario)
    case ("RMa") % Case for Rural Macro-cell scenario
        h = 5; % Average building height in meters
        w = 20; % Average street width in meters

        % Calculate breakpoint distance based on the frequency and heights of transmitter and receiver
        dBP = 2 * pi * c.txHeight * c.rxHeight * (c.freq * 1e9 / 3e8);

        % Calculate path loss based on different models
        pl1 = pathL(receiver.position.distance3d, h); % Path loss based on 3D distance and building height
        pl2 = pathL(dBP, h) + 40 * log10(receiver.position.distance3d / dBP); % Path loss considering breakpoint distance

        % Non-line-of-sight (NLOS) path loss calculation
        pln = 161.04 - 7.1 * log10(w) + 7.5 * log10(h) - (24.37 - 3.7 * (h / c.txHeight)^2) * log10(c.txHeight) + (43.42 - 3.1 * log10(c.txHeight)) * (log10(receiver.position.distance3d) - 3) + 20 * log10(c.freq) - (3.2 * (log10(11.75 * c.rxHeight))^2 - 4.97);

        % Select the appropriate path loss based on the distance
        if (receiver.position.distance < dBP)
            pl = pl1;
        else
            pl = pl2;
        end

        % If the Line-of-Sight (LOS) is not available, use the maximum path loss
        if (receiver.position.los == 0)
            pl = max(pl, pln);
        end

    case("UMa") % Case for Urban Macro-cell scenario
        hBS = c.txHeight - 1; % Height of Base Station minus 1
        hR = c.rxHeight - 1; % Height of Receiver minus 1

        % Calculate breakpoint distance
        dBP = 4 * hBS * hR * (c.freq * 1e9 / 3e8);

        % LOS path loss models
        pl1 = 28 + 22 * log10(receiver.position.distance3d) + 20 * log10(c.freq);
        pl2 = 28 + 40 * log10(receiver.position.distance3d) + 20 * log10(c.freq) - 9 * log10((dBP)^2 + (hBS - hR)^2);

        % NLOS path loss model
        pln = 13.54 + 39.08 * log10(receiver.position.distance3d) + 20 * log10(c.freq) - 0.6 * (hR - 1.5);

        % Select the appropriate path loss
        if (receiver.position.distance < dBP)
            pl = pl1;
        else
            pl = pl2;
        end

        % If no Line-of-Sight, choose the maximum path loss
        if (receiver.position.los == 0)
            pl = max(pl, pln);
        end

    case("UMi") % Case for Urban Micro-cell scenario
        hBS = c.txHeight - 1; % Height of Base Station minus 1
        hR = c.rxHeight - 1; % Height of Receiver minus 1

        % Calculate breakpoint distance
        dBP = 4 * hBS * hR * (c.freq * 1e9 / 3e8);

        % LOS path loss models
        pl1 = 32.4 + 21 * log10(receiver.position.distance3d) + 20 * log10(c.freq);
        pl2 = 32.4 + 40 * log10(receiver.position.distance3d) + 20 * log10(c.freq) - 9.5 * log10((dBP)^2 + (hBS - hR)^2);

        % NLOS path loss model
        pln = 35.3 * log10(receiver.position.distance3d) + 22.4 + 21.3 * log10(c.freq) - 0.3 * (hR - 1.5);

        % Select the appropriate path loss
        if (receiver.position.distance < dBP)
            pl = pl1;
        else
            pl = pl2;
        end

        % If no Line-of-Sight, choose the maximum path loss
        if (receiver.position.los == 0)
            pl = max(pl, pln);
        end
end
end

% Helper function to calculate path loss
function [pl] = pathL(distance, h)
c = myPackageConstant(); % Retrieve constants for the scenario

% Path loss calculation based on frequency, distance, and building height
pl = 20 * log10(40 * pi * distance * c.freq / 3) + min(0.03 * h^1.72, 10) * log10(distance) - min(0.044 * h^1.72, 14.77) + 0.002 * log10(h) * distance;
end
