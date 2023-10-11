% Function definition for 'transmit', which takes two arguments: 'common' and 'receiver'.
% 'common' appears to contain parameters common to the transmission, while 'receiver' contains information about the receiver's status.
function [error] = transmit(common, receiver)

    % Initialize a structure array 'parameters' to store various transmission settings.
    parameters = struct([]);

    % Start setting up the transmission parameters, which include:
    % - numLayers: Likely the number of layers in a MIMO (Multiple Input, Multiple Output) setup.
    % - perfectEstimization: A flag indicating whether perfect channel estimation is assumed.
    % - cr1, cr2, cr3: Possibly 'coding rates' for three different FEC (Forward Error Correction) schemes or layers.
    % - slots: The number of time slots available for transmission.
    % - nTxAnts: The number of transmitting antennas (relevant for MIMO setups).
    % - channel: The channel over which the transmission is occurring.
    % - g2, g3: Could be gain values for different transmission paths or channels.
    % - modulationUL, modulationLL, modulationL3: Modulation schemes for different layers or types of transmission.
    % - nRxAnts: The number of receiving antennas (also relevant for MIMO).
    % - snr: The Signal-to-Noise Ratio at the receiver.
    % - delayProfile: A specific profile for channel delay, which differs based on whether the Line-of-Sight (LOS) is present.
    parameters(1).numLayers = 1;
    parameters(1).perfectEstimization = true;
    parameters(1).cr1 = 490 / 1024;
    parameters(1).cr2 = 490 / 1024;
    parameters(1).cr3 = 490 / 1024;
    parameters(1).slots = common.slots;
    parameters(1).nTxAnts = common.nTxAnts;
    parameters(1).channel = common.channel;
    parameters(1).g2 = common.g2;
    parameters(1).g3 = common.g3;
    parameters(1).modulationUL = common.modulationUL;
    parameters(1).modulationLL = common.modulationLL;
    parameters(1).modulationL3 = common.modulationL3;
    parameters(1).nRxAnts = common.nRxAnts;
    parameters(1).snr = receiver.rxSnr;

    % Choose the delay profile based on whether the transmission is Line-of-Sight (LOS)
    if receiver.position.los == 0
        parameters(1).delayProfile = "TDL-A"; % Non-LOS scenario
    else
        parameters(1).delayProfile = "TDL-D"; % LOS scenario
    end

    % Initialize a structure 'error' to store the results of the simulation.
    error = struct([]);

    % Call a function 'physicalLayer5GNThroughput' which likely simulates the transmission based on the 'parameters' and returns some measure of success or error rate (e.g., throughput, bit error rate).
    error = physicalLayer5GNThroughput(parameters);

    % Access a specific field 'UL' from the 'error' structure, possibly indicating 'Uplink' errors or throughput. However, this line doesn't seem to affect the output or any variable and might be redundant or part of a debugging process.
    error.UL;
end
