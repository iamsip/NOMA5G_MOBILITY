function [error]=physicalLayer5GNThroughput(parameters)
    %%
    % Simulation Parameters
    
    SNRdB = parameters.snr;                % SNR in dB
    totalNoSlots = parameters.slots;         % Number of slots to simulate
    perfectEstimation = parameters.perfectEstimization; % Perfect synchronization and channel estimation
    rng("default");            % Set default random number generator for repeatability
    %%
    % Carrier Configuration
    
    carrier = nrCarrierConfig;
    
    %%
    % PDSCH and DM-RS Configuration
    
    pdschUL = nrPDSCHConfig;
    pdschUL.Modulation = parameters.modulationUL;
    pdschUL.NumLayers = parameters.numLayers ;
    pdschUL.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation
    
    %%
    %Set the DM-RS parameters. To improve channel estimation, add an additional DM-RS position
    
    pdschUL.DMRS.DMRSAdditionalPosition = 1;
    pdschUL.DMRS.DMRSConfigurationType = 1;
    pdschUL.DMRS.DMRSLength = 2;
    
    pdschLL=pdschUL;
    pdschL3=pdschUL;
    pdschLL.Modulation = parameters.modulationLL;
    pdschL3.Modulation = parameters.modulationL3;
    %%
    % DL-SCH Configuration
    
    NHARQProcesses = 16;     % Number of parallel HARQ processes
    rvSeq = [0 2 3 1];
    
    % Coding rate
    if pdschLL.NumCodewords == 1
        codeRate = parameters.cr1;
        codeRatell = parameters.cr2;
        codeRatel3 = parameters.cr3;
    else
        codeRate = [490 490]./1024;
    end
    

    % Create DL-SCH encoder object
    encodeDLSCH = nrDLSCH;
    encodeDLSCH.MultipleHARQProcesses = true;
    encodeDLSCH.TargetCodeRate = codeRate;

    encodeDLSCHll=encodeDLSCH;
    encodeDLSCHll.TargetCodeRate=codeRatell;

    encodeDLSCHl3=encodeDLSCH;
    encodeDLSCHl3.TargetCodeRate=codeRatel3;
    
    % Create DLSCH decoder object
    decodeDLSCH = nrDLSCHDecoder;
    decodeDLSCH.MultipleHARQProcesses = true;
    decodeDLSCH.TargetCodeRate = codeRate;
    decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
    decodeDLSCH.MaximumLDPCIterationCount = 6;

    decodeDLSCHll=decodeDLSCH;
    decodeDLSCHll.TargetCodeRate=codeRatell;
    decodeDLSCHl3=decodeDLSCH;
    decodeDLSCHl3.TargetCodeRate=codeRatel3;
    

    %%
    % HARQ Management
    
    harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdschUL.NumCodewords);
    
    %%
    % Channel Configuration
    % Check that the number of layers is valid for the number of antennas
    
    nTxAnts = parameters.nTxAnts;
    nRxAnts = parameters.nRxAnts;
    if pdschUL.NumLayers > min(nTxAnts,nRxAnts)
        errorTransmit("The number of layers ("+string(pdschUL.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")")
    end
    
    %%
    % Create a channel object
    channel = parameters.channel;
    channel.DelayProfile = parameters.delayProfile;
    channel.NumTransmitAntennas = nTxAnts;
    channel.NumReceiveAntennas = nRxAnts;
    
    %%
    % Set the channel sample rate to that of the OFDM signal. To obtain the sampling rate of the OFDM signal, use the nrOFDMInfo function.
    
    ofdmInfo = nrOFDMInfo(carrier);
    channel.SampleRate = ofdmInfo.SampleRate;
    
    %%
    % Transmission and Reception
    
    constPlot = comm.ConstellationDiagram;                                          % Constellation diagram object
    constPlot.ReferenceConstellation = getConstellationRefPoints(pdschUL.Modulation); % Reference constellation values
    constPlot.EnableMeasurements = 1;                                               % Enable EVM measurements
    
    % Initial timing offset
    offset = 0;
    
    estChannelGrid = getInitialChannelEstimate(channel,carrier);
    newPrecodingWeight = getPrecodingMatrix(pdschLL.PRBSet,pdschLL.NumLayers,estChannelGrid);
    throughputUL=0;
    throughputLL=0;
    throughputL3=0;
    errorRateUL=0;
    errorRateLL=0;
    errorRateL3=0;
    for nSlot = 1:totalNoSlots
        % New slot
        %   carrier.NSlot = nSlot;
        %%
        % Calculate Transport Block Size
        % Generate PDSCH indices info, which is needed to calculate the transport block size
        [pdschIndicesLL,pdschInfoLL] = nrPDSCHIndices(carrier,pdschLL);
        
        % Calculate transport block sizes
        Xoh_PDSCH = 0;
        trBlkSizesLL = nrTBS(pdschLL.Modulation,pdschUL.NumLayers,numel(pdschUL.PRBSet),pdschInfoLL.NREPerPRB,codeRatell,Xoh_PDSCH);
       
        %%
        %enhance layer data
        for cwIdx = 1:pdschUL.NumCodewords
            if harqEntity.NewData(cwIdx)
                % Create and store a new transport block for transmission
                trBlkLL = randi([0 1],trBlkSizesLL(cwIdx),1);
                setTransportBlock(encodeDLSCHll,trBlkLL,cwIdx-1,harqEntity.HARQProcessID);
        
                % If the previous RV sequence ends without successful decoding, flush the soft buffer
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCHll,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end
        %%
        %DL-SCH Encoding
        
        codedTrBlockLL = encodeDLSCHll(pdschLL.Modulation,pdschUL.NumLayers,pdschInfoLL.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        %%
        %PDSCH Modulation and MIMO Precoding
        
        pdschSymbolsLL = nrPDSCH(carrier,pdschLL,codedTrBlockLL);
                %%
        %L3 data 

        [pdschIndicesL3,pdschInfoL3] = nrPDSCHIndices(carrier,pdschL3);
        trBlkSizesL3 = nrTBS(pdschL3.Modulation,pdschL3.NumLayers,numel(pdschL3.PRBSet),pdschInfoL3.NREPerPRB,codeRatel3,Xoh_PDSCH);
         for cwIdx = 1:pdschUL.NumCodewords
            if harqEntity.NewData(cwIdx)
                % Create and store a new transport block for transmission
                trBlkL3 = randi([0 1],trBlkSizesL3(cwIdx),1);
                setTransportBlock(encodeDLSCHl3,trBlkL3,cwIdx-1,harqEntity.HARQProcessID);
        
                % If the previous RV sequence ends without successful decoding, flush the soft buffer
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCHll,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
         end

        codedTrBlockL3 = encodeDLSCHl3(pdschL3.Modulation,pdschUL.NumLayers,pdschInfoL3.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        pdschSymbolsL3 = nrPDSCH(carrier,pdschL3,codedTrBlockL3);
          %%
        % HARQ Processing (Buffer Management)
        % Get new transport blocks and flush decoder soft buffer, as required
        [pdschIndicesUL,pdschInfoUL] = nrPDSCHIndices(carrier,pdschUL);
        trBlkSizesUL = nrTBS(pdschUL.Modulation,pdschUL.NumLayers,numel(pdschUL.PRBSet),pdschInfoUL.NREPerPRB,codeRate,Xoh_PDSCH);


        for cwIdx = 1:pdschUL.NumCodewords
            if harqEntity.NewData(cwIdx)
                % Create and store a new transport block for transmission
                trBlkUL = randi([0 1],trBlkSizesUL(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlkUL,cwIdx-1,harqEntity.HARQProcessID);
        
                % If the previous RV sequence ends without successful decoding, flush the soft buffer
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end
        
        %%
        %DL-SCH Encoding
        
        codedTrBlockUL = encodeDLSCH(pdschUL.Modulation,pdschUL.NumLayers,pdschInfoUL.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        
        %%
        %PDSCH Modulation and MIMO Precoding
        
        pdschSymbolsUL = nrPDSCH(carrier,pdschUL,codedTrBlockUL);
        
        %%
        % LDM integration
        G3=parameters.g3; %LDM power ratio dB
        g3=10^(G3/10);
        pdschSymbolsLLC=LDMSuperposition(pdschSymbolsLL,pdschSymbolsL3,g3);
        G2=parameters.g2; %LDM power ratio dB
        g2=10^(G2/10);
        pdschSymbolsL=LDMSuperposition(pdschSymbolsUL,pdschSymbolsLLC,g2);

    %%
        precodingWeights = newPrecodingWeight;
        pdschSymbolsPrecoded = pdschSymbolsL*precodingWeights;
        
        %%
        %PDSCH DM-RS Generation
        
        dmrsSymbols = nrPDSCHDMRS(carrier,pdschUL);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdschUL);
        
        %%
        %Mapping to Resource Grid
        
        pdschGrid = nrResourceGrid(carrier,nTxAnts);
        
        [~,pdschAntIndices] = nrExtractResources(pdschIndicesUL,pdschGrid);
        pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;
        
        % PDSCH DM-RS precoding and mapping
        for p = 1:size(dmrsSymbols,2)
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
            pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*precodingWeights(p,:);
        end
        
        %%
        %OFDM Modulation
        
        [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
        
        %%
        %Propagation Channel
        %Pad the input signal with enough zeros to ensure that the generated signal is flushed out of the channel filter
    
        chInfo = info(channel);
        maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
        txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];
        
        %Send the signal through the channel and add noise
        
        [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
        noise = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(rxWaveform));
        rxWaveform = rxWaveform + noise;
        
        %%
        %Timing Synchronization
        
        if perfectEstimation
            % Get path filters for perfect timing estimation
            pathFilters = getPathFilters(channel);
            [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
        else
            [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
            offset = hSkipWeakTimingOffset(offset,t,mag);
        end
        rxWaveform = rxWaveform(1+offset:end,:);
        
        %%
        %OFDM Demodulation
        
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
        
        %%
        %Channel Estimation
        
        if perfectEstimation
            % Perform perfect channel estimation between transmit and receive
            % antennas.
            estChGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
        
            % Get perfect noise estimate (from noise realization)
            noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));
        
            % Get precoding matrix for next slot
            newPrecodingWeight = getPrecodingMatrix(pdschUL.PRBSet,pdschUL.NumLayers,estChGridAnts);
        
            % Apply precoding to estChGridAnts. The resulting estimate is for
            % the channel estimate between layers and receive antennas.
            estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
        else
            % Perform practical channel estimation between layers and receive
            % antennas.
            [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdschUL.DMRS.CDMLengths);
        
            % Remove precoding from estChannelGrid before precoding
            % matrix calculation
            estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
        
            % Get precoding matrix for next slot
            newPrecodingWeight = getPrecodingMatrix(pdschUL.PRBSet,pdschUL.NumLayers,estChGridAnts);
        end
        
    %     figure (01)
    %     mesh(abs(estChGridLayers(:,:,1,1)));
    %     title('Channel Estimate');
    %     xlabel('OFDM Symbol');
    %     ylabel("Subcarrier");
    %     zlabel("Magnitude");
        
        
        %%
        % Equalization
        
        [pdschRx,pdschHest] = nrExtractResources(pdschIndicesLL,rxGrid,estChGridLayers);
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
        
    %     figure(02)
    %     constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
    %     constPlot.ShowLegend = true;
    %     % Constellation for the first layer has a higher SNR than that for the
    %     % last layer. Flip the layers so that the constellations do not mask
    %     % each other.
    %     constPlot(fliplr(pdschEq));
        
        %%
        % PDSCH Decoding
        
        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdschUL,pdschEq,noiseEst);
        % Scale LLRs by CSI
        csi = nrLayerDemap(csi);                                    % CSI layer demapping
        for cwIdx = 1:pdschUL.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale
        end
        
        %%
        % DL-SCH Decoding
        
        decodeDLSCH.TransportBlockLength = trBlkSizesUL;
        [decbitsUL,blkerrUL] = decodeDLSCH(dlschLLRs,pdschUL.Modulation,pdschUL.NumLayers, ...
            harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
        %%
        % Core layer reconstruction
        setTransportBlock(encodeDLSCH,decbitsUL,cwIdx-1,harqEntity.HARQProcessID);
        rCodedTrBlockUL = encodeDLSCH(pdschUL.Modulation,pdschUL.NumLayers,pdschInfoUL.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        rpdschSymbolsUL = nrPDSCH(carrier,pdschUL,rCodedTrBlockUL);
        
        %%
        % LDM SIC
        rpdschSymbolsLL= LDMSubstraction(pdschEq,rpdschSymbolsUL,g2);
        
        %%
        % LDM DECODING
       
        
        % PDSCH Decoding
        
        [dlschLLRsE,rxSymbols] = nrPDSCHDecode(carrier,pdschLL,rpdschSymbolsLL,noiseEst);
        % Scale LLRs by CSI
        
        %%
        % DL-SCH Decoding
        
        decodeDLSCHll.TransportBlockLength = trBlkSizesLL;
        [decbitsLL,blkerrLL] = decodeDLSCHll(dlschLLRsE,pdschLL.Modulation,pdschUL.NumLayers, ...
            harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
        %%
        %l3 DECODING
        % LL reconstruction
        setTransportBlock(encodeDLSCHll,decbitsLL,cwIdx-1,harqEntity.HARQProcessID);
        rCodedTrBlockLL = encodeDLSCHll(pdschLL.Modulation,pdschUL.NumLayers,pdschInfoLL.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);;
        rCodedTrSymbolLL = nrPDSCH(carrier,pdschLL,rCodedTrBlockLL);

        % LDM SIC
        rpdschSymbolsL3= LDMSubstraction(rpdschSymbolsLL,rCodedTrSymbolLL,g3);
        [dlschLLRsE,rxSymbols] = nrPDSCHDecode(carrier,pdschL3,rpdschSymbolsL3,noiseEst);
        decodeDLSCHl3.TransportBlockLength = trBlkSizesL3;
        [decbitsL3,blkerrL3] = decodeDLSCHl3(dlschLLRsE,pdschL3.Modulation,pdschUL.NumLayers, ...
            harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        %%
        % error calculation
        errorUL=nnz(xor(decbitsUL, trBlkUL));
        throughputUL=throughputUL+(trBlkSizesUL-errorUL)*100/trBlkSizesUL;
        errorRateUL=errorRateUL+errorUL/trBlkSizesUL;

        errorLL=nnz(xor(decbitsLL, trBlkLL));
        throughputLL=throughputLL+(trBlkSizesLL-errorLL)*100/trBlkSizesLL;
        errorRateLL=errorRateLL+errorLL/trBlkSizesLL;

        errorL3=nnz(xor(decbitsL3, trBlkL3));
        throughputL3=throughputL3+(trBlkSizesL3-errorL3)*100/trBlkSizesL3;
        errorRateL3=errorRateL3+errorL3/trBlkSizesL3;
        %errorRateUL=errorRateUL+((trBlkSizesUL-errorUL)+(trBlkSizesLL-errorLL))/trBlkSizesLL;
        %%
        % HARQ Process Update
        
       % statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);
       % disp("Slot "+(nSlot));
        
    end % for nSlot = 1:totalNoSlots
    throughputUL=throughputUL./totalNoSlots;
    throughputLL=throughputLL./totalNoSlots;
    throughputL3=throughputL3./totalNoSlots;
    
    errorRateUL=errorRateUL./totalNoSlots;
    errorRateLL=errorRateLL./totalNoSlots;
    errorRateL3=errorRateL3./totalNoSlots;

    error=struct([]);
    error(1).TUL=throughputUL;
    error(1).TLL=throughputLL;
    error(1).TL3=throughputL3;
    error(1).UL=errorRateUL;
    error(1).LL=errorRateLL;
    error(1).L3=errorRateL3;
    release(channel)
end
%%
% Local Functions
function noise = generateAWGN(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
% Generate AWGN for a given value of SNR in dB (SNRDB), which is the
% receiver SNR per RE and antenna, assuming the channel does
% not affect the power of the signal. NRXANTS is the number of receive
% antennas. NFFT is the FFT size used in OFDM demodulation. SIZERXWAVEFORM
% is the size of the receive waveform used to calculate the size of the
% noise matrix.

    % Normalize noise power by the IFFT size used in OFDM modulation, as
    % the OFDM modulator applies this normalization to the transmitted
    % waveform. Also normalize by the number of receive antennas, as the
    % channel model applies this normalization to the received waveform by
    % default. The SNR is defined per RE for each receive antenna (TS
    % 38.101-4).
    SNR = 10^(SNRdB/10); % Calculate linear noise gain
    N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
    noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
end
    
function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% Calculate precoding matrix given an allocation and a channel estimate
    
    % Allocated subcarrier indices
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);
    
    % Average channel estimate
    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);
    
    % SVD decomposition
    [~,~,V] = svd(Hest);
    
    wtx = V(:,1:NLayers).';
    wtx = wtx/sqrt(NLayers); % Normalize by NLayers
end

function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Obtain an initial channel estimate for calculating the precoding matrix.
% This function assumes a perfect channel estimate

    % Clone of the channel
    chClone = channel.clone();
    chClone.release();

    % No filtering needed to get channel path gains
    chClone.ChannelFiltering = false;    
    
    % Get channel path gains
    [pathGains,sampleTimes] = chClone();
    
    % Perfect timing synchronization
    pathFilters = getPathFilters(chClone);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end

function refPoints = getConstellationRefPoints(mod)
% Calculate the reference constellation points for a given modulation
% scheme.
    switch mod
        case "QPSK"
            nPts = 4;
        case "16QAM"
            nPts = 16;
        case "64QAM"
            nPts = 64;
        case "256QAM"
            nPts = 256;            
    end
    binaryValues = int2bit(0:nPts-1,log2(nPts));
    refPoints = nrSymbolModulate(binaryValues(:),mod);
end

function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate.

    % Linearize 4-D matrix and reshape after multiplication
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid = estChannelGrid*W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end
