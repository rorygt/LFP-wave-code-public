% Script to quickly visualize the results of the MAIN function

loadNew = 1;
phaseCmap = pmkmp_new(256, 'ostwald_o');
skipLength = 5;
startIndex = 3800;
spikeSmoothSpan = 21;
vScale = 200;

% Load files
if loadNew
    % Load spikes
    load(outputFile1Name, 'spikeMatrix', 'Fs')
    % Make sure that spikes are represented by ones
    spikeMatrix = spikeMatrix ~= 0;
    % Smooth spikes with a moving average filter to obtain spike rate
    spikeRate = zeros(size(spikeMatrix));
    for ii = 1:size(spikeMatrix, 1)
        spikeRate(ii,:) = smooth(spikeMatrix(ii,:), spikeSmoothSpan);
    end
    % Convert to grid and SI units
    spikeRate = vector2grid(spikeRate) * Fs;
    
    % Load phase maps, velocity fields and pattern information
    load(outputFile2Name, 'phase', 'experiment', 'file', 'discardTimeSteps')
    phase = vector2grid(phase);
    load(outputFile3Name, 'velocityX', 'velocityY', 'badChannelsGrid', ...
        'pwActive', 'syActive', 'patternCentreStruct')
    % Convert critical point structure to more usable form
    centres = combinePatternCentres(patternCentreStruct);
    centres = sortStruct(centres, 'time');
    % Trim spikes to correct time values to match other data
    spikeRate = spikeRate(:,:,discardTimeSteps:(size(spikeMatrix,2) - ...
        discardTimeSteps));
%     % Calculate binary array of spike wave activity
%     [phi, v0, v_direction] = orderParameter(spikevx, spikevy);
%     [~, spwActive] = addToPatternsStructure(...
%         'planeWave', phi>=planeWaveThreshold, [], params);
%     spwActive = spwActive(discardTimeSteps:(end - discardTimeSteps));
    
end


for index = startIndex:skipLength:size(phase, 3)
    
    % View phase
    subplot(2,2,1)
    displayGrid(phase(:,:,index), phaseCmap, [-pi pi], 0, 1);
    title(sprintf('%s-%d phase: %0.3f s, step %d', ...
        experiment, file, index/Fs, index))
    axis square
    
    % View velocity field with pattern types
    subplot(2,2,3)
    quiver(0.5:9.5, 0.5:9.5, vScale*flipud(velocityX(:,:,index)), ...
        vScale*flipud(-velocityY(:,:,index)), 0)
    set(gca, 'XTick', 0:2.5:10, 'YTick', 0:2.5:10)%, ...
       % 'XTickLabel', [], 'YTickLabel', []);
    axis([0 10 0 10])
    axis square
    title('Phase velocity field')
    
    % Show criticial point centres
    hold on
    critIndex = find(centres.time == index);
    for ii = 1:length(critIndex)
        icrit = critIndex(ii);
        % Choose marker colour/type based on pattern
        switch char(centres.name(icrit))
            case 'saddle'
                mcol = 'k';
                mtype = 'x';
            case 'sink'
                mcol = 'r';
                mtype = 'filled';
            case 'source'
                mcol = 'r';
                mtype = 'o';
            case 'spiralIn'
                mcol = 'b';
                mtype = 'filled';
            case 'spiralOut'
                mcol = 'b';
                mtype = 'o';
        end
        
        % Coordinates are in row, column form so convert to cartesian
        % coordinates for plotting
        scatter(centres.coords(icrit,2)-0.5, 10.5-centres.coords(icrit,1), 80, ...
            mcol, mtype, 'LineWidth', 2)
        
    end
    hold off
    
    % View multi-unit firing rate at each channel
    subplot(2,2,2)
    displayGrid(spikeRate(:,:,index), jet, [0, Fs], 0, 1);
    title(sprintf('Multi-unit firing rate, %d ms window', spikeSmoothSpan));
    
    % Add indicator for synchronized activity
    subplot(6,4,19)
    imagesc([0 1], [0 1], syActive(index), [0 1])
    colormap gray
    set(gca, 'XTick', [], 'YTick', [])
    title('Synchrony')
    
    % Add indicator for plane wave
    subplot(6,4,20)
    imagesc([0 1], [0 1], pwActive(index), [0 1])
    colormap gray
    set(gca, 'XTick', [], 'YTick', [])
    title('Plane wave')
    
%     % TEMP: Add indicator for spike wave
%     subplot(4,4,15)
%     imagesc([0 1], [0 1], spwActive(index), [0 1])
%     colormap gray
%     title('Spike plane wave')
    
    drawnow
    w = waitforbuttonpress;
    
end
