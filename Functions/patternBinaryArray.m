function [patternsActive, patternNames] = patternBinaryArray(patterns, dataLength, patternNames)
% Turns the patterns structure into a TxP binary array, where there are P
% unique types of patterns in PATTERNS.TYPE and T is the maximum time step
% in PATTERNS.ENDTIME. PATTERNSACTIVE contains a 1 in [m, n] if pattern n
% is active at time m. PATTERNNAMES is a 1xP cell array containing the
% names of the P types of patterns in the order that they are indexed.

% The 3 columns are:
% 1) Current time step, 2) Row coordinate, 3) Column coordinate

if nargin < 3
    patternNames = sort(unique(patterns.type))';
end

p = length(patternNames);
if nargin == 1
    t = max(patterns.endTime);
else
    t = dataLength;
end
patternsActive = zeros(t, p);


% Loop over all pattern types
for ip = 1:p
    ipatterns = cropStructure(patterns, patternNames{ip});
    
    % Loop over every separate occurance of the pattern
    for ij = 1:length(ipatterns.type)
        patternsActive(ipatterns.startTime(ij):ipatterns.endTime(ij), ip) = 1;
    end  
end

patternsActive = patternsActive == 1;