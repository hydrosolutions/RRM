function out = onlineAssimilationTimer(in)
% USAGE: out = onlineAssimilationTimer(in)
%
% DESCRIPTION: This is the top-level function for the online assimilation model that
% starts timers and performs the consecutive model steps.
%
% VERSION:  07/01/2013, Version 1.0 (created)
%
% LAST MODIFIED: 07/01/2013
%
% Copyright (C) 2013 Tobias Siegfried
%
% This file is part of the iMoMo climate and flow forecasting and data 
% assimilation toolbox.

if isempty(timerfind('Tag','onlineAssimilationModel'))
    % start the assimilation time
    tIMOMO=timer('TimerFcn','onlineAssimilationModel','Tag','startOnlineAssimilation');
    set(tIMOMO,'ExecutionMode','fixedRate')
    set(tIMOMO,'Period',120) % just try to get it running every minute
else
    % do nothing (so far)
end


