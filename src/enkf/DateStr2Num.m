function D = DateStr2Num(S, F)
% DATESTR2NUM - Fast conversion of DATESTR to DATENUM
% The builtin DATENUM command is very powerful, but if you know the input
% format exactly, a specific MEX can be much faster: e.g. for single strings
% DateStr2Num takes 0.96% of the processing time of DATENUM (with specified
% date format), for a {1 x 10000} cell string, DateStr2Num needs 1.8% to 4.2%
% of the DATENUM time. (Matlab 2009a, MSVC 2008 using SSE2 optimization).
%
% D = DateStr2Num(S, F)
% INPUT:
%   S: String or cell string in DATESTR(F) format.
%      In opposite to DATENUM the validity of the input string is not checked
%      (e.g. 1 <= month <= 12).
%   F: Integer number defining the input format. Accepted:
%          0:  'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17
%          1:  'dd-mmm-yyyy'            01-Mar-2000
%         29:  'yyyy-mm-dd'             2000-03-01
%         30:  'yyyymmddTHHMMSS'        20000301T154517
%         31:  'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17
%      Including the milliseconds (not a valid DATEFORM number in DATESTR!):
%        300:  'yyyymmddTHHMMSS.FFF'    20000301T154517.123
%      Optional, default: 0.
%
% OUTPUT:
%   D: Serial date number. If S is a cell, D has is same size.
%
% EXAMPLES:
%   C = {'2010-06-29 21:59:13', '2010-06-29 21:59:13'};
%   D = DateStr2Num(C, 31)
%   >> [734318.916122685, 734318.916122685]
%   Equivalent Matlab command (but a column vector is replied ever):
%   D = datenum(C, 'yyyy-mm-dd HH:MM:SS')
%
% NOTES: The parsing of the strings works for clean ASCII characters only:
%   '0' must have the key code 48!
%   Month names must be English with the 2nd and 3rd charatcer in lower case.
%
% COMPILATION:
%   Windows: mex -O DateStr2Num.c
%   Linux:   mex -O CFLAGS="\$CFLAGS -std=c99" DateStr2Num.c
%   Precompiled Mex: http://www.n-simon.de/mex
%
% Tested: Matlab 6.5, 7.7, 7.8, 32bit, WinXP
% Compiler: LCC 2.4/3.8, BCC 5.5, Open Watcom 1.8, MSVC 2008
%         Compatibility to MacOS, Linux, 64 bit is assumed, but not tested.
% Author: Jan Simon, Heidelberg, (C) 2010-2011 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also DATESTR, DATENUM, DATEVEC.
% FEX: DateConvert 25594 (Jan Simon)

% $JRev: R-f V:006 Sum:b0uK5gDiBBdO Date:25-Mar-2011 14:33:32 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_DateStr2Num $
% $File: Published\DateStr2Num\DateStr2Num.m $
% History:
% 001: 29-Jun-2010 22:07, First version.
% 002: 30-Jun-2010 16:08, Accept formats 0, 1, 29, 30, 31.

% ==============================================================================
% The MEX is not found?! Perhaps the parent folder of this function is the
% current folder...
warning('JSimon:DateStr2Num:NoMex', 'Using slow M-version of DateStr2Num');

% Compared to DATENUM, this M-version is 20 times faster than DATENUM(S, 31) for
% a string and 1.1 to 1.4 times slower for a {1 x 10000} cell string.
if ischar(S)
   S = {S};
elseif ~iscellstr(S)
   error('JSimon:DateStr2Num:BadInputType', ...
      'Input must be a string or cell string.');
end

if nargin == 1
   F = 0;
end

nS = numel(S);
D  = zeros(size(S));
switch F
   case 0
      for iS = 1:nS
         aS    = S{iS};
         x     = aS - '0';
         month = (strfind('janfebmaraprmayjunjulaugsepoctnovdec', ...
            lower(aS(4:6))) + 2) / 3;
         year  = x(8) * 1000 + x(9) * 100 + x(10) * 10 + x(11);
         
         D(iS) = datenummx(year, month, x(1) * 10 + x(2), ...
            x(13) * 10 + x(14), x(16) * 10 + x(17), ...
            x(19) * 10 + x(20));
      end
      
   case 1
      for iS = 1:nS
         aS    = S{iS};
         x     = aS - '0';
         month = (strfind('janfebmaraprmayjunjulaugsepoctnovdec', ...
            lower(aS(4:6))) + 2) / 3;
         year  = x(8) * 1000 + x(9) * 100 + x(10) * 10 + x(11);
         
         D(iS) = datenummx(year, month, x(1) * 10 + x(2));
      end
      
   case 29
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(6)*10 + x(7), ...  % Month
            x(9)*10 + x(10));     % Day
      end
      
   case 30
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(5)*10 + x(6), ...  % Month
            x(7)*10 + x(8), ...  % Day
            x(10)*10 + x(11), x(12)*10 + x(13), x(14)*10 + x(15));
      end
      
   case 31
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(6)*10 + x(7), ...  % Month
            x(9)*10 + x(10), ... % Day
            x(12)*10 + x(13), x(15)*10 + x(16), x(18)*10 + x(19));
      end

   case 300  % yyyymmddTHHMMSS.FFF
      for iS = 1:nS
         x     = S{iS} - '0';
         D(iS) = datenummx( ...
            x(1)*1000 + x(2)*100 + x(3)*10 + x(4), ... % Year
            x(5)*10 + x(6), ...  % Month
            x(7)*10 + x(8), ...  % Day
            x(10)*10 + x(11), x(12)*10 + x(13), x(14)*10 + x(15) + ...
            x(17) / 10 + x(18) / 100 + x(19) / 1000);
      end
      
   otherwise
      error('JSimon:DateStr2Num:BadFormat', 'Bad format');
end

return;
