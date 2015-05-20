function uTest_DateStr2Num(doSpeed)
% Automatic test: DateStr2Num
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_DateStr2Num(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-g V:006 Sum:LpKrpLBPoJzf Date:25-Mar-2011 14:33:32 $
% $License: BSD $
% $File: Published\DateStr2Num\uTest_DateStr2Num.m $

% Initialize: ==================================================================
if nargin == 0
   doSpeed = true;
end

NRandTest = 10000;

isMatlab7 = sscanf(version, '%f', 1) >= 7.0;

% Do the work: =================================================================
% Hello:
disp(['==== Test DateStr2Num:  ', datestr(now, 0)]);
disp(['  File: ', which('DateStr2Num')]);

fprintf('\n== Check current date:\n');

% Number format:
NowNum  = now;
NowTime = datenum(floor(datevec(NowNum)));  % Integer seconds only!
NowDate = floor(NowNum);

N0   = DateStr2Num(datestr(NowTime, 0), 0);
N1   = DateStr2Num(datestr(NowTime, 1), 1);
N29  = DateStr2Num(datestr(NowTime, 29), 29);
N30  = DateStr2Num(datestr(NowTime, 30), 30);
N31  = DateStr2Num(datestr(NowTime, 31), 31);
N300 = DateStr2Num([datestr(NowTime, 30), '.123'], 300);

if isequal(NowTime, N0) == 0
   error(['*** ', mfilename, ': DATESTR(0) format failed.']);
end
if isequal(NowDate, N1) == 0
   error(['*** ', mfilename, ': DATESTR(1) format failed.']);
end
if isequal(NowDate, N29) == 0
   error(['*** ', mfilename, ': DATESTR(29) format failed.']);
end
if isequal(NowTime, N30) == 0
   error(['*** ', mfilename, ': DATESTR(30) format failed.']);
end
if isequal(NowTime, N31) == 0
   error(['*** ', mfilename, ': DATESTR(31) format failed.']);
end
if round(((NowTime - N300) * 86400 + 0.123) * 1000) ~= 0
   error(['*** ', mfilename, ': yyyymmddTHHMMSS.FFF format failed.']);
end

disp('  ok: Successful for 0, 1, 29, 30, 31, 300 format');


disp('== Random test data...');
drawnow;
maxDate = datenum('31-Dec-9999');
minDate = datenum('01-Jan-0100');
NumPool = minDate + rand(NRandTest, 1) * (maxDate - minDate);
NumPool = datenummx(floor(datevec(NumPool)));
for iFormat = [0, 1, 29, 30, 31]
   if iFormat == 1 || iFormat == 29  % Without time:
      Want = floor(NumPool);
   else                              % With time
      Want = NumPool;
   end
   StrPool = cellstr(datestr(Want, iFormat));
   ConvC   = DateStr2Num(StrPool, iFormat);
   ConvS   = zeros(size(Want));
   for iP = 1:NRandTest
      ConvS(iP) = DateStr2Num(StrPool{iP}, iFormat);
   end
   
   if any(abs(Want - ConvS) > 131072 * eps)
      error(['*** ', mfilename, ': Failed to convert: "', ...
         StrPool{iP}, '" in format ' sprintf('%d', iFormat)]);
   end
   if ~isequal(ConvS, ConvC)
      error(['*** ', mfilename, ...
         ': Conversion differs for string an cell strings: "', char(10), ...
         StrPool{iP}, '" in format ' sprintf('%d', iFormat)]);
   end
end
fprintf('  ok: converted %d dates to formats 0, 1, 29, 30, 31\n', NRandTest);

% Speed: -----------------------------------------------------------------------
disp([char(10), '== Speed tests...']);
formatList = { ...
   0,   'dd-mmm-yyyy HH:MM:SS'; ...
   1,   'dd-mmm-yyyy'; ...
   29,  'yyyy-mm-dd'; ...
   30,  'yyyymmddTHHMMSS'; ...
   31,  'yyyy-mm-dd HH:MM:SS'; ...
   300, 'yyyymmddTHHMMSS.FFF'};

% Find a suiting number of loops:
C = datestr(floor(datevec(now)), 0);
if doSpeed
   iLoop     = 0;
   startTime = cputime;
   while cputime - startTime < 1.0
      v = datenum(C);     %#ok<NASGU>
      clear('v');
      iLoop = iLoop + 1;
   end
   nLoops = 100 * ceil(iLoop / ((cputime - startTime) * 50));
else
   disp('  Single loop => displayed times are not meaningful!');
   nLoops = 1;
end

fprintf('  Single string, %d loops on this machine.\n', nLoops);
fprintf([blanks(25), 'DATENUM   DateStr2Num\n']);
if isMatlab7
   knownFormats = [0, 1, 29, 30, 31, 300];
else  % Matlab 6:
   knownFormats = [0, 1];
end

for iFormat = knownFormats
   thisFormat = formatList{[formatList{:, 1}] == iFormat, 2};
   
   if iFormat < 32
      C = datestr(floor(datevec(now)), iFormat);
   elseif isMatlab7
      C = datestr(now, thisFormat);
   else  % Matlab 6:
      tmpNow = now;
      tmp = round(mod(tmpNow * 86400, 1) * 1000);
      C = [datestr(tmpNow, 30), sprintf('.%d', tmp)];
   end
   
   if isMatlab7
      tic;
      for i = 1:nLoops
         v = datenum(C, thisFormat);  %#ok<NASGU>
         clear('v');
      end
      tDatenum = toc + eps;
   else
      tic;
      for i = 1:nLoops
         v = datenum(C);  %#ok<NASGU>, no Format in Matlab 6
         clear('v');
      end
      tDatenum = toc + eps;
   end
   
   tic;
   for i = 1:nLoops * 10
      v = DateStr2Num(C, iFormat);  %#ok<NASGU>
      clear('v');
   end
   tNum = toc / 10;
   
   disp([sprintf('  %-20s', thisFormat), ...
      sprintf('   %-7.3f', tDatenum), '   ', sprintf('%-7.3g', tNum), ...
      '     => ', sprintf('%.2f', 100 * tNum / tDatenum), '%']);
end

% Cell string:
fprintf('\n');
if doSpeed
   nDate = 10000;
else
   nDate = 1000;
end

% Find a suiting number of loops:
if doSpeed
   nowValue = now;
   C = cell(nDate, 1);
   for iC = 1:nDate
      C{iC, 1} = datestr(nowValue, 0);
   end

   iLoop      = 0;
   startTime  = cputime;
   thisFormat = formatList{[formatList{:, 1}] == 0, 2};
   while cputime - startTime < 1.0
      v = datenum(C, thisFormat);     %#ok<NASGU>
      clear('v');
      iLoop = iLoop + 1;
   end
   nLoops = max(2, ceil(iLoop / (cputime - startTime)));
else
   disp('  Single loop => displayed times are not meaningful!');
   nLoops = 2;
end

fprintf('  {1 x %d} cell string, %d loops on this machine.\n', nDate, nLoops);
fprintf([blanks(25), 'DATENUM   DateStr2Num\n']);
if isMatlab7
   knownFormats = [0, 1, 29, 30, 31, 300];
else
   knownFormats = 0;
end

C = cell(nDate, 1);
for iFormat = knownFormats
   nowValue = now;
   for iC = 1:nDate
      if iFormat < 32
         C{iC} = datestr(nowValue, iFormat);
      elseif isMatlab7
         C{iC} = datestr(nowValue, 'yyyymmddTHHMMSS.FFF');
      else
         tmp = round(mod(nowValue * 86400, 1) * 1000);
         C{iC} = [datestr(nowValue, 30), sprintf('.%d', tmp)];
      end
   end
   
   formatIndex = [formatList{:, 1}] == iFormat;
   
   if isMatlab7
      thisFormat = formatList{formatIndex, 2};
      tic;
      for i = 1:nLoops
         v = datenum(C, thisFormat);  %#ok<NASGU>
         clear('v');
      end
      tDatenum = toc + eps;
   else
      tic;
      for i = 1:nLoops
         v = datenum(C);  %#ok<NASGU>, no Format in Matlab 6
         clear('v');
      end
      tDatenum = toc + eps;
   end
   
   tic;
   for i = 1:nLoops * 10   % Increase precision
      v = DateStr2Num(C, iFormat);  %#ok<NASGU>
      clear('v');
   end
   tNum = toc / 10;
   
   disp([sprintf('  %-20s', thisFormat), ...
      '   ', sprintf('%-7.3f', tDatenum), ...
      '   ', sprintf('%-7.3g', tNum), ...
      '     => ', sprintf('%.2f', 100 * tNum / tDatenum), '%']);
end

% Goodbye:
disp([char(10), 'DateStr2Num passed the tests.']);

return;
