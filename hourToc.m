function hourToc(tStart)
% FORMAT hourToc(tStart)
% A great way to see how long your scripts are taking when you're testing
% and streamlining them.
%
% tStart:   Result of tic (use tStart = tic to start timer).
% output:   Will print a string of the time elapsed since tStart was
%           created in hours, minutes, and seconds.

tEnd = toc(tStart);

secs = floor(rem(rem(tEnd, 3600), 60));

if tEnd > 60
    mins = floor(rem(tEnd, 3600) / 60);
end

if tEnd > 3600
    hrs = floor(tEnd / 3600);
    fprintf('%d hours, %d minutes, and %d seconds\n', hrs, mins, secs);
elseif tEnd > 60
    fprintf('%d minutes and %d seconds\n', mins, secs);
else
    fprintf('%d seconds\n', secs);
end
end