function r = wait_for_file(fn, wt)

r = 0;
counter = 0;
while ~exist(fn,'file') && counter<wt

pause(1);
counter = counter+1;
disp(counter);
end
if exist(fn,'file'), r = 1;end