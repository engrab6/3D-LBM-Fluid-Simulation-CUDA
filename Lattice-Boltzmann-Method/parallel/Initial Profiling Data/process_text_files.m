% preamble; 
a = dir('output_*.txt')

filenames = {a(:).name}'

for (i = 1:length(filenames))
   currentFile = filenames{i};
   N(i) = str2num(strtok(currentFile, 'output_N=')); 
   d = importdata(currentFile);
   numPoints(i) = mean(d(:,1)); 
   runtime(i) = min(d(:,2)); % s
   speed(i) = max(d(:,3)); % Mlups
   bandwidth(i) = max(d(:,4)); % GiB/s 
end
% Data from nvprof for run_lbm; run_number percUtilization totalRuntime numberOfCalls maxRuntime avgRuntime minRuntime
run_lbm = [1 70.80  0.19360 1200  0.16134 0.15751  0.16560
2 64.46  0.24262 1452  0.16710 0.16106  0.18144
3 63.02  0.42119 1728  0.24375 0.22429  0.28448
4 64.22  0.76611 2028  0.37777 0.36388  0.42871
5 63.99  1.19449 2352  0.50786 0.49044  0.58129
6 64.24  1.88675 2700  0.69880 0.67181  0.79182
7 65.02  2.85056 3072  0.92792 0.89898  1.0311
8 62.73  3.68819 3468  1.0635 1.0255 1.1929
9 58.23  4.15292 3888  1.0681 1.0353 1.2166
10 58.76  5.74177 4332  1.3254 1.2972 1.4931
11 61.00  8.03512 4800  1.6740 1.6348 1.8557
12 62.16  10.8633 5292  2.0528 1.9948 2.2173
13 62.19  13.9420 5808  2.4005 2.2989 2.5817
14 54.33  12.0335 6348  1.8956 1.8576 2.1759
15 56.41  16.4731 6912  2.3833 2.3259 2.5960
16 57.43  21.2734 7500  2.8364 2.7765 3.0319
17 59.35  27.4881 8112  3.3886 3.3292 3.5825
18 60.43  34.6293 8748  3.9585 3.8919 4.1516
19 60.34  42.0560 9408  4.4702 4.3778 4.6891
20 61.07  52.1786 10092  5.1703 5.0563 5.3980
21 61.39  63.3723 10800  5.8678 5.7919 6.0652];
% Data from nvprof for set_f_temp; N percUtilization totalRuntime numberOfCalls maxRuntime avgRuntime minRuntime

set_f_temp = [   10  28.89 0.079012  1200 0.065843 0.064064 0.075072
                 11  35.24 0.13264  1452 0.091351 0.082625  0.11357
                 12  36.77 0.24572  1728  0.14220  0.13418  0.16685
                 13   35.63 0.42496  2028  0.20955  0.19914  0.24176
                 14  35.89 0.66998  2352  0.28485  0.27370  0.33044
                 15  35.67 1.04774  2700  0.38805  0.37069  0.43815
                 16  34.91 1.53058  3072  0.49823  0.47543  0.55975
                 17  37.20 2.18707  3468  0.63064  0.60196  0.70337
                 18  41.70 2.97435  3888  0.76501  0.73533  0.85150
                 19  41.18 4.02411  4332  0.92893  0.90436  1.0298
                 20  38.94 5.12934  4800  1.0686  1.0181  1.1997
                 21  37.80 6.60584  5292  1.2483  1.2011  1.3892
                 22  37.77 8.46856  5808  1.4581  1.4123  1.6081
                 23  45.62 10.1050  6348  1.5918  1.5437  1.7481
                 24  43.55 12.7191  6912  1.8401  1.7721  1.9939
                 25  42.53 15.7527  7500  2.1004  2.0466  2.2592
                 26  40.62 18.8126  8112  2.3191  2.2763  2.4805
                 27  39.54 22.6624  8748  2.5906  2.5346  2.7488
                 28  39.63 27.6228  9408  2.9361  2.8776  3.1031
                 29  38.91 33.2416 10092  3.2939  3.2307  3.4572
                 30  38.58 39.8285 10800  3.6878  3.5961  3.8662];
run_lbm(:,1) = run_lbm(:,1) + 9; 

dfo = defaultColorOrder();     

setDefaultFont('Ubuntu', 18); 
fullfig(); 
subplot(2,2,1)
setAxisDefaults(gca); hold on; 
plot(N.^3, runtime,'-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'white', 'MarkerSize', 10);
plot(N.^3, [run_lbm(:,3)./12 set_f_temp(:,3)./12], '-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'white', 'MarkerSize', 10);
xlabel('Number of Points Per Side'); 
ylabel('Runtime per Simulation (s)'); 
legend({'total', 'run\_lbm', 'set\_f\_temp'}, 'Color', 'none');
% saveCurrentFigure('runtime_vs_N.png'); 
subplot(2,2,2)
% fullfig(); 
setAxisDefaults(gca); hold on; 
plot(N.^3, speed, '-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'white', 'MarkerSize', 10); 
xlabel('Number of Points Per Side'); 
ylabel('Simulation Speed (Mlups)'); 
% saveCurrentFigure('speed_vs_N.png'); 
a = subplot(2,2,3);
% fullfig(); 
setAxisDefaults(gca); hold on; 
plot(N.^3, [run_lbm(:,6) set_f_temp(:,6)], '-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'white', 'MarkerSize', 10); 
xlabel('Number of Points Per Side'); 
ylabel('Average Kernel Runtime (ms)'); 
legend({'run\_lbm', 'set\_f\_temp'}, 'Color', 'none');
% saveCurrentFigure('avgTime_vs_N.png'); 
subplot(2,2,4)
setAxisDefaults(gca); hold on; 
plot(N.^3, bandwidth, '-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'white', 'MarkerSize', 10); 
xlabel('Number of Points Per Side'); 
ylabel('Bandwidth (GiB/s)'); 