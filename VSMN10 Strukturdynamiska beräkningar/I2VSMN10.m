%%% Project 2 VSMN10 - Hugo Persson, Max Kristensen, Isabella McAuley Skriver

%% Half-power bandwidth

clear all; clc; close all;

%File data
filename = 'Transient_group_47.xlsx';
data = readmatrix(filename, 'Sheet', 'Sheet1');
time = data(:,1);
acc1 = data(:,2);
acc2 = data(:,3);
acc3 = data(:,4);
force = data(:,5);
accData = [acc1, acc2, acc3];

%FFT
[freqF, Force_fft] = fftoperator(time, force);
zeta = zeros(1, 3);

figure(1)

for i = 1:3
    acc = accData(:, i);
    [freqA, Acc_fft] = fftoperator(time, acc);
    FRF_mag = Acc_fft ./ Force_fft; %Magnitude
    [FRFmax, idx_max] = max(FRF_mag); %Peak
    f_res = freqA(idx_max); %Resonance frequency
    HP_mag = FRFmax / sqrt(2); %Half-power magnitude

    %Half-power bandwidth method
    idx_left = find(FRF_mag(1:idx_max) <= HP_mag, 1, 'last');
    idx_right = find(FRF_mag(idx_max:end) <= HP_mag, 1, 'first') + (idx_max - 1);

    %Calculating damping ratio
    f1 = freqA(idx_left);
    f2 = freqA(idx_right);
    zeta(i) = (f2 - f1) / (2 * f_res);

    %Plotting
    subplot(3,1,i);
    plot(freqA, FRF_mag, 'b-'); 
    hold on;
    plot(f_res, FRFmax, 'ro');
    plot([f1 f2], [HP_mag HP_mag],'g*');
    xlabel('Frequency (Hz)');
    title(sprintf('Accelerometer %d: \\zeta = %.5f', i, zeta(i)));
    axis([1000 2500 0 30])
    grid on;
end

%Results
for i = 1:3
    fprintf('Accelorometer %d : d = %.5f\n', i, zeta(i));
end

%% Compare accelerations

clear all; clc; close all;

%Experimental
filename = 'Transient_group_47.xlsx';
datae = readmatrix(filename, 'Sheet', 'Sheet1');
timee = datae(:,1);
acc1e = datae(:,2);
acc2e = datae(:,3);
acc3e = datae(:,4);
forcee = datae(:,5);
accDatae = [acc1e, acc2e, acc3e];

%Model 1 - Transient Loading
filename = 'Transientdata.xlsx';
data1 = readmatrix(filename, 'Sheet', 'Sheet1');
time1 = data1(:,1);
acc11 = data1(:,2);
acc21 = data1(:,3);
acc31 = data1(:,4);
accData1 = [acc11, acc21, acc31];

%Model 2 - 50 modes
filename = '50 moder.xlsx';
data2 = readmatrix(filename, 'Sheet', 'Sheet1');
time2 = data2(:,1);
acc12 = data2(:,2);
acc22 = data2(:,3);
acc32 = data2(:,4);
accData2 = [acc12, acc22, acc32];

%Model 3 - 100 modes
filename = '100 moder.xlsx';
data3 = readmatrix(filename, 'Sheet', 'Sheet1');
time3 = data3(:,1);
acc13 = data3(:,2);
acc23 = data3(:,3);
acc33 = data3(:,4);
accData3 = [acc13, acc23, acc33];

%Model 4 - 200 modes
filename = '200 moder.xlsx';
data4 = readmatrix(filename, 'Sheet', 'Sheet1');
time4 = data4(:,1);
acc14 = data4(:,2);
acc24 = data4(:,3);
acc34 = data4(:,4);
accData4 = [acc14, acc24, acc34];

%Low-pass filter
fs=10000;
cutoff=100;
order=1;
[b,a]=butter(order, cutoff/(fs/2), 'low');
acc_filtered=filter(b,a,accDatae);

%Root-Mean-Square
Experimental=mean(abs(accDatae(24:3278+24,:)))
FilterdExperimental=mean(abs(acc_filtered(24:3278+24,:)))
TransientLoading=mean(abs(accData1))
TruncatedModel50=mean(abs(accData2))
TruncatedModel100=mean(abs(accData3))
TruncatedModel200=mean(abs(accData4))

ExperimentalTot=mean(mean(abs(accDatae(24:3278+24,:))))
FilterdExperimentalTot=mean(mean(abs(acc_filtered(24:3278+24,:))))
TransientLoadingTot=mean(mean(abs(accData1)))
TruncatedModel50Tot=mean(mean(abs(accData2)))
TruncatedModel100Tot=mean(mean(abs(accData3)))
TruncatedModel200Tot=mean(mean(abs(accData4)))

figure(1) %Experimental vs Filtered Experimental
for i = 1:3
    subplot(3,1,i);
    plot(timee, accDatae(:,i), 'r-'); %Experimental
    hold on;
    plot(timee, acc_filtered(:,i), 'b-'); %Filtered Experimental
    hold on
    title(sprintf('Accelerometer %d', i));
    grid on;
    legend('Experimental','Filtered Experimental','Location','northeast'); 
    axis([0 0.4 -100 100])
end
xlabel('Time [s]')

figure(3) %Experimental vs Transient Loading
for i = 1:3
    subplot(3,1,i);
    plot(timee, accDatae(:,i), 'r-'); %Experimental
    hold on;
    plot(time1, accData1(:,i), 'b-'); %Transient Loading
    hold on;
    title(sprintf('Accelerometer %d', i));
    grid on;
    legend('Experimental','Transient Loading','Location','northeast'); 
    axis([0 0.15 -120 120])
end
xlabel('Time [s]')

figure(2) %Filterd Experimental vs Transient Loading
for i = 1:3
    subplot(3,1,i);
    plot(timee, acc_filtered(:,i), 'r-'); %Experimental acceleration
    hold on;
    plot(time1, accData1(:,i), 'b-'); %Transient Loading
    hold on;
    title(sprintf('Accelerometer %d', i));
    grid on;
    legend('Experimental Filtered','Transient Loading','Location','northeast'); 
    axis([0 0.15 -30 30])
end
xlabel('Time [s]')