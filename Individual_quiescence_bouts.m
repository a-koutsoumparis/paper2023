clc;
close all;
clear;

% Read position data
a = readtable('Sample traces.xlsx','Sheet',1);
speed = table2array(a(:,1));
% If each worm has a different number of frames, matlab fills in NaN
speed = speed(~isnan(speed));

% Time interval between images
dt1 = 5;

% Duration of the movie in frames
movie_frames = size(speed,1);

t1 = movie_frames * dt1;

% Set time axis
time = (dt1:dt1:t1)';

% Smooth speed according to frame density
if 1/t1 <= 10^-4
    s_speed = smooth(speed,0.01,'rloess');
else
    s_speed = speed;
end

% Normalize speed at the 0.5th quantile
n_speed = zeros(length(s_speed),1);
for j=1:size(s_speed,1)
        n_speed(j,1) = (s_speed(j,1) - quantile(s_speed,0.005))/(quantile(s_speed,0.995)-quantile(s_speed,0.005));
end

% Create list of indexes where normalized speed is under a threshold
index = zeros(size(n_speed,1),1);
threshold = 0.4; % Arbitrary, but never above 50%
for k=1:size(n_speed,1)
    if n_speed(k) <= threshold
        index(k) = k;
    else 
        index(k) = 0;
    end
end

% Get the time points that correspond to the indexes
time_points = zeros(size(index,1),1);
for l=1:size(index,1)
    if index(l) ~= 0
        time_points(l) = time(l);
    end
end

% Specify quiescence bout time period of at least 120 sec
quiet_frames = 120/dt1;
quiescence_bouts = zeros(size(time_points,1)-quiet_frames,1);
for m=1:size(time_points,1)-quiet_frames
    if sum(time_points(m:(m+quiet_frames))) == sum(time(m:(m+quiet_frames))) % Continuity condition
        quiescence_bouts(m:(m+quiet_frames)) = time_points(m:(m+quiet_frames));
    end
end

% Add extra 0 at the end of the interval to allow for limit calculations
quiescence_bouts = [quiescence_bouts;0];

% Find the limits of quiescence bout intervals
limits = zeros(size(quiescence_bouts,1),1);

% In case there is a quiescence bout in the beginning
if quiescence_bouts(1) > 0
    limits(1) = quiescence_bouts(1);
end

for n=2:size(quiescence_bouts,1)
    if quiescence_bouts(n) > 0 && quiescence_bouts(n-1) == 0
        limits(n) = quiescence_bouts(n);
    elseif quiescence_bouts(n) > 0 && quiescence_bouts(n+1) == 0
        limits(n) = quiescence_bouts(n);
    end
end

% Select for non-zero values - put help column at end
limits_elements = nonzeros(limits);

% Reshape into matrix (Row 1 = Beginning of bout, Row 2 = End of bout, Columns = Different bouts)
limits_elements = reshape(limits_elements, [2,length(limits_elements)/2]);

% Add extra 0 column at the end to allow plotting
limits_elements = [limits_elements zeros(2,1)];

% Total quiescence time (minutes)
quiet_time = zeros(size(limits_elements,2)-1,1);
if size(limits_elements,2) > 1
    for q=1:size(limits_elements,2)-1
        quiet_time(q) = (limits_elements(2,q) - limits_elements(1,q))/60;
    end
else quiet_time = 0; % If there are no quiescence bouts
end

% Total wake time (minutes)
wake_time = zeros(size(limits_elements,2)-2,1);
if size(limits_elements,2) > 2
    for z=1:size(limits_elements,2)-2
        wake_time(z) = (limits_elements(1,z+1) - limits_elements(2,z))/60;
    end
else wake_time = NaN; % If there are no quiescence bouts
end

if mean(mean(limits_elements)) == 0 % In case there are no quiescence bouts
    wake_time = t1/60;
elseif max(time)-limits_elements(2,end-1) == 0 % In case movie ends in a quiescence bout
   wake_time = [(limits_elements(1,1)-min(time))/60;wake_time];
else
    wake_time = [(limits_elements(1,1)-min(time))/60;wake_time;(max(time)-limits_elements(2,end-1))/60];
end

% Quiescence fraction
Fraction_sleep = (sum(quiet_time)/(t1/60)) * 100;

% Quiescent bout frequency
if quiet_time == 0
    Frequency_quiet = 0;
else
    Frequency_quiet = size(quiet_time,1);
end

% Wake bout frequency
Frequency_wake = size(wake_time,1);

% Average quiescence duration (minutes)
Quiet_average = mean(quiet_time,'omitnan');

% Average wake duration (minutes)
Wake_average = mean(wake_time,'omitnan');

% Average speed in quiescence bouts
if size(limits_elements,2) > 1
    for i=1:size(limits_elements,2)-1
        m_s{i} = s_speed(limits_elements(1,i)/dt1:limits_elements(2,i)/dt1);
    end
else m_s = {NaN}; % If there are no quiescence bouts (Quiet time = 0 --> Quiet speed = NaN)
end

for i=1:size(m_s,2)
    m_s_a(i,1) = mean(m_s{i});
end

avg_sleep_speed = mean(m_s_a,'omitnan');

% Average speed in wake
if size(limits_elements,2) > 2
    for z=1:size(limits_elements,2)-2
        m_w{z} = s_speed(limits_elements(2,z)/dt1:limits_elements(1,z+1)/dt1);
    end
else m_w = {NaN}; % If there are no or only 1 quiescence bout
end

for z=1:size(m_w,2)
    m_w_a(z,1) = mean(m_w{z});
end

% In case there are no quiescence bouts
if mean(mean(limits_elements)) == 0 % In case there are no quiescence bouts
    m_w_a = mean(s_speed);
elseif mean(s_speed(limits_elements(2,end-1)/dt1:max(time)/dt1)) == s_speed(end) % In case movie ends in a quiescence bout
   m_w_a = [mean(s_speed(min(time)/dt1:limits_elements(1,1)/dt1));m_w_a;];
else
   m_w_a = [mean(s_speed(min(time)/dt1:limits_elements(1,1)/dt1));m_w_a;mean(s_speed(limits_elements(2,end-1)/dt1:max(time)/dt1))];
end

avg_wake_speed = mean(m_w_a,'omitnan');

% Speed matrix - Concatenate 2 matrices of different size
lenA = length(m_w_a)*2;
lenB = length(m_s_a)*2;
C = NaN((lenA + lenB),2);
C(1:2:lenA,1) = m_w_a;
C(1:2:lenB,2) = m_s_a;
for i=1:size(C,1)/2
    I(i,1) = 2*i;
end
C(I,:) = [];

% Time matrix - Concatenate 2 matrices of different size
lenD = length(wake_time)*2;
lenE = length(quiet_time)*2;
F = NaN((lenD + lenE),2);
F(1:2:lenD,1) = wake_time;
F(1:2:lenE,2) = quiet_time;
for i=1:size(F,1)/2
    J(i,1) = 2*i;
end
F(J,:) = [];

% Plot the data
figure(1)
% plot(time./3600,speed,'kx')
% hold on
plot(time./3600,s_speed,'r-','LineWidth',1)
xlim([0 max(time)./3600])
xticks(0:0.5:max(time)/3600)
% ylim([min(speed) max(speed)])
ylim([0 500])
yticks(0:50:500)
xlabel('Time (hours)')
ylabel('ImSub intensity (arb. unit)')
title('Wild type speed during L1 arrest (48h) - Sample trace')

% Add a patch
shade = [];
for p=1:size(limits_elements,2)
    shade = patch([limits_elements(1,p)/3600 limits_elements(2,p)/3600 limits_elements(2,p)/3600 limits_elements(1,p)/3600],[0+2 0+2 500-2 500-2],[0.1 1 0.1],'EdgeColor',[0.9 1 0.9]);
    alpha(shade, 0.1)
end

TABLE = [Fraction_sleep avg_wake_speed avg_sleep_speed Wake_average Quiet_average Frequency_wake Frequency_quiet];