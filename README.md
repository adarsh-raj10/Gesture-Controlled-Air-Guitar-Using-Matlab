# Gesture-Controlled-Air-Guitar-Using-Matlab
Gesture Controlled Air Guitar Using Matlab

go = true;
Fs = 44100;
eff = 0;
while go
a = webread('http://192.168.43.77/');
b = webread('http://192.168.43.248/');
j = a.variables.jerk;
c = b.variables.chord;
e = b.variables.effect;
if e ~= 0
eff = eff + 1;
disp("Changing effect...")
if eff == 5
eff = 0;
end
end
if c ~= 0
if j == 1
disp("Create sound...")
y = create_sound(c);
if eff ~=0
disp("Adding effects...")
y = add_effect(y, Fs, eff);
end
disp("Playing music...")
sound(y, Fs);
disp("...Done")
end
end
end
Create Sound
function y = create_music(c)
Fs = 44100;
A = 110; % The A string of a guitar is normally
tuned to 110 Hz
Eoffset = -5;
Doffset = 5;
Goffset = 10;
Boffset = 14;
E2offset = 19;
F = linspace(1/Fs, 1000, 2^12);
x = zeros(Fs*2, 1);
%For barre chords
if c == 1
% A major
fret = [5 7 7 6 5 5];
elseif c == 2
% F major
fret = [1 3 3 2 1 1];
elseif c == 3
% B major
fret = [7 9 9 8 7 7];
elseif c == 4
% E minor
fret = [0 2 2 0 0 0];
elseif c == 5
% A minor
fret = [5 7 7 5 5 5];
elseif c == 6
% C major
fret = [8 10 10 9 8 8];
elseif c == 7
% D major
fret = [10 12 12 11 10 10];
elseif c == 8
% G major
fret = [3 5 5 4 3 3];
end
delay = [round(Fs/(A*2^((fret(1)+Eoffset)/12))), ...
round(Fs/(A*2^(fret(2)/12))), ...
round(Fs/(A*2^((fret(3)+Doffset)/12))), ...
round(Fs/(A*2^((fret(4)+Goffset)/12))), ...
round(Fs/(A*2^((fret(5)+Boffset)/12))), ...
round(Fs/(A*2^((fret(6)+E2offset)/12)))];
b = cell(length(delay),1);
a = cell(length(delay),1);
H = zeros(length(delay),4096);
note = zeros(length(x),length(delay));
for indx = 1:length(delay)
% Build a cell array of numerator and
denominator coefficients.
b{indx} = firls(42, [0 1/delay(indx) 2/delay(indx)
1], [0 0 1 1]).';
a{indx} = [1 zeros(1, delay(indx)) -0.5 -0.5].';
% Populate the states with random numbers and
filter the input zeros.
zi = rand(max(length(b{indx}),length(a{indx}))-
1,1);
note(:, indx) = filter(b{indx}, a{indx}, x, zi);
% Make sure that each note is centered on zero.
note(:, indx) = note(:, indx)-mean(note(:, indx));
[H(indx,:),W] = freqz(b{indx}, a{indx}, F, Fs);
end
offset = 35;
offset = ceil(offset*Fs/1000);
for indx = 1:size(note, 2)
note(:, indx) = [zeros(offset*(indx-1),1); ...
note((1:end-offset*(indx-1)), indx)];
end
combinedNote = sum(note,2);
combinedNote = combinedNote/max(abs(combinedNote)); y
= combinedNote;
end
Add Effects
function y = add_effect(in, Fs, effect)
Fs = 44100;
distortion_gain = 11;
distortion_tone = 0.4;
flanger_mix = 0.8;
flanger_delay = 5;
flanger_width = 5;
flanger_rate = 0.4;
tremblo_rate = 4;
tremblo_depth = 0.7;
tremblo_lag = 100;
tremblo_osc = 'tri';
digdelay_depth = 0.8;
digdelay_delay = 1/3;
digdelay_feedback = 0.25;
switch effect
case 1
y =
distortion(in,distortion_gain,distortion_tone,Fs);
case 2
y =
flanger(in,flanger_mix,flanger_delay,flanger_width,flanger_
rate,Fs);
case 3
y =
tremblo(in,Fs,tremblo_rate,tremblo_depth,tremblo_lag,trembl
o_osc);
case 4
y =
digdelay(in,digdelay_depth,digdelay_delay,digdelay_feedback ,Fs);
otherwise
y = in
end
Distortion
function [out] = distortion(in,gain,tone,fs)
%DISTORTION simulates a distorted guitar
% IN - guitar input sound vector (Nx1)
% FS - sampling rate of IN
% GAIN - distortion level
% TONE - 0 -> bass , 1 -> treble
B=in;
%B=B/max(abs(B));
B=gain*B;
B=erf(B);
% From research, we chose to lowpass at 1.2kHz and highpass
at 265Hz
[ZL,PL,KL]=butter(1,1200/(fs/2),'low');
[ZH,PH,KH]=butter(1,265/(fs/2),'high');
[BL,AL]=zp2tf(ZL,PL,KL);
[BH,AH]=zp2tf(ZH,PH,KH);
% OUTPUT - filter based on position/value of tone control
out=(1-tone)*filter(BL,AL,B)+tone*filter(BH,AH,B);
end
Flanger
function [out] = flanger(in,mix,delay,width,rate,fs)
%FLANGER simulates a guitar flanger effect pedal
% IN - input vector
% MIX - depth - amt of delayed signal added to IN (0 to
1)
% DELAY - min delay time - 100usec to 10msec (in
msec) 0.1 to 10
% WIDTH - sweep depth - how wide sweep is (100nsec
to 10msec)
% (in msec, 0.0001)
% RATE - frequency of LFO - 0.05 to 5 Hz
in=interp1(1:length(in),in,1:.25:length(in))
; fsn=fs*4;
minDelaySamp=ceil(delay*fsn/1000); %convert to msec, then
samples
maxDelaySamp=ceil((delay+width)*fsn/1000); %convert to
msec, then samples
n=(1:length(in)+maxDelaySamp)'; %how long to extend in by
for LFO
LFO=sawtooth(2*pi*rate/(fsn)*n,.5); %sawtooth more commonly
used in flangers
delayTimeSamples=(delay+width/2+width/2*LFO)*fsn/1000;
% instantaneous delay in samples (computed by looking
at graph from class
% PDF)
out=zeros(length(in)+minDelaySamp,1); %initialized
output vec
out(1:maxDelaySamp)=in(1:maxDelaySamp);
% copy front of signal before min delay
for i=maxDelaySamp+1:length(in) % starting from next sample
delaySamples=ceil(delayTimeSamples(i)); %whole number
of current delay
out(i)=in(i)+mix*out(i-delaySamples); %add input and
fraction of delay
end
out=downsample(out,4);
end
Tremblo
function [ out ] = tremblo( in, fs, rate, depth, lag, LFO )
%TREMOLO simulates a guitar tremolo stomp box
% IN - input guitar vector
% FS - sampling rate of IN
% RATE - frequency of LFO
% DEPTH - amplitude of modulating LFO
% LAG - stereo lag in msec
% LFO - low frequency oscillator - sin, tri, or squ
dbstop if error
% depth 0 to 1 (amplitude 0 to 1)
% rate 0.05 to 5 Hz
lagSamples=ceil(lag/1000*fs); %converts ms of lag to #
samples
n = (1:length(in))';
argWaveLeft = 2*pi*rate/fs * n;
argWaveRight = 2*pi*rate/fs * (n+lagSamples);
% delay implemented as phase difference of the
right channel with respect
% to the left
switch LFO
case{'squ','square'}
wave_left = depth*square(argWaveLeft);
wave_right = depth*square(argWaveRight);
case{'tri','triangle'}
wave_left = depth*sawtooth(argWaveLeft,.5);
wave_right = depth*sawtooth(argWaveRight,.5);
case{'sin','sine'}
wave_left = depth*sin(argWaveLeft);
wave_right = depth*sin(argWaveRight);
end
left = in.*(1+wave_left); right =
in.*(1+wave_right); out=[left,right]; %
creating a stereo output
end
Digital Delay
function [out]=digdelay(in,depth,delay,feedback,fs)
%DIGDELAY simulates a digital delay guitar effects box
% IN - input guitar audio vector
% DEPTH - 0 (no delay sig added) to 1 (equal
amplitude delay signal)
% DELAY - delay time from 0.1 msec (0.0001 sec) to 8 sec
(in seconds)
% FEEDBACK - feedback gain - amplitude of delayed signal
to be fedback
dbstop if error
sampleDelay = ceil(delay*fs); % delay in sec * samp/sec
= samples
delayedSignal =
[zeros(sampleDelay,1);in(1:end);zeros(64*fssampleDelay,
1)];
% first delay by padding beginning with proper # zeros
% then zero pad the end of the delayed signal for the rest
of the delays
% that come from feedback
% below, in is also padded with the same # zeros to handle
longer delay
% times
in=[in;zeros(64*fs,1)];
n=1; %iterator variable for # delays
while feedback >= 0.0001 %going to continually
decrease feedback inside loop
delayedSignal=delayedSignal+...
feedback*[zeros(sampleDelay*n,1);...
delayedSignal(1:end-sampleDelay*n)];
feedback=feedback^2;
n=n+1;
end
out = in + depth*delayedSignal;
% removing excess zeros from the
end content = find(out,1,'last');
out=out(1:content);
end
