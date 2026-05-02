% Data inladen

% Zonder de .[...] krijg je een 1x1 struct, waarop je fieldnames()
% kan oproepen. Die fieldnames onthullen dat de structs maar 1 fieldname
% hebben, namelijk die met de naam van het bestand zelf.

dayssn = load('dayssn.mat').dayssn; 
yearssn = load('yearssn.mat').yearssn;

%% Day data

% Eerste kolom: jaartal
% Tweede kolom: data

Xday = fft(dayssn(:,2));
Xday = abs(Xday);                   % Compute amplitude
nyqday = floor(length(Xday)/2);     % Nyquist frequency

%%
% Plot frequency spectrum (lower than Nyquist frequency)

% For real inputs, the FFT output is SYMMETRIC!!
%
%       X_{n-k}  =  X_k^*
%      |X_{n-k}| = |X_k^*|
%
% The second half of the (complex) coefficients vector 
% is just a mirror image of the first half (with a complex conjugate twist). 
% 
% That’s why we usually throw away the second half and 
% only look up to the Nyquist frequency.


% Plot up to the Nyquist freq.
figure;
plot(Xday(2:nyqday)); % Only frequencies up to Nyquist frequency (no DC)
title('Frequency spectrum of day data (up to Nyquist)');
xlabel('Frequency');
ylabel('Amplitude');


% Het laatste datapunt (n-1!) komt overeen met
% de eerste AC-component!
%
%       |X_{n-k}| = |X_k^*|
%    => |X_{n-1}| = |X_1^*|
%
% De DC-component wordt dus NIET gespiegeld.


% Plot WHOLE frequency spectrum
figure;
plot(Xday(2:end-1)); % Only frequencies up to Nyquist frequency (no DC)
title('WHOLE Frequency spectrum of day data');
xlabel('Frequency');
ylabel('Amplitude');


%% Vind periode

numMeasuredDays = length(Xday)
numMeasuredYears = numMeasuredDays / 365

% Find period of max ampl.
[~, indexday] = max(Xday(2:nyqday));    % Find index of dominant frequency

% Over 100 jaar 14 cycli
numCyclesCompletedOverMeasuredDays = indexday

daysPerCycle = numMeasuredDays / numCyclesCompletedOverMeasuredDays; 
yearsPerCycle = daysPerCycle / 365;

fprintf('Period computed from day data is %.2f years.\n', yearsPerCycle); % 10.63

%% Year data

% Eerste kolom: jaartal
% Tweede kolom: data per jaar

Xyear = fft(yearssn(:,2));
Xyear = abs(Xyear);                   % Compute amplitude
nyqyear = floor(length(Xyear)/2);     % Nyquist frequency

%% Plot frequency spectrum (Year data)

figure;
plot(Xyear(2:nyqyear)); 
title('Frequency spectrum of year data (up to Nyquist)');
xlabel('Frequency (Cycles per total dataset length)');
ylabel('Amplitude');

%% Vind periode (Year data)

% LET OP DE EENHEID: Elk datapunt is nu een jaar, geen dag!
numMeasuredYears_YearData = length(Xyear); 

% Find period of max ampl. (Sla DC-component over)
[~, indexyear] = max(Xyear(2:nyqyear));    

% Dit is weer het aantal cycli dat in de totale dataset past:
numCyclesCompletedOverMeasuredYears = indexyear;

% Bereken lengte van 1 cyclus direct in jaren:
yearsPerCycle_YearData = numMeasuredYears_YearData / numCyclesCompletedOverMeasuredYears; 

fprintf('Period computed from year data is %.2f years.\n', yearsPerCycle_YearData);