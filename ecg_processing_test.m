% Corify Care | Javier Soriano Zorio

load("ecgConditioningExample.mat");

%% CARACTERIZACIÓN

tic
n_channels = size(ecg,2); % nº canales ECG
ecg_points = size(ecg,1); % nº puntos ECG
ecg_format = class(ecg); % formato de los datos
expected_format = 'double'; % comprobación de formato
if ~isa(ecg, expected_format)
    error('ECG format is not expected. Must be double');
end

discon_channels = false(1,n_channels); % canales desconectados
mean_channels = mean(ecg);
mean_threshold_up = 0.2*median(mean_channels);
mean_threshold_dw = 2*median(mean_channels);
for i=1:n_channels
    if mean_channels(i)<mean_threshold_up || mean_channels(i)>mean_threshold_dw
        discon_channels(i) = true;
        fprintf('WARNING: Channel %d is disconnected \n',i);
    end
end
ecg(:,discon_channels) = [];
n_channels = size(ecg,2);

% Se ha utilizado la media porque la varianza no detecta los canales
% desconectados con valores muy elevados. 
% Después de identificar los canales desconectados, se han eliminado para
% mejorar el tiempo computacional lo máximo posible.
% Se ha generado el vector de canales donde se pueden ver los canales
% desconectados además del texto en la consola.


%% FILTRADO

spike_threshold = 5;  % recording spikes
for i = 1:n_channels
    ecg(:, i) = filloutliers(ecg(:, i), 'linear', 'ThresholdFactor', spike_threshold);
end

lp_freq = 60;  % low-pass filter
[lpB,lpA] = butter(4,lp_freq/(fs/2),'low');  
ecg_lp = filtfilt(lpB,lpA,ecg);

notch_freq = 50;  % notch filter
q_factor = 30;   
wo = notch_freq/(fs/2);
bw = wo/q_factor;
[nB, nA] = iirnotch(wo, bw);
ecg_notch = filtfilt(nB, nA, ecg_lp);

hp_freq = 0.5;  % baseline wander removal
[hpB, hpA] = butter(2, hp_freq/(fs/2), 'high');  
ecg_bw = filtfilt(hpB, hpA, ecg_notch);

% En la eliminación de spikes se ha utilizado el valor 5 en el umbral porque 
% es un valor estándar y da un resultado muy bueno.
% Una alternativa para mejorar el tiempo computacional es hacer el paso bajo en valores 
% menores de 50Hz y no así ahorrarse el notch, pero se ha preferido hacer los dos porque 
% estaba en las instrucciones.
% La eliminación de la línea base se puede hacer también restando la media
% pero el paso alto tenía mejor resultado.

%% VERIFICACIÓN

for i=1:n_channels % representación de resultados del pipeline
    figure
    grid on
    subplot(2,2,1)
    plot(ecg(:,i))
    title('Raw ECG data')
    subplot(2,2,2)
    plot(ecg_lp(:,i))
    title('Low-pass filtered ECG')
    subplot(2,2,3)
    plot(ecg_notch(:,i))
    title('Low-pass + notch ECG')
    subplot(2,2,4)
    plot(ecg_bw(:,i))
    title('filtered + baseline wander removed ECG')
end

comp_time = toc; % medida del coste computacional
fprintf('Computational time = %.2f \n',comp_time);
fprintf('Data format = %s \n', ecg_format);




