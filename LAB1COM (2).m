A = 1; %AMPLITUD
fc = 1000; %FRECUENCIA
T = 1/100000; %TIEMPO DE MUESTRAS
t = 0:T:0.003; %VECTOR DE TIEMPO
L = 200;

x = A * sin(2 * pi * fc * t); %GENERACION ONDA SINU

 

fs = 5000; %FRECUENCIA MUESTREO
TPAM = 1/5000; %PERIODO MUESTREO PAM
d = 0.5; %CICLO TRABAJO

%GENERAR SEÑAL PAM NATURAL
pam = (mod(t, TPAM) < (d * TPAM));
x_pam = x .* pam;

%GENERAL SEÑAL PAM INSTANaTNEA
fs = 5000;
y = [0, diff(floor(square(2*pi*fs*t,100*d)+1)/2)] == 1;
INS = x .* y;


tiledlayout(3,1);

figure;

hold on;
xlabel('Segundos');
ylabel('Amplitud');
title('Señales');
plot(t,x, 'b', 'LineWidth', 0.5);
plot(t, x_pam, 'g', 'LineWidth', 0.5);
plot(t, INS, 'r', 'LineWidth', 0.5);
legend('Onda Seno', 'Señal PAM Nat.', 'Señal PAM Inst.'); 
hold off;


%FFT DE ONDAS
%FFT ONDA SENO ORIGINAL
L = length(t);

XFFT = fft(x);
magnitud = abs(XFFT / L);
magnitud = magnitud(1:L/2+1);
%FFT MUESTREO NATURAL
FFT_pam = fft(x_pam);
TNAT = abs(FFT_pam / L);
TNAT = TNAT(1:L/2+1);
%FFT MUESTREO INSTANTANEO
FFT_ins = fft(INS);
TINS = abs(FFT_ins / L);
TINS = TINS(1:L/2+1);

Fs = 1/T;         
f = Fs * (0:(L/2)) / L;  


figure;
subplot(3,1,1);
plot(f, TNAT);
title('T Fourier MUESTREO NATURAL');
xlabel('Frecuencia');
ylabel('Magnitud');

subplot(3,1,2);
plot(f, magnitud);
title('T. Fourier Original');
xlabel('Frecuencia');
ylabel('Magnitud');

subplot(3,1,3);
plot(f, TINS);
title('T. Fourier INSTANTANEA');
xlabel('Frecuencia');
ylabel('Magnitud');

%PARTE 2 PCM de señal INSTANTANEA
N = 64; % Número de bits para PCM
pcm_levels = 2^N; % Total de niveles PCM
% Cuantizar la señal instantánea usando PCM 
%SEÑAL ORIGINAL: x NATURAL: x_pam INSTANTANEA: INS
pcm_signal_INS = round((INS + 1) * (pcm_levels - 1) / 2); % Cuantización

% Normaliza las señales para que estén en el mismo rango de amplitud (0 a 1)
x_norm = (x - min(x)) / (max(x) - min(x));
INS_norm = (INS - min(INS)) / (max(INS) - min(INS));
pcm_signal_INS_norm = (pcm_signal_INS - min(pcm_signal_INS))/(max(pcm_signal_INS) - min(pcm_signal_INS));
quantization_error_INS = INS - ((2 * pcm_signal_INS / (pcm_levels - 1)) -1);


% Graficar la señal original en azul
figure;
plot(t, x_norm, 'b', 'LineWidth', 1.5);
hold on;

plot(t, INS_norm, 'r', 'LineWidth', 1.5);

stem(t, pcm_signal_INS_norm, 'g', 'Marker', 'o', 'LineWidth', 1.5);

xlabel('Tiempo (s)');
ylabel('Amplitud Normalizada');
title('Señal Original, Señal PAM Instantánea y Señal PAM Cuantificada (PCM)');
legend('Señal Original', 'Señal PAM Instantánea', 'Señal PAM Cuantificada (PCM)');
grid on;

figure;

plot(t, quantization_error_INS, 'k--', 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Error de Cuantización');
title('Error de Cuantización (PCM)');
grid on;