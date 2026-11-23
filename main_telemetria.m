clear; clc; close all;
addpath('src');

%% 1. CAPTURA DE DATOS (Contrato de interfaz Hardware ESP32-S3 <-> MATLAB)
fuente = generar_fuente_cardiaca('sintetico', '');

%% 2. ANALISIS ESPECTRAL Y FILTRADO PPG (Capitulo II)
[ppg_filtrada, W_hz] = bloque_espectral(fuente);

%% 3. MUESTREO DE NYQUIST, CUANTIZACION PCM Y BANDA BASE (Capitulo III)
[bits_pcm, tabla_sqnr]          = bloque_pcm(ppg_filtrada, fuente.Fs_sim, W_hz);
[f_psd, psd_nrz_db, psd_man_db] = bloque_banda_base(fuente.trama_bits);

%% 4. ENTROPIA DE SHANNON Y CODIFICACION HUFFMAN (Capitulo IV)
[bits_huffman, eficiencia, tabla_huff] = bloque_huffman(fuente.estados);

%% 5. CAPACIDAD DE CANAL Y CODIGO CORRECTOR HAMMING (7,4) (Capitulos IV y V)
[eb_n0_db, ber_sin_cod, ber_hamming, tabla_ber] = evaluacion_canal(fuente.trama_bits);

fprintf('\n================ RESUMEN FINAL DEL SISTEMA ================\n');
fprintf('Ancho de banda PPG (W): %.2f Hz | Nyquist min: %.2f Hz\n', W_hz, 2*W_hz);
fprintf('Eficiencia Huffman estados: %.2f%% (%d bits vs 4000 fijos)\n', eficiencia, length(bits_huffman));
disp(tabla_ber);
