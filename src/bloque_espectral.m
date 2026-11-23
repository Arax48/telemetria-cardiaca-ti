function [ppg_filtrada, W_hz] = bloque_espectral(fuente)
    Fs = fuente.Fs_sim;
    x  = fuente.senal_contaminada;
    N  = length(x);

    % 1. Espectro unilateral mediante FFT
    X = fft(x);
    P2 = abs(X / N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0:floor(N/2))' / N;

    % 2. Calculo de Ancho de Banda Efectivo W (99% de potencia en banda cardiaca 0.5-40 Hz)
    idx_banda = (f >= 0.5) & (f <= 40);
    potencia = P1(idx_banda).^2;
    f_banda  = f(idx_banda);
    pot_acum = cumsum(potencia) / sum(potencia);
    idx_W    = find(pot_acum >= 0.99, 1, 'first');
    W_hz     = f_banda(idx_W);

    % 3. Filtro Pasa-Banda Butterworth fase cero (0.5 Hz a W_hz) en frecuencia
    f_dos_lados = [0:floor(N/2)-1, -ceil(N/2):-1]' * (Fs / N);
    f_abs = abs(f_dos_lados);
    orden = 4;
    H_hp = zeros(size(f_abs));
    nz = f_abs > 0;
    H_hp(nz) = 1 ./ (1 + (0.5 ./ f_abs(nz)).^(2*orden));
    H_lp = 1 ./ (1 + (f_abs ./ max(W_hz, 8)).^(2*orden));
    H_bp = H_hp .* H_lp;

    ppg_filtrada = real(ifft(X .* H_bp));

    % Espectro de la senal filtrada
    Xf = fft(ppg_filtrada);
    P1_f = abs(Xf(1:floor(N/2)+1) / N);
    P1_f(2:end-1) = 2 * P1_f(2:end-1);

    % 4. Grafica comparativa Tiempo y Frecuencia
    fig = figure('Name', 'Paso 2 - Analisis Espectral y Filtrado', 'Position', [100 100 950 520]);

    subplot(2,1,1);
    plot(fuente.t, x, 'Color', [0.65 0.65 0.65]); hold on;
    plot(fuente.t, ppg_filtrada, 'b', 'LineWidth', 1.3);
    xlim([0 4]); grid on;
    title('Senal PPG en el Tiempo: Contaminada vs. Filtrada');
    xlabel('Tiempo (s)'); ylabel('Amplitud');
    legend('Cruda (con deriva 0.2 Hz + ruido)', 'Filtrada Pasa-Banda');

    subplot(2,1,2);
    plot(f, P1, 'Color', [0.65 0.65 0.65]); hold on;
    plot(f, P1_f, 'r', 'LineWidth', 1.3);
    xline(W_hz, '--k', sprintf('W = %.2f Hz (99%% pot.)', W_hz), 'LineWidth', 1.2);
    xlim([0 15]); grid on;
    title('Espectro Unilateral |X(f)|');
    xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
    legend('Espectro Crudo', 'Espectro Filtrado');

    saveas(fig, 'figs/paso2_espectral.png');
end
