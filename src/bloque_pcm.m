function [bits_pcm, tabla_sqnr] = bloque_pcm(ppg_filtrada, Fs_sim, W_hz)
    % 1. Muestreo valido (Fs = 50 Hz > 2*W_hz) vs Submuestreo con Aliasing (Fs_alias = 4 Hz < 2*W_hz)
    Fs_ok    = 50;
    Fs_alias = 4;

    paso_ok    = round(Fs_sim / Fs_ok);
    paso_alias = round(Fs_sim / Fs_alias);

    t_orig  = (0:length(ppg_filtrada)-1)' / Fs_sim;
    x_ok    = ppg_filtrada(1:paso_ok:end);
    t_ok    = t_orig(1:paso_ok:end);

    x_alias = ppg_filtrada(1:paso_alias:end);
    t_alias = t_orig(1:paso_alias:end);

    % 2. Cuantizacion uniforme y SQNR para b = [4, 6, 8, 10, 12] bits
    bits_lista = [4, 6, 8, 10, 12]';
    sqnr_exp   = zeros(size(bits_lista));
    sqnr_teo   = 6.02 * bits_lista + 1.76;

    x_min = min(x_ok);
    x_max = max(x_ok);
    P_senal = mean(x_ok.^2);

    bits_pcm = [];

    for i = 1:length(bits_lista)
        b = bits_lista(i);
        L = 2^b;
        delta = (x_max - x_min) / L;

        % Indices de cuantizacion [0, L-1]
        idx_q = floor((x_ok - x_min) / delta);
        idx_q(idx_q >= L) = L - 1;
        idx_q(idx_q < 0)  = 0;

        % Nivel reconstruido (punto medio del intervalo)
        x_q = x_min + (idx_q + 0.5) * delta;

        error_q = x_ok - x_q;
        P_error = mean(error_q.^2);
        sqnr_exp(i) = 10 * log10(P_senal / P_error);

        % Guardar flujo binario PCM para b = 8 bits
        if b == 8
            matriz_bits = dec2bin(idx_q, 8) - '0';
            bits_pcm = reshape(matriz_bits', [], 1);
        end
    end

    tabla_sqnr = table(bits_lista, sqnr_exp, sqnr_teo, ...
        'VariableNames', {'Bits_b', 'SQNR_Exp_dB', 'SQNR_Teorica_dB'});

    % 3. Graficas de Muestreo (Nyquist vs Aliasing) y Curva SQNR
    fig = figure('Name', 'Paso 3 - Muestreo Nyquist y Cuantizacion PCM', 'Position', [120 120 950 520]);

    subplot(2,1,1);
    plot(t_orig, ppg_filtrada, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
    stem(t_ok, x_ok, 'b', 'MarkerSize', 4, 'LineWidth', 1);
    plot(t_alias, x_alias, 'r-o', 'LineWidth', 1.4, 'MarkerFaceColor', 'r');
    xlim([0 4]); grid on;
    title(sprintf('Teorema de Nyquist (2W = %.2f Hz): Muestreo Valido vs. Aliasing', 2*W_hz));
    xlabel('Tiempo (s)'); ylabel('Amplitud');
    legend('PPG Continua', sprintf('Fs = %d Hz (Sin Aliasing)', Fs_ok), ...
        sprintf('Fs = %d Hz (Con Aliasing)', Fs_alias));

    subplot(2,1,2);
    plot(bits_lista, sqnr_teo, 'k--s', 'LineWidth', 1.3, 'MarkerSize', 7); hold on;
    plot(bits_lista, sqnr_exp, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on;
    title('Desempeno de Cuantizacion Uniforme PCM: SQNR Experimental vs. Teorica');
    xlabel('Resolucion del ADC (bits b)'); ylabel('SQNR (dB)');
    xticks(bits_lista);
    legend('SQNR Teorica (6.02b + 1.76)', 'SQNR Experimental PPG', 'Location', 'northwest');

    saveas(fig, 'figs/paso3_pcm.png');
end
