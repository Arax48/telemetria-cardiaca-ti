function [f_psd, psd_nrz_db, psd_man_db] = bloque_banda_base(trama_bits)
    Rb   = 1000;          % Tasa de bit de referencia (1 kbps)
    sps  = 20;            % Muestras por bit
    Fs   = Rb * sps;      % Frecuencia de muestreo en banda base (20 kHz)
    A    = 1;             % Amplitud +-1 V

    N_bits = min(length(trama_bits), 10000);
    bits   = trama_bits(1:N_bits);

    % 1. Generacion de formas de onda NRZ Polar y Manchester
    nrz = repelem(2*bits - 1, sps) * A;

    man = zeros(N_bits * sps, 1);
    mitad = sps / 2;
    for k = 1:N_bits
        idx = (k-1)*sps + 1 : k*sps;
        if bits(k) == 1
            man(idx(1:mitad))     =  A;
            man(idx(mitad+1:end)) = -A;
        else
            man(idx(1:mitad))     = -A;
            man(idx(mitad+1:end)) =  A;
        end
    end

    % 2. Densidad Espectral de Potencia (PSD) por promedio de periodogramas
    L_seg = 1024;
    win   = 0.5 - 0.5*cos(2*pi*(0:L_seg-1)'/(L_seg-1)); % Ventana Hann
    U     = sum(win.^2);
    paso  = L_seg / 2;
    n_seg = floor((length(nrz) - L_seg) / paso) + 1;

    P_nrz = zeros(L_seg/2 + 1, 1);
    P_man = zeros(L_seg/2 + 1, 1);

    for s = 1:n_seg
        ini = (s-1)*paso + 1;
        seg_n = nrz(ini:ini+L_seg-1) .* win;
        seg_m = man(ini:ini+L_seg-1) .* win;

        Fn = abs(fft(seg_n)).^2 / (Fs * U);
        Fm = abs(fft(seg_m)).^2 / (Fs * U);

        Fn_1 = Fn(1:L_seg/2 + 1); Fn_1(2:end-1) = 2 * Fn_1(2:end-1);
        Fm_1 = Fm(1:L_seg/2 + 1); Fm_1(2:end-1) = 2 * Fm_1(2:end-1);

        P_nrz = P_nrz + Fn_1 / n_seg;
        P_man = P_man + Fm_1 / n_seg;
    end

    f_psd      = (0:L_seg/2)' * (Fs / L_seg);
    psd_nrz_db = 10 * log10(P_nrz + 1e-12);
    psd_man_db = 10 * log10(P_man + 1e-12);

    % 3. Canal limitado en banda (fc = 0.75*Rb) + AWGN -> Interferencia Intersimbolo (ISI)
    N_tot = length(nrz);
    f_2lados = abs([0:floor(N_tot/2)-1, -ceil(N_tot/2):-1]' * (Fs / N_tot));
    fc = 0.75 * Rb;
    H_canal = 1 ./ (1 + (f_2lados / fc).^6);
    nrz_canal = real(ifft(fft(nrz) .* H_canal)) + 0.08 * randn(size(nrz));

    % 4. Graficas: Formas de onda, PSD y Diagrama de Ojo
    fig = figure('Name', 'Paso 4 - Codificacion de Linea, PSD y Diagrama de Ojo', ...
        'Position', [140 80 980 680]);

    % Subplot 1: Primeros 12 bits
    subplot(3,1,1);
    n_vis = 12 * sps;
    t_ms  = (0:n_vis-1)' / Fs * 1000;
    plot(t_ms, nrz(1:n_vis), 'b', 'LineWidth', 1.5); hold on;
    plot(t_ms, man(1:n_vis), 'r--', 'LineWidth', 1.3);
    ylim([-1.4 1.4]); grid on;
    title(sprintf('Codigos de Linea (Primeros 12 bits de trama: %s)', num2str(bits(1:12)')));
    xlabel('Tiempo (ms)'); ylabel('Voltaje (V)');
    legend('NRZ Polar', 'Manchester (Autosincronizable)');

    % Subplot 2: PSD
    subplot(3,1,2);
    plot(f_psd / Rb, psd_nrz_db, 'b', 'LineWidth', 1.4); hold on;
    plot(f_psd / Rb, psd_man_db, 'r', 'LineWidth', 1.4);
    xlim([0 4]); ylim([-70 5]); grid on;
    title('Densidad Espectral de Potencia (PSD) Normalizada');
    xlabel('Frecuencia Normalizada (f / R_b)'); ylabel('PSD (dB/Hz)');
    legend('NRZ Polar (Pico en DC, nulo en 1 R_b)', 'Manchester (Sin DC, nulo en 2 R_b)');

    % Subplot 3: Diagrama de Ojo (ventanas de 2 Tb)
    subplot(3,1,3);
    win_ojo = 2 * sps;
    n_trazas = 150;
    t_ojo = linspace(0, 2, win_ojo);
    hold on;
    for k = 2:n_trazas
        traza = nrz_canal((k-1)*sps + 1 : (k-1)*sps + win_ojo);
        plot(t_ojo, traza, 'Color', [0 0.45 0.74 0.25], 'LineWidth', 1);
    end
    xline(0.5, '--k'); xline(1.5, '--k');
    ylim([-1.5 1.5]); grid on;
    title('Diagrama de Ojo en Canal Limitado en Banda (f_c = 0.75 R_b + AWGN)');
    xlabel('Tiempo Normalizado (t / T_b)'); ylabel('Amplitud (V)');

    saveas(fig, 'figs/paso4_banda_base.png');
end
