function [eb_n0_db, ber_sin_cod, ber_hamming, tabla_ber] = evaluacion_canal(trama_bits)
    % 1. Matrices del Codigo Lineal Sistematico Hamming (7,4)
    P = [1 1 0;
         1 0 1;
         0 1 1;
         1 1 1];
    G = [eye(4), P];        % Matriz generadora (4x7)
    H = [P', eye(3)];       % Matriz de verificacion de paridad (3x7)
    Ht = H';                % Tabla de sindromes (7x3)
    R = 4 / 7;              % Tasa de codigo

    % Ajustar trama_bits a multiplo de 4 y repetir para estadistica BER solida
    n_bloques_base = floor(length(trama_bits) / 4);
    bits_base = trama_bits(1 : n_bloques_base * 4);
    bits_tx   = repmat(bits_base, 12, 1); % >100,000 bits para curvas BER limpias
    N_bits    = length(bits_tx);
    n_bloques = N_bits / 4;

    % Codificacion Hamming (7,4): u (n_bloques x 4) -> c (n_bloques x 7)
    U_mat = reshape(bits_tx, 4, n_bloques)';
    C_mat = mod(U_mat * G, 2);
    bits_cod = reshape(C_mat', [], 1);

    % 2. Barrido Eb/N0 sobre Canal AWGN con modulacion BPSK antipodal (0 -> -1, 1 -> +1)
    eb_n0_db    = (0:1:10)';
    ber_sin_cod = zeros(size(eb_n0_db));
    ber_hamming = zeros(size(eb_n0_db));
    ber_teorica = zeros(size(eb_n0_db));

    tx_sim_unc = 2 * bits_tx  - 1;
    tx_sim_cod = 2 * bits_cod - 1;

    for i = 1:length(eb_n0_db)
        eb_n0_lin = 10^(eb_n0_db(i) / 10);

        % --- Canal sin codificar (energia por simbolo = Eb) ---
        sigma_unc = sqrt(1 / (2 * eb_n0_lin));
        rx_unc    = tx_sim_unc + sigma_unc * randn(size(tx_sim_unc));
        dec_unc   = rx_unc > 0;
        ber_sin_cod(i) = mean(dec_unc ~= bits_tx);
        ber_teorica(i) = 0.5 * erfc(sqrt(eb_n0_lin));

        % --- Canal con Hamming (7,4) (energia por bit codificado Ec = R * Eb) ---
        sigma_cod = sqrt(1 / (2 * R * eb_n0_lin));
        rx_cod    = tx_sim_cod + sigma_cod * randn(size(tx_sim_cod));
        bits_rx   = double(rx_cod > 0);

        % Decodificacion por sindrome S = R_mat * H' mod 2
        R_mat = reshape(bits_rx, 7, n_bloques)';
        S_mat = mod(R_mat * Ht, 2);

        % Correccion de 1 error por bloque de 7 bits
        for col = 1:7
            match = (S_mat(:,1) == Ht(col,1)) & ...
                    (S_mat(:,2) == Ht(col,2)) & ...
                    (S_mat(:,3) == Ht(col,3));
            R_mat(match, col) = mod(R_mat(match, col) + 1, 2);
        end

        U_dec = R_mat(:, 1:4);
        bits_dec = reshape(U_dec', [], 1);
        ber_hamming(i) = mean(bits_dec ~= bits_tx);
    end

    % Evitar ceros exactos en escala logaritmica si a 10 dB hay 0 errores
    ber_hamming_plot = max(ber_hamming, 1 / N_bits);

    tabla_ber = table(eb_n0_db, ber_teorica, ber_sin_cod, ber_hamming, ...
        'VariableNames', {'Eb_N0_dB', 'BER_Teorica_BPSK', 'BER_Sin_Codificar', 'BER_Hamming_7_4'});

    % 3. Capacidad de Canal de Shannon: C/B = log2(1 + SNR)
    snr_db  = -5:0.5:20;
    snr_lin = 10.^(snr_db / 10);
    cap_esp = log2(1 + snr_lin); % bits/s/Hz

    % 4. Graficas: Curvas BER vs Eb/N0 y Capacidad de Shannon
    fig = figure('Name', 'Paso 6 - Correccion de Errores Hamming (7,4) y Capacidad de Shannon', ...
        'Position', [160 110 980 500]);

    subplot(1,2,1);
    semilogy(eb_n0_db, ber_teorica, 'k--', 'LineWidth', 1.3); hold on;
    semilogy(eb_n0_db, ber_sin_cod, 'r-o', 'LineWidth', 1.4, 'MarkerFaceColor', 'r');
    semilogy(eb_n0_db, ber_hamming_plot, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on; ylim([1e-5 1]); xlim([0 10]);
    title('Desempeno BER en Canal AWGN (BPSK)');
    xlabel('E_b / N_0 (dB)'); ylabel('Tasa de Error de Bit (BER)');
    legend('BPSK Teorica Sin Cod.', 'BPSK Simulada Sin Cod.', 'Hamming (7,4) Corregido', ...
        'Location', 'southwest');

    subplot(1,2,2);
    plot(snr_db, cap_esp, 'b', 'LineWidth', 1.8); hold on;
    yline(R, '--r', sprintf('Tasa Hamming R = 4/7 (%.2f bps/Hz)', R), 'LineWidth', 1.3);
    snr_min_db = 10*log10(2^R - 1);
    xline(snr_min_db, ':k', sprintf('SNR min = %.2f dB', snr_min_db), 'LineWidth', 1.2);
    grid on; xlim([-5 20]); ylim([0 7]);
    title('Capacidad de Canal de Shannon (C/B)');
    xlabel('SNR (dB)'); ylabel('Eficiencia Espectral C/B (bits/s/Hz)');
    legend('Limite de Shannon log_2(1+SNR)', 'Tasa de Operacion Hamming', 'Location', 'northwest');

    saveas(fig, 'figs/paso6_canal_hamming.png');
end
