function fuente = generar_fuente_cardiaca(modo, ruta_csv)
    if strcmpi(modo, 'sintetico')
        Fs_sim = 1000;
        T_dur  = 10;
        t      = (0:1/Fs_sim:T_dur-1/Fs_sim)';
        
        bpm_base = 75;
        f0 = bpm_base / 60;
        
        % Senal PPG + deriva de movimiento (0.2 Hz) + ruido
        ppg_limpia    = 1.0*sin(2*pi*f0*t) + 0.45*sin(4*pi*f0*t - pi/3) + 0.15*sin(6*pi*f0*t);
        artefacto_mov = 0.35 * sin(2*pi*0.2*t);
        ruido_adc     = 0.08 * randn(size(t));
        
        fuente.t                 = t;
        fuente.Fs_sim            = Fs_sim;
        fuente.senal_limpia      = ppg_limpia;
        fuente.senal_contaminada = ppg_limpia + artefacto_mov + ruido_adc;
        
        % 2000 tramas con los 4 estados: 1=NORMAL, 2=SIN_DETECCION, 3=ALERTA_BAJA, 4=ALERTA_ALTA
        N_tramas = 2000;
        probs = [0.78, 0.05, 0.09, 0.08];
        u = rand(N_tramas, 1);
        estados = ones(N_tramas, 1);
        estados(u > probs(1)) = 2;
        estados(u > sum(probs(1:2))) = 3;
        estados(u > sum(probs(1:3))) = 4;
        
        bpm_serie = zeros(N_tramas, 1);
        for k = 1:N_tramas
            switch estados(k)
                case 1, bpm_serie(k) = randi([60, 100]);
                case 2, bpm_serie(k) = 0;
                case 3, bpm_serie(k) = randi([40, 59]);
                case 4, bpm_serie(k) = randi([101, 150]);
            end
        end
        
        bateria_pct = round(linspace(100, 65, N_tramas))';
        
        % Empaquetado binario (17 bits por trama)
        bits_estado = dec2bin(estados - 1, 2) - '0';
        bits_bpm    = dec2bin(bpm_serie, 8) - '0';
        bits_bat    = dec2bin(bateria_pct, 7) - '0';
        trama_bits  = reshape([bits_estado, bits_bpm, bits_bat]', [], 1);
        
        fuente.estados     = estados;
        fuente.bpm_serie   = bpm_serie;
        fuente.bateria_pct = bateria_pct;
        fuente.trama_bits  = trama_bits;
        
    elseif strcmpi(modo, 'hardware')
        tabla = readtable(ruta_csv);
        fuente.t                 = tabla.tiempo_s;
        fuente.Fs_sim            = round(1 / mean(diff(fuente.t)));
        fuente.senal_contaminada = tabla.adc_crudo;
        fuente.senal_limpia      = tabla.adc_filtrado;
        fuente.estados           = tabla.estado_id;
        fuente.bpm_serie         = tabla.bpm;
        fuente.bateria_pct       = tabla.bateria;
        
        bits_estado = dec2bin(fuente.estados - 1, 2) - '0';
        bits_bpm    = dec2bin(fuente.bpm_serie, 8) - '0';
        bits_bat    = dec2bin(fuente.bateria_pct, 7) - '0';
        fuente.trama_bits = reshape([bits_estado, bits_bpm, bits_bat]', [], 1);
    end
end
