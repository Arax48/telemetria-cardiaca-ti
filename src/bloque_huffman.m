function [bits_huffman, eficiencia, tabla_huff] = bloque_huffman(estados)
    nombres = {'1: NORMAL', '2: SIN_DETECCION', '3: ALERTA_BAJA', '4: ALERTA_ALTA'};
    N = length(estados);
    simbolos = (1:4)';

    % 1. Probabilidades empiricas y Entropia de Shannon H(X)
    conteos = zeros(4, 1);
    for k = 1:4
        conteos(k) = sum(estados == k);
    end
    probs = conteos / N;

    p_nz = probs(probs > 0);
    H_X  = -sum(p_nz .* log2(p_nz)); % bits/simbolo

    % 2. Construccion del arbol binario de Huffman (MATLAB base)
    codigos = cell(4, 1);
    nodos = cell(4, 1);
    for k = 1:4
        nodos{k} = struct('prob', probs(k), 'syms', k);
    end

    while numel(nodos) > 1
        p_actual = cellfun(@(s) s.prob, nodos);
        [~, idx_ord] = sort(p_actual, 'ascend');
        nodos = nodos(idx_ord);

        n1 = nodos{1};
        n2 = nodos{2};

        for s = n1.syms
            codigos{s} = ['0', codigos{s}];
        end
        for s = n2.syms
            codigos{s} = ['1', codigos{s}];
        end

        nuevo = struct('prob', n1.prob + n2.prob, 'syms', [n1.syms, n2.syms]);
        nodos = [{nuevo}; nodos(3:end)];
    end

    % 3. Codificacion del flujo de estados y verificacion sin perdida
    cadena_bits = strjoin(codigos(estados), '');
    bits_huffman = (cadena_bits' - '0');

    % Decodificacion para verificar integridad 100% sin perdida
    estados_dec = zeros(N, 1);
    buf = '';
    pos = 1;
    for b = 1:length(cadena_bits)
        buf = [buf, cadena_bits(b)];
        for s = 1:4
            if strcmp(buf, codigos{s})
                estados_dec(pos) = s;
                pos = pos + 1;
                buf = '';
                break;
            end
        end
    end
    assert(isequal(estados, estados_dec), 'Error en decodificacion Huffman');

    % 4. Metricas de desempeno (Longitud media L, Eficiencia eta, Tasa compresion)
    longitudes = cellfun(@length, codigos);
    L_prom     = sum(probs .* longitudes);
    eficiencia = (H_X / L_prom) * 100;
    bits_fijos = N * 2; % 2 bits por estado en codificacion fija
    ahorro_pct = (1 - length(bits_huffman) / bits_fijos) * 100;

    tabla_huff = table(nombres', probs, codigos, longitudes, ...
        'VariableNames', {'Estado', 'Probabilidad', 'Codigo_Huffman', 'Longitud_Bits'});

    % 5. Grafica de probabilidades, palabras codigo y ahorro de bits
    fig = figure('Name', 'Paso 5 - Entropia de Shannon y Codificacion Huffman', ...
        'Position', [150 100 950 500]);

    subplot(1,2,1);
    bar(simbolos, probs, 0.55, 'FaceColor', [0.15 0.45 0.75]);
    grid on; ylim([0 0.95]);
    xticks(1:4);
    xticklabels({'NORMAL', 'SIN\_DET', 'AL\_BAJA', 'AL\_ALTA'});
    title(sprintf('Distribucion de Estados | Entropia H(X) = %.3f bits/simb', H_X));
    ylabel('Probabilidad Empirica p(x_i)');
    for k = 1:4
        text(k, probs(k) + 0.03, sprintf('%.1f%%\n[%s]', probs(k)*100, codigos{k}), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
    end

    subplot(1,2,2);
    categorias = categorical({'Limite Shannon\nN·H(X)', 'Huffman\n(Variable)', 'Trama Fija\n(2 bits/estado)'});
    categorias = reordercats(categorias, {'Limite Shannon\nN·H(X)', 'Huffman\n(Variable)', 'Trama Fija\n(2 bits/estado)'});
    valores_bits = [N * H_X, length(bits_huffman), bits_fijos];
    b_bar = bar(categorias, valores_bits, 0.5);
    b_bar.FaceColor = 'flat';
    b_bar.CData = [0.2 0.65 0.35; 0.15 0.45 0.75; 0.75 0.3 0.3];
    grid on; ylim([0 bits_fijos * 1.15]);
    title(sprintf('Eficiencia = %.2f%% | Ahorro vs Fijo = %.1f%%', eficiencia, ahorro_pct));
    ylabel('Bits Totales Transmitidos (N = 2000 tramas)');
    for k = 1:3
        text(k, valores_bits(k) + 120, sprintf('%d bits', round(valores_bits(k))), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    saveas(fig, 'figs/paso5_huffman.png');
end
