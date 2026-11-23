# Sistema de Telemetría y Monitoreo Cardíaco — Teoría de la Información
Cadena completa de simulación en MATLAB y documentación en LaTeX integrada con el proyecto de Tecnología Digital (ESP32-S3 Mini / ESP32):

- **Cap. I:** Contrato de interfaz Hardware <-> MATLAB (`src/generar_fuente_cardiaca.m`).
- **Cap. II:** Análisis espectral FFT, ancho de banda efectivo y filtrado pasa-banda PPG (`src/bloque_espectral.m`).
- **Cap. III:** Muestreo de Nyquist, cuantización PCM (`SQNR`), códigos de línea NRZ/Manchester y diagrama de ojo (`src/bloque_pcm.m`, `src/bloque_banda_base.m`).
- **Cap. IV:** Entropía de Shannon y compresión sin pérdida Huffman sobre estados clínicos (`src/bloque_huffman.m`).
- **Cap. V:** Capacidad de canal y corrección de errores Hamming (7,4) en canal AWGN (`src/evaluacion_canal.m`).
