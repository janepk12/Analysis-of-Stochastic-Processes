clc;
clear;
close all;
format long g;

%% ---------------------- Load TLS and IBIS-S data ------------------------
tls_data  = load('tls.txt');   % [time, displacement]
ibis_data = load('ibis.txt');  % [time, displacement]

t_tls = tls_data(:,1);
y_tls = tls_data(:,2);

t_ibis = ibis_data(:,1);
y_ibis = ibis_data(:,2);

% sampling interval (both ~0.005 s)
dt_tls  = t_tls(2)  - t_tls(1);
dt_ibis = t_ibis(2) - t_ibis(1);

%% ---------------------- Figure 1: raw TLS and IBIS ----------------------
figure(1); clf;
plot(t_tls,  y_tls,  '.-');
hold on;
plot(t_ibis, y_ibis, '.-');
legend('TLS observations','IBIS-S observations','Location','best');
title('TLS and IBIS-S displacement time series');
xlabel('Time [s]');
ylabel('Displacement [m]');
grid on;

%% --------- Cross-correlation TLS vs IBIS to find time shift ------------
% Demeaned and normalized
y_tls_dm  = y_tls  - mean(y_tls);
y_ibis_dm = y_ibis - mean(y_ibis);

[xcf_tls_ibis, lags_tls_ibis] = xcorr(y_tls_dm, y_ibis_dm, 'normalized');

[~, idx_max]  = max(abs(xcf_tls_ibis));
lag_tls_ibis  = lags_tls_ibis(idx_max);  % in samples
corr_tls_ibis = xcf_tls_ibis(idx_max);
time_shift_tls = lag_tls_ibis * dt_tls;  % in seconds

fprintf('TLS vs IBIS-S:\n');
fprintf('  lag  = %d samples\n', lag_tls_ibis);
fprintf('  corr = %.3f\n', corr_tls_ibis);
fprintf('  time shift = %.3f s\n\n', time_shift_tls);

%% --------------- Figure 2: cross-correlation TLS–IBIS -------------------
figure(2); clf;
plot(lags_tls_ibis, xcf_tls_ibis, '.');
xlabel('Shift / lags [samples]');
ylabel('Cross-correlation value');
title('Cross-correlation between TLS and IBIS-S');
grid on;
hold on;
plot(lag_tls_ibis, corr_tls_ibis, 'ko', 'MarkerFaceColor','k');
hold off;

%% --------------- Correct TLS time axis using time_shift_tls -------------
t_tls_corr = t_tls + time_shift_tls;

figure(3); clf;
plot(t_ibis, y_ibis, '-','DisplayName','IBIS-S');
hold on;
plot(t_tls,      y_tls, ':','DisplayName','TLS original');
plot(t_tls_corr, y_tls, '-','DisplayName','TLS corrected');
xlabel('Time [s]');
ylabel('Displacement [m]');
title('TLS corrected by cross-correlation (reference: IBIS-S)');
legend('Location','best');
grid on;
hold off;

%% ---------------------- Load accelerometer data -------------------------
acc_data = load('acceleration.txt');   % [time, acceleration in g]
t_acc = acc_data(:,1);
y_acc = acc_data(:,2);

dt_acc = t_acc(2) - t_acc(1);

%% ------ Figure 4: IBIS-S, corrected TLS and raw accelerometer ----------
figure(4); clf;

yyaxis left;
plot(t_ibis, y_ibis, '-', 'DisplayName','IBIS-S');
hold on;
plot(t_tls_corr, y_tls, '-', 'DisplayName','TLS corrected');
ylabel('Displacement [m]');

% --- Removing negative times of accelerometer ---
idx_acc_plot = t_acc >= 0;
yyaxis right;
plot(t_acc(idx_acc_plot), y_acc(idx_acc_plot), '-', ...
    'DisplayName','Accelerometer (raw)');
ylabel('Acceleration [g]');

xlabel('Time [s]');
title('IBIS-S, corrected TLS and raw accelerometer measurements');
legend('Location','best');
grid on;
hold off;


%% -------- Time shift of accelerometer relative to IBIS-S ----------------
% Work only on the overlapping time interval
t_start = max( t_ibis(1), t_acc(1) );
t_end   = min( t_ibis(end), t_acc(end) );

idx_acc  = (t_acc  >= t_start) & (t_acc  <= t_end);
idx_ibis = (t_ibis >= t_start) & (t_ibis <= t_end);

t_acc_sub  = t_acc(idx_acc);
y_acc_sub  = y_acc(idx_acc);

% Resample IBIS-S onto the accelerometer time grid (no extrapolation)
y_ibis_resamp = interp1(t_ibis(idx_ibis), y_ibis(idx_ibis), t_acc_sub, 'linear');

% Invert accelerometer axis to be consistent with TLS/IBIS
y_acc_sub_inv = -y_acc_sub;

% Demean both signals
y_acc_dm      = y_acc_sub_inv      - mean(y_acc_sub_inv);
y_ibis_res_dm = y_ibis_resamp      - mean(y_ibis_resamp);

[xcf_acc_ibis, lags_acc_ibis] = xcorr(y_acc_dm, y_ibis_res_dm, 'normalized');

[~, idx_acc_max] = max(abs(xcf_acc_ibis));
lag_acc_ibis     = lags_acc_ibis(idx_acc_max);   % samples on acc grid
corr_acc_ibis    = xcf_acc_ibis(idx_acc_max);
time_shift_acc   = lag_acc_ibis * dt_acc;       % seconds

fprintf('Accelerometer vs IBIS-S (overlapping interval):\n');
fprintf('  lag  = %d samples (accelerometer grid)\n', lag_acc_ibis);
fprintf('  corr = %.3f\n', corr_acc_ibis);
fprintf('  time shift = %.3f s\n\n', time_shift_acc);

%% ----------- Figure 5: cross-correlation ACC–IBIS (resampled) ----------
figure(5); clf;
plot(lags_acc_ibis * dt_acc, xcf_acc_ibis, '.');
xlabel('Shift / lags [s]');
ylabel('Cross-correlation value');
title('Cross-correlation between accelerometer and IBIS-S (resampled)');
grid on;
hold on;
plot(lag_acc_ibis * dt_acc, corr_acc_ibis, 'ko', 'MarkerFaceColor','k');
hold off;

%% -------- Correct accelerometer time axis using time_shift_acc ---------
t_acc_corr = t_acc + time_shift_acc;
y_acc_corr = -y_acc;   % keep inverted sign also for full series

%% -------- Figure 6: All three time series corrected --------------------
figure(6); clf;

yyaxis left;
plot(t_ibis,     y_ibis, '-', 'DisplayName','IBIS-S');
hold on;
plot(t_tls_corr, y_tls,  '-', 'DisplayName','TLS corrected');
ylabel('Displacement [m]');

% ---- Removing negative times from corrected acc ----
idx_acc_plot2 = t_acc_corr >= 0;

yyaxis right;
plot(t_acc_corr(idx_acc_plot2), y_acc_corr(idx_acc_plot2), '-', ...
     'DisplayName','Accelerometer corrected');
ylabel('Acceleration [g] (sign inverted)');

xlabel('Time [s]');
title('IBIS-S, TLS corrected and accelerometer corrected');
legend('Location','best');
grid on;
hold off;