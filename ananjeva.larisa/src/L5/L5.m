clear; clc; close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'materials', 'figs');
resdir = fullfile(root, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end

K = 200;
T = 0.2;
w0 = 1/T;
Mreq = 1.5;
Sreq = 37;

num0 = [K]; den0 = [T 1 0];

wn0 = sqrt(K/T); zeta0 = 1/(2*sqrt(K*T));
M0 = 1/(2*zeta0*sqrt(1 - zeta0^2));
sigma0 = 100*exp(-pi*zeta0/sqrt(1 - zeta0^2));

zetaD = sqrt((1 - sqrt(1 - 1/Mreq^2))/2);
b = 4*K*zetaD^2;
natFreqD = sqrt(K*b);
numD = [K]; denD = [1/b 1 0];

numK = [T 1];
denK = [1/b 1];
alpha = b*T; Khf = alpha;
Tzero = T; Tpole = 1/b;

W0 = tf(num0, den0);
F0 = feedback(W0, 1);
Wd = tf(numD, denD);
Fd = feedback(Wd, 1);
Wk = tf(numK, denK);
Wc = minreal(Wk*W0);
Fcorr = minreal(feedback(Wc, 1));

openSys = {W0, Wd, Wc}; closedSys = {F0, Fd, Fcorr};
systemsNames = {'original', 'desired', 'corrected'};

n = numel(openSys);
gainMargin = zeros(1, n);
phaseMargin = zeros(1, n);
phaseCrossFreq = nan(1, n);
gainCrossFreq = zeros(1, n);

for k = 1:n
    [gainMargin(k), phaseMargin(k), phaseCrossFreq(k), gainCrossFreq(k)] = margin(openSys{k});
end
gainMargindB = 20*log10(gainMargin);

stable = zeros(1, n); Mpeak = zeros(1, n); wPeak = zeros(1, n);
Tsettle = zeros(1, n); Sig = zeros(1, n);
for k = 1:n
    pk = pole(closedSys{k});
    stable(k) = double(all(real(pk) < 0));
    [Mpeak(k), wPeak(k)] = getPeakGain(closedSys{k});
    info = stepinfo(closedSys{k});
    Tsettle(k) = info.SettlingTime; Sig(k) = info.Overshoot;
end

fid = fopen(fullfile(resdir, 'summary.csv'), 'w');
fprintf(fid, 'system,K,stable,wc,pm_deg,gm_db,M,wPeak,Tsettle,sigma\n');
for k = 1:n
    fprintf(fid, '%s,%.6g,%d,%.6f,%.6f,%.6g,%.6f,%.6f,%.6f,%.6f\n', ...
        systemsNames{k}, K, stable(k), gainCrossFreq(k), phaseMargin(k), gainMargindB(k), Mpeak(k), wPeak(k), Tsettle(k), Sig(k));
end
fclose(fid);

poles = {pole(W0), pole(F0), pole(Fd), pole(Fcorr)};
ptag = {'original_open', 'original_closed', 'desired_closed', 'corrected_closed'};
fid = fopen(fullfile(resdir, 'poles.csv'), 'w');
fprintf(fid, 'system,index,re,im\n');
for k = 1:numel(poles)
    pk = poles{k}(:);
    [~, idx] = sortrows([-real(pk), -abs(imag(pk))]);
    pk = pk(idx);
    for j = 1:numel(pk)
        re = real(pk(j)); im = imag(pk(j));
        if abs(re) < 1e-9, re = 0; end
        if abs(im) < 1e-9, im = 0; end
        fprintf(fid, '%s,%d,%.6f,%.6f\n', ptag{k}, j, re, im);
    end
end
fclose(fid);

fid = fopen(fullfile(resdir, 'corrector.csv'), 'w');
fprintf(fid, 'section,zero_rad_s,pole_rad_s,T_zero,T_pole,alpha,K_hf\n');
fprintf(fid, '1,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', w0, b, Tzero, Tpole, alpha, Khf);
fclose(fid);

wBreak = [0.1 w0 b 1e4]; sa = {[-20 -40 -40], [-20 -20 -40], [0 20 0]};
La = zeros(3, numel(wBreak));
La(:, 1) = [20*log10(K/wBreak(1)); 20*log10(K/wBreak(1)); 0];
for k = 1:3
    for j = 1:3
        La(k, j+1) = La(k, j) + sa{k}(j)*log10(wBreak(j+1)/wBreak(j));
    end
end
wGrid = logspace(log10(wBreak(1)), log10(wBreak(end)), 1200);
ampInterp = zeros(3, numel(wGrid));
for k = 1:3
    ampInterp(k, :) = interp1(log10(wBreak), La(k, :), log10(wGrid));
end
wcAsym = nan(1, 2);
for k = 1:2
    j = find(La(k, 1:end-1) > 0 & La(k, 2:end) <= 0, 1);
    wcAsym(k) = wBreak(j)*10^(La(k, j)/abs(sa{k}(j)));
end

wVec = logspace(-1, 3, 1500).'; H = cell(1, n);
for k = 1:n
    H{k} = squeeze(freqresp(openSys{k}, wVec));
end
Hk = squeeze(freqresp(Wk, wVec)); Hc = squeeze(freqresp(Fcorr, wVec));

t0 = linspace(0, 1.5*Tsettle(1), 4000).'; y0 = step(F0, t0);
[ym1, i1] = max(y0); t1 = t0(i1);
f = canvas(19, 7);
h = plot(t0, y0, 'LineWidth', 1.4); hold on;
plot(t1, ym1, 'o', 'Color', h.Color, 'MarkerFaceColor', h.Color, 'MarkerSize', 5);
yline(1, ':', 'LineWidth', 0.9);
text(t1 + 0.25, ym1 - 0.18, sprintf('\\sigma = %.1f %%', Sig(1)), 'FontSize', 9);
grid on; box on; xlim([0, t0(end)]);
xlabel('t, с'); ylabel('y(t)');
title('Переходная характеристика исходной замкнутой системы', 'FontWeight', 'normal');
savepng(f, 'fig01_step_original.png', figdir);

f = canvas(19, 12);
subplot(2, 1, 1);
hk = cell(1, 2);
for k = 1:2
    hk{k} = semilogx(wVec, 20*log10(abs(H{k})), 'LineWidth', 1.3); hold on;
end
yline(0, ':'); xline(gainCrossFreq(1), ':', '\omega_{c1}'); xline(gainCrossFreq(2), ':', '\omega_{c2}');
grid on; box on; xlim([wVec(1), wVec(end)]);
xlabel('\omega, рад/с'); ylabel('L(\omega), дБ');
title('ЛАЧХ разомкнутых систем', 'FontWeight', 'normal');
legend([hk{:}], {sprintf('исходная, \\omega_c = %.2f рад/с', gainCrossFreq(1)), ...
    sprintf('желаемая, \\omega_c = %.2f рад/с', gainCrossFreq(2))}, ...
    'Location', 'southwest'); legend boxoff;
subplot(2, 1, 2);
for k = 1:2
    semilogx(wVec, unwrap(angle(H{k}))*180/pi, 'LineWidth', 1.3); hold on;
end
yline(-180, ':'); xline(gainCrossFreq(1), ':'); xline(gainCrossFreq(2), ':');
grid on; box on; xlim([wVec(1), wVec(end)]); ylim([-270, -60]);
xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ разомкнутых систем', 'FontWeight', 'normal');
savepng(f, 'fig02_bode_open.png', figdir);

f = canvas(19, 7);
h1 = semilogx(wGrid, ampInterp(1, :), 'LineWidth', 1.5); hold on;
h2 = semilogx(wGrid, ampInterp(2, :), 'LineWidth', 1.5);
for j = 2:3
    xline(wBreak(j), ':', sprintf('%.4g рад/с', wBreak(j)));
end
grid on; box on; xlim([wBreak(1), wBreak(end)]);
xlabel('\omega, рад/с'); ylabel('L(\omega), дБ');
title('Асимптотические ЛАЧХ исходной и желаемой систем', 'FontWeight', 'normal');
legend([h1 h2], {'L0 исходной', 'Lж желаемой'}, 'Location', 'southwest'); legend boxoff;
savepng(f, 'fig03_lach_asymptotic.png', figdir);

f = canvas(19, 7);
h1 = semilogx(wGrid, ampInterp(3, :), 'LineWidth', 1.6); hold on;
h2 = semilogx(wVec, 20*log10(abs(Hk)), '--', 'LineWidth', 1.2);
for j = 2:3
    xline(wBreak(j), ':', 'LineWidth', 0.9);
end
grid on; box on; xlim([wBreak(1), wBreak(end)]);
xlabel('\omega, рад/с'); ylabel('Lку(\omega), дБ');
title('ЛАЧХ последовательного корректирующего устройства', 'FontWeight', 'normal');
legend([h1 h2], {'асимптотическая Lку', 'точная Lку'}, 'Location', 'southeast');
legend boxoff;
savepng(f, 'fig04_lku.png', figdir);

tc = linspace(0, 1.6*Tsettle(3), 4000).'; yc = step(Fcorr, tc);
[ym3, i3] = max(yc); t3 = tc(i3);
f = canvas(19, 7);
h = plot(tc, yc, 'LineWidth', 1.5); hold on;
plot(t3, ym3, 'o', 'Color', h.Color, 'MarkerFaceColor', h.Color, 'MarkerSize', 5);
yline(1.02, '--'); yline(0.98, '--');
xline(Tsettle(3), ':', 'HandleVisibility', 'off');
text(Tsettle(3) + 0.004, 0.25, sprintf('T_{пп} = %.4g с', Tsettle(3)), 'FontSize', 9);
text(t3, ym3 + 0.05, sprintf('максимум %.3f', ym3), 'FontSize', 9);
grid on; box on; xlim([0, tc(end)]); ylim([0, 1.45]);
xlabel('t, с'); ylabel('y(t)');
title('Переходная характеристика скорректированной системы', 'FontWeight', 'normal');
savepng(f, 'fig05_step_corrected.png', figdir);

f = canvas(19, 12);
subplot(2, 1, 1);
semilogx(wVec, 20*log10(abs(Hc)), 'LineWidth', 1.3); hold on;
yline(20*log10(Mpeak(3)), ':', sprintf('M = %.3g', Mpeak(3)));
xline(wPeak(3), ':', sprintf('\\omega_r = %.3g рад/с', wPeak(3)));
grid on; box on; xlim([wVec(1), wVec(end)]);
xlabel('\omega, рад/с'); ylabel('L(\omega), дБ');
title('ЛАЧХ скорректированной замкнутой системы', 'FontWeight', 'normal');
subplot(2, 1, 2);
semilogx(wVec, unwrap(angle(Hc))*180/pi, 'LineWidth', 1.3);
grid on; box on; xlim([wVec(1), wVec(end)]);
xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ скорректированной замкнутой системы', 'FontWeight', 'normal');
savepng(f, 'fig06_bode_closed.png', figdir);

function f = canvas(w, h)
f = figure('Visible', 'off', 'Color', 'white', 'Units', 'centimeters', ...
    'Position', [1 1 w h]);
set(f, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
end

function savepng(f, name, dir)
exportgraphics(f, fullfile(dir, name), 'Resolution', 300);
close(f);
end
