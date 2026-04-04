% adaptacio.m
% Impedance matching network design using Smith chart (Gamma plane) operations
%
% Network 1 — Shunt susceptance + transmission line:
%   Port 1 ---+--- [line Z0, ell] --- ZL
%             |
%            jB
%            GND
%
% Network 2 — Transmission line + series reactance:
%   Port 1 --- [line Z0, ell] --- [+jX] --- ZL
%
% Network 3 — Transmission line + open-circuit stub (shunt):
%   Port 1 ---+--- [line Z0, ell1] --- ZL
%             |
%        [OC stub, ell_s]
%            GND
%
% Network 4 — Transmission line + short-circuit stub (shunt):
%   Port 1 ---+--- [line Z0, ell1] --- ZL
%             |
%        [SC stub, ell_s]
%            GND
%
% Network 5 — Lambda/4 transformer (Z0') + transmission line (Z0):
%   Port 1 --- [line Z0', lambda/4] --- [line Z0, ell1] --- ZL
%
% Network 6 — Lambda/4 transformer (Z0') + series reactance:
%   Port 1 --- [line Z0', lambda/4] --- [+jX] --- ZL
%
% Network 7 — Lambda/4 transformer (Z0') + shunt susceptance:
%   Port 1 ---+--- [line Z0', lambda/4] ---+--- ZL
%             |                             |
%            (Z0)                          jB
%                                          GND

Z0 = 50;   % system characteristic impedance [ohm]

fprintf('=== Impedance matching network design (Z0 = %.0f ohm) ===\n\n', Z0);

R   = input('Normalized load resistance  R = Re(ZL/Z0): ');
X   = input('Normalized load reactance   X = Im(ZL/Z0): ');

fprintf('\nAvailable networks:\n');
fprintf('  [1] Shunt jB + transmission line\n');
fprintf('  [2] Transmission line + series jX\n');
fprintf('  [3] Transmission line + open-circuit stub (shunt)\n');
fprintf('  [4] Transmission line + short-circuit stub (shunt)\n');
fprintf('  [5] Lambda/4 transformer + transmission line\n');
fprintf('  [6] Lambda/4 transformer + series jX\n');
fprintf('  [7] Lambda/4 transformer + shunt jB\n');

net = input('Network type (1-7): ');

if ~ismember(net, 1:7)
    fprintf('Invalid network type. Choose 1 to 7.\n');
    return;
end

zL = R + 1j*X;
yL = 1 / zL;

% -------------------------------------------------------------------------
% Load reflection coefficient
% -------------------------------------------------------------------------
GammaL = (zL - 1) / (zL + 1);
rho    = abs(GammaL);
phi    = angle(GammaL);   % rad, in (-pi, pi]

fprintf('\nNormalized load : zL = %.4f %+.4fj\n', real(zL), imag(zL));
fprintf('Normalized admit: yL = %.4f %+.4fj\n', real(yL), imag(yL));
fprintf('GammaL          : |GammaL| = %.6f,  angle = %+.2f deg\n\n', ...
        rho, rad2deg(phi));

if rho < 1e-9
    fprintf('Load is already matched to Z0. No network needed.\n');
    return;
end

% =========================================================================
switch net

% =========================================================================
% NETWORK 1: Shunt jB + transmission line
% =========================================================================
case 1
    fprintf('--- Network 1: Shunt susceptance jB + Transmission line ---\n\n');

    % Rotate GammaL until it hits G=1 circle: cos(psi) = -|GammaL|
    cos_psi = -rho;
    psi_base = acos(cos_psi);
    psi_sols = [psi_base, -psi_base];

    for k = 1:2
        psi = psi_sols(k);
        ell_lam = mod((phi - psi) / (4*pi), 0.5);
        Gamma_in = rho * exp(1j * psi);

        Yin_line  = (1 - Gamma_in) / (1 + Gamma_in);
        B_norm    = -imag(Yin_line);
        B_real    =  B_norm / Z0;

        % Verification
        Yin_total = Yin_line + 1j*B_norm;
        Gamma_total = (1 - Yin_total) / (1 + Yin_total);

        fprintf('Solution %d:\n', k);
        fprintf('  ell/lambda            = %.6f   (%.2f deg)\n', ell_lam, ell_lam*360);
        fprintf('  ell (real, Z0=%d ohm) = %.6f * lambda\n', Z0, ell_lam);
        fprintf('  B (normalized, *Y0)   = %+.6f\n', B_norm);
        fprintf('  B (real)              = %+.6f S   (%+.4f mS)\n', B_real, B_real*1e3);
        fprintf('  Gamma_total           = %.6f %+.6fj  [target: 0+0j]\n', ...
                real(Gamma_total), imag(Gamma_total));
        fprintf('\n');
    end

% =========================================================================
% NETWORK 2: Transmission line + series jX
% =========================================================================
case 2
    fprintf('--- Network 2: Transmission line + Series reactance jX ---\n\n');

    % Rotate GammaL until it hits R=1 circle: cos(psi) = +|GammaL|
    cos_psi = +rho;
    psi_base = acos(cos_psi);
    psi_sols = [psi_base, -psi_base];

    for k = 1:2
        psi = psi_sols(k);
        ell_lam = mod((phi - psi) / (4*pi), 0.5);
        Gamma_in = rho * exp(1j * psi);

        Zin_line  = (1 + Gamma_in) / (1 - Gamma_in);
        X_norm    = -imag(Zin_line);
        X_real    =  X_norm * Z0;

        % Verification
        Zin_total = Zin_line + 1j*X_norm;
        Gamma_total = (Zin_total - 1) / (Zin_total + 1);

        fprintf('Solution %d:\n', k);
        fprintf('  ell/lambda            = %.6f   (%.2f deg)\n', ell_lam, ell_lam*360);
        fprintf('  ell (real, Z0=%d ohm) = %.6f * lambda\n', Z0, ell_lam);
        fprintf('  X (normalized, *Z0)   = %+.6f\n', X_norm);
        fprintf('  X (real)              = %+.6f ohm\n', X_real);
        fprintf('  Gamma_total           = %.6f %+.6fj  [target: 0+0j]\n', ...
                real(Gamma_total), imag(Gamma_total));
        fprintf('\n');
    end

% =========================================================================
% NETWORK 3: Transmission line + open-circuit stub (shunt)
% =========================================================================
case 3
    fprintf('--- Network 3: Line + Open-circuit stub (shunt) ---\n\n');

    % Same rotation as Network 1: find ell1 and B, then convert B to OC stub
    cos_psi = -rho;
    psi_base = acos(cos_psi);
    psi_sols = [psi_base, -psi_base];

    for k = 1:2
        psi = psi_sols(k);
        ell1_lam = mod((phi - psi) / (4*pi), 0.5);
        Gamma_in = rho * exp(1j * psi);

        Yin_line = (1 - Gamma_in) / (1 + Gamma_in);
        B_norm   = -imag(Yin_line);

        % Open-circuit stub: Yin_stub = -j*cot(beta*ell) => B = -cot(beta*ell)
        % ell_stub/lambda = mod(atan(-1/B) / (2*pi), 0.5)
        ell_stub_lam = mod(atan2(-1, B_norm) / (2*pi), 0.5);

        % Verification: stub susceptance
        B_stub = -1 / tan(2*pi * ell_stub_lam);
        Yin_total = Yin_line + 1j*B_stub;
        Gamma_total = (1 - Yin_total) / (1 + Yin_total);

        fprintf('Solution %d:\n', k);
        fprintf('  ell1/lambda (series line)  = %.6f   (%.2f deg)\n', ell1_lam, ell1_lam*360);
        fprintf('  B required (normalized)    = %+.6f\n', B_norm);
        fprintf('  ell_stub/lambda (OC stub)  = %.6f   (%.2f deg)\n', ell_stub_lam, ell_stub_lam*360);
        fprintf('  B_stub verification        = %+.6f  [target: %+.6f]\n', B_stub, B_norm);
        fprintf('  Gamma_total                = %.6f %+.6fj  [target: 0+0j]\n', ...
                real(Gamma_total), imag(Gamma_total));
        fprintf('\n');
    end

% =========================================================================
% NETWORK 4: Transmission line + short-circuit stub (shunt)
% =========================================================================
case 4
    fprintf('--- Network 4: Line + Short-circuit stub (shunt) ---\n\n');

    % Same rotation as Network 1: find ell1 and B, then convert B to SC stub
    cos_psi = -rho;
    psi_base = acos(cos_psi);
    psi_sols = [psi_base, -psi_base];

    for k = 1:2
        psi = psi_sols(k);
        ell1_lam = mod((phi - psi) / (4*pi), 0.5);
        Gamma_in = rho * exp(1j * psi);

        Yin_line = (1 - Gamma_in) / (1 + Gamma_in);
        B_norm   = -imag(Yin_line);

        % Short-circuit stub: Yin_stub = j*tan(beta*ell) => B = tan(beta*ell)
        % ell_stub/lambda = mod(atan(B) / (2*pi), 0.5)
        ell_stub_lam = mod(atan2(B_norm, 1) / (2*pi), 0.5);

        % Verification: stub susceptance
        B_stub = tan(2*pi * ell_stub_lam);
        Yin_total = Yin_line + 1j*B_stub;
        Gamma_total = (1 - Yin_total) / (1 + Yin_total);

        fprintf('Solution %d:\n', k);
        fprintf('  ell1/lambda (series line)  = %.6f   (%.2f deg)\n', ell1_lam, ell1_lam*360);
        fprintf('  B required (normalized)    = %+.6f\n', B_norm);
        fprintf('  ell_stub/lambda (SC stub)  = %.6f   (%.2f deg)\n', ell_stub_lam, ell_stub_lam*360);
        fprintf('  B_stub verification        = %+.6f  [target: %+.6f]\n', B_stub, B_norm);
        fprintf('  Gamma_total                = %.6f %+.6fj  [target: 0+0j]\n', ...
                real(Gamma_total), imag(Gamma_total));
        fprintf('\n');
    end

% =========================================================================
% NETWORK 5: Lambda/4 transformer (Z0') + transmission line (Z0)
% =========================================================================
case 5
    fprintf('--- Network 5: Lambda/4 transformer + Transmission line ---\n\n');

    % Rotate GammaL until Im(Zin) = 0 (real axis of Smith chart)
    % psi such that Zin is real => psi = -phi or psi = pi - phi
    psi_sols = [-phi, pi - phi];

    for k = 1:2
        psi = psi_sols(k);
        ell1_lam = mod((phi - psi) / (4*pi), 0.5);
        Gamma_in = rho * exp(1j * psi);

        Zin_at_junction = (1 + Gamma_in) / (1 - Gamma_in);
        Rin = real(Zin_at_junction);

        if Rin <= 0
            fprintf('Solution %d: invalid (Rin = %.4f <= 0)\n\n', k, Rin);
            continue;
        end

        % Lambda/4 transformer: Z0' = sqrt(Z0 * Zin_real) normalized => sqrt(Rin)
        Z0p_norm = sqrt(Rin);
        Z0p_real = Z0p_norm * Z0;

        % Verification: lambda/4 transforms Rin to Z0'^2/Rin
        Zin_after_quarter = Z0p_norm^2 / Rin;  % should be 1
        Gamma_total = (Zin_after_quarter - 1) / (Zin_after_quarter + 1);

        fprintf('Solution %d:\n', k);
        fprintf('  ell1/lambda (Z0 line)      = %.6f   (%.2f deg)\n', ell1_lam, ell1_lam*360);
        fprintf('  Zin at junction (norm)     = %.6f %+.6fj  [should be real]\n', ...
                real(Zin_at_junction), imag(Zin_at_junction));
        fprintf('  Z0'' (normalized)           = %.6f\n', Z0p_norm);
        fprintf('  Z0'' (real)                 = %.4f ohm\n', Z0p_real);
        fprintf('  Lambda/4 length            = 0.250000 * lambda   (90.00 deg)\n');
        fprintf('  Gamma_total                = %.6f %+.6fj  [target: 0+0j]\n', ...
                real(Gamma_total), imag(Gamma_total));
        fprintf('\n');
    end

% =========================================================================
% NETWORK 6: Lambda/4 transformer (Z0') + series jX
% =========================================================================
case 6
    fprintf('--- Network 6: Lambda/4 transformer + Series reactance jX ---\n\n');

    % Add series jX to make ZL real: X_added = -Im(zL)
    X_norm = -imag(zL);
    zL_real = real(zL) + 1j*(imag(zL) + X_norm);  % = real(zL) + 0j
    Rin = real(zL_real);

    if Rin <= 0
        fprintf('Cannot match: Re(ZL) = %.4f <= 0\n', Rin);
        return;
    end

    % Lambda/4 transformer: Z0' = sqrt(Rin)
    Z0p_norm = sqrt(Rin);
    Z0p_real = Z0p_norm * Z0;
    X_real   = X_norm * Z0;

    % Verification
    Zin_after_quarter = Z0p_norm^2 / Rin;
    Gamma_total = (Zin_after_quarter - 1) / (Zin_after_quarter + 1);

    fprintf('Solution (unique):\n');
    fprintf('  Series X (normalized)      = %+.6f\n', X_norm);
    fprintf('  Series X (real)            = %+.4f ohm\n', X_real);
    fprintf('  ZL after series X (norm)   = %.6f %+.6fj  [should be real]\n', ...
            real(zL_real), imag(zL_real));
    fprintf('  Z0'' (normalized)           = %.6f\n', Z0p_norm);
    fprintf('  Z0'' (real)                 = %.4f ohm\n', Z0p_real);
    fprintf('  Lambda/4 length            = 0.250000 * lambda   (90.00 deg)\n');
    fprintf('  Gamma_total                = %.6f %+.6fj  [target: 0+0j]\n', ...
            real(Gamma_total), imag(Gamma_total));
    fprintf('\n');

% =========================================================================
% NETWORK 7: Lambda/4 transformer (Z0') + shunt jB
% =========================================================================
case 7
    fprintf('--- Network 7: Lambda/4 transformer + Shunt susceptance jB ---\n\n');

    % Add shunt jB to make YL real: B_added = -Im(yL)
    B_norm = -imag(yL);
    yL_real = real(yL) + 1j*(imag(yL) + B_norm);  % = real(yL) + 0j
    Gin = real(yL_real);

    if Gin <= 0
        fprintf('Cannot match: Re(YL) = %.4f <= 0\n', Gin);
        return;
    end

    % Lambda/4 transformer: Z0' = 1/sqrt(Gin) (normalized)
    Z0p_norm = 1 / sqrt(Gin);
    Z0p_real = Z0p_norm * Z0;
    B_real   = B_norm / Z0;

    % Verification: lambda/4 with Z0' transforms Y_real to Z0'^2 * Y_real
    % In normalized terms: Zin = Z0'^2 * Gin = (1/Gin) * Gin = 1
    Zin_after_quarter = Z0p_norm^2 * Gin;
    Gamma_total = (Zin_after_quarter - 1) / (Zin_after_quarter + 1);

    fprintf('Solution (unique):\n');
    fprintf('  Shunt B (normalized, *Y0)  = %+.6f\n', B_norm);
    fprintf('  Shunt B (real)             = %+.6f S   (%+.4f mS)\n', B_real, B_real*1e3);
    fprintf('  YL after shunt B (norm)    = %.6f %+.6fj  [should be real]\n', ...
            real(yL_real), imag(yL_real));
    fprintf('  Z0'' (normalized)           = %.6f\n', Z0p_norm);
    fprintf('  Z0'' (real)                 = %.4f ohm\n', Z0p_real);
    fprintf('  Lambda/4 length            = 0.250000 * lambda   (90.00 deg)\n');
    fprintf('  Gamma_total                = %.6f %+.6fj  [target: 0+0j]\n', ...
            real(Gamma_total), imag(Gamma_total));
    fprintf('\n');

end
