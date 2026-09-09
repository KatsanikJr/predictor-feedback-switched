clc; 
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Example V.A
%
% Reproduces the numerical example in Section V.A of the paper:
%
% "Predictor-Feedback Stabilization of Linear Switched Systems with
%  State-Dependent Switching or Delayed Switching Input"
%
% This example illustrates Theorem 3.1 for the state-dependent switching
% system with a constant delay in the input.
%
% Main paper references:
%   Plant dynamics:          Eq. (1)
%   Predictor state:         Eq. (8)
%   Predictor-feedback law:  Eq. (9)
%   Numerical plant data:    Eq. (64)
%   Gains/Lyapunov matrices: Eq. (65)
%   Switching rule:          Eqs. (66)-(68)
%
% Run this script FIRST. The simulation data are saved in
% "myworkspace_IV_A.mat".
%
% Then run:
%       plots_IV_A.m
% to generate the corresponding state/control plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- System Setup: Eq. (64) in the paper ---
A{1} = [1  -5;  0  3];
A{2} = [3  0;  5  1];

B{1} = [0; 1];
B{2} = [1; 0];

% Matrices S_1 and S_2 in Eq. (65) of the paper.
S_1 = [1.0432 1.0031 ; 1.0031 2.9722];
S_2 = [2.9722 -1.0031 ; -1.0031 1.0432];

% State-feedback gains K_1 and K_2 in Eq. (65)
K{1} = [0 -2];
K{2} = [-2 0];

%% --- Simulation Parameters ---

% Parameters used in Section V.A of the paper:
% D = 2 s, X(0) = [2 -1]^T, U(s) = 0 for s in [-D,0),
% and predictor/simulation step h = Delta t = 10^(-3) s.
D        = 2;
dt       = 0.001;
h_step   = dt;
T        = 10;

time     = 0:dt:T;
numSteps = length(time)-1;
delaySteps = round(D/dt);
N        = round(D/h_step);

% Preallocate
X       = zeros(2, numSteps + 1);
P       = zeros(2, delaySteps + numSteps + 1);
U       = zeros(1, delaySteps + numSteps + 1);
sigmaX  = zeros(1, numSteps + 1);
sigmaP  = zeros(1, delaySteps + numSteps + 1);

% Initial conditions
X(:,1)      = [2; -1];  %X(0) = [2 -1]^T
sigmaX(1) = sigma_paper(X(:,1), S_1, S_2);
P(:,1)      = X(:,1);
sigmaP(1)   = sigmaX(1);
U(1:delaySteps) = 0;

%% --- Main simulation loop ---
for i = 1:numSteps
    sigX = sigmaX(i); 

    % Numerical computation of the implicit predictor in Eq. (8).
    %
    % Starting from the current plant state X(t), the predictor is
    % propagated over one delay interval of length D using the
    % left-endpoint (Euler) integration rule.
    %
    % Since the switching signal depends on the predictor state, the
    % active mode is recomputed at every predictor integration step.

    % Predictor window from t-D to t
    modeP  = sigmaX(i);
    P_prev = X(:, i);
   
    
    for j = 1:N
        idx_old = i + j - 1;   
        idx_new = i + j;

        % Left-endpoint approximation of Eq. (8)
        f0      = A{modeP}*P_prev + B{modeP}*U(idx_old);
        P_next  = P_prev + h_step*f0;

        % Evaluate the state-dependent switching rule, Eq. (67)
        mode_next   = sigma_paper(P_next, S_1, S_2);

        P(:, idx_new)    = P_next;
        sigmaP(idx_new)  = mode_next;

        P_prev = P_next;
        modeP  = mode_next;
    end

    % Predictor state at the end of the prediction horizon,
    % corresponding to P(t) in the feedback law (9)
    Pt     = P(:, i + delaySteps );
    sigPt  = modeP;

    % Predictor-feedback law (9)
    U(i + delaySteps) = K{sigPt} * Pt;

    % Euler update of the plant (1)
    X(:, i+1)         = X(:, i) + dt*( A{sigX}*X(:,i) + B{sigX}*U(i) );

    % Actual state-dependent switching signal sigma(X(t))
    sigmaX(i+1)       = sigma_paper(X(:,i+1),S_1, S_2);
end

% --- Save Simulation Data ---
% The plotting script plots_IV_A.m loads this workspace.
% Therefore, run this simulation script before running plots_IV_A.m.
save myworkspace_IV_A.mat

function mode_next = sigma_paper(x, P1, P2)

% State-dependent switching law used in Eq. (67).
    %
    % E_i(x) = x' S_i x, as defined in Eq. (66).
    % For the two-mode case:
    %   mode 1 is selected if E1 > E2;
    %   otherwise mode 2 is selected.
    %
    % Thus equality is assigned to mode 2, aligned with the
    % max-index tie-breaking rule in Eq. (67).
    E1 = x.'*P1*x;
    E2 = x.'*P2*x;

    if E1 > E2
        mode_next = 1;
    else
        mode_next = 2;
    end
end

