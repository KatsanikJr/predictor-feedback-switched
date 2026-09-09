clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Example V.B
%
% Reproduces the numerical example in Section V.B of the paper:
%
% "Predictor-Feedback Stabilization of Linear Switched Systems with
%  State-Dependent Switching or Delayed Switching Input"
%
% This example illustrates Theorem 4.3 for the switched system with a
% delayed switching input.
%
% Main paper references:
%   Plant dynamics:             Eq. (42)
%   Nominal switching law:      Eq. (47)
%   Predictor-feedback law:     Eq. (49)
%   Predictor state:            Eqs. (50)-(52)
%   Numerical plant data:       Eq. (69)
%   Lyapunov matrices:          Eq. (70)
%   Switching regions:          Eq. (71)
%
% Run this script FIRST. The simulation data are saved in
% "myworkspace_IV_B.mat".
%
% Then run:
%       plots_IV_B.m
% to generate the corresponding state/switching-input plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- System Setup: Eqs. (69)-(70) in the paper ---
A{1} = [1 -5; 0 1];
A{2} = [1  0; 5 1];

% Lyapunov matrices used in the paper, Eq. (70).
% These matrices are generated/verified from find_LMI_IV_B.m 
S1 = [ 3.2487   6.9864; 6.9864  31.3624];
S2 = [ 46.1415  -12.6366; -12.6366    7.1167];

%% Simulation parameters
D  = 2;
dt = 0.0001;
T  = 30;

time       = 0:dt:T;
numSteps   = length(time)-1;
delaySteps = round(D/dt);

% State and switching-input arrays
X = zeros(2,numSteps+1);
U = zeros(1,delaySteps+numSteps+1);

% Initial state
X(:,1) = [2; -1];

% Initial switching-input history on [-D,0):
% mode 1 on the first half and mode 1 on the second half.
halfDelay = round(delaySteps/2);
U(1:halfDelay) = 1;
U(halfDelay+1:delaySteps) = 1;


%% --- Main simulation loop ---
for i = 1:numSteps

    % Euler update of the plant (42)
    X(:,i+1) = X(:,i) + dt*(A{U(i)}*X(:,i));

    % Construct the explicit predictor in Eq. (52) using the stored
    % switching-input history over the interval [t-D,t).
    % Stored switching-input history over one delay interval

    % Current delay window 
    window = U(i:i+delaySteps-1);

    % Switching instants in the interval [t-D,t)
    switches_D   = find(diff(window) ~= 0) + i;
    num_switches = length(switches_D);

   % Store values before each switch and the last value     
    mn = U(i);
    sn = i;

    if num_switches > 0
        for n = 1:num_switches
            mn = [mn, U(switches_D(n))];
            sn = [sn, switches_D(n)];
        end
    end

    sn = [sn, i+delaySteps];

    % Explicit predictor in Eq. (52)
    prodX = eye(2);

    for n = 1:num_switches+1
        flow = expm(A{mn(n)}*(sn(n+1)-sn(n))*dt);
        prodX = flow*prodX;
    end

    Pt = prodX*X(:,i);

    % Predictor-feedback switching law, Eqs. (47), (49)
    U(i+delaySteps) = sigma_paper(Pt,S1,S2);

end

% --- Save Simulation Data ---
% The plotting script plots_IV_B.m loads this workspace.
% Therefore, run this simulation script before running plots_IV_B.m.save myworkspace_IV_B.mat
save myworkspace_IV_B.mat


function mode_next = sigma_paper(x,S1,S2)
% State-dependent switching law used in Eq. (47).
    %
    % E_i(x) = x' S_i x, as defined in Eq. (46).
    % For the two-mode case:
    %   mode 1 is selected if E1 > E2;
    %   otherwise mode 2 is selected.
    %
    % Thus equality is assigned to mode 2, aligned with the
    % max-index tie-breaking rule in Eq. (47).

    E1 = x.'*S1*x;
    E2 = x.'*S2*x;

    if E1 > E2
        mode_next = 1;
    else
        mode_next = 2;   
    end
end
