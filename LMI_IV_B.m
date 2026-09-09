clc;
clear all;
close all;
cvx_clear;

% Grid-based BMI/LMI search for the numerical example in Section V.B.
%
% The script searches for matrices S1, S2, Q1 and Q2 satisfying
% Assumption 4.1. The search verifies the regional Lyapunov decrease
% condition in Eq. (43) with the no-sliding conditions in Eqs. (44)-(45).
%
% The S-procedure multipliers g12 and g21 are scanned over a logarithmic
% grid. At each fixed pair (g12,g21), the problem is solved as an
% LMI feasibility problem using CVX.
%
% The feasible matrices obtained are the ones reported in Eq. (70) and used in Sim_IV_B.m.
%
% CVX toolbox for MATLAB is required to run this script.

%% System matrices from Eq. (69)
A{1} = [1 -5; 0 1];
A{2} = [1  0; 5 1];

n = 2;

% Small positive bounds used to ensure strict feasibility
epsS = 1e-6;
epsL = 1e-7;

% Grid used for the S-procedure multipliers
gammas = logspace(-4,4,41);

found = false;

% Numerical bounds used in the no-sliding feasibility checks
LAMBDA_MAX = 1e4;
KAPPA_MAX  = 1e4;

for g12 = gammas
    for g21 = gammas

        % S-procedure LMI used to enforce the regional Lyapunov
        % decrease condition in Eq. (43).
        %
        % For mode i, Eq. (43) requires
        %
        %   x'*(Ai'*Si + Si*Ai + Qi)*x <= 0, for x'*(Si-Sj)*x >= 0.
        % A sufficient S-procedure condition is constructed below using
        % the nonnegative multiplier gamma_ij.

        cvx_begin sdp quiet
            cvx_precision high

            variable S1(n,n) symmetric
            variable S2(n,n) symmetric
            variable Q1(n,n) symmetric
            variable Q2(n,n) symmetric

            S1 >= epsS*eye(n);
            S2 >= epsS*eye(n);
            Q1 >= epsS*eye(n);
            Q2 >= epsS*eye(n);
            
            % Mode 1: sufficient LMI for Eq. (43) in Omega_1
            A{1}'*S1 + S1*A{1} + Q1 + g12*(S1-S2) <= -epsL*eye(n);
            
            % Mode 2: sufficient LMI for Eq. (43) in Omega_2
            A{2}'*S2 + S2*A{2} + Q2 + g21*(S2-S1) <= -epsL*eye(n);

            minimize(0)
        cvx_end

        if ~contains(cvx_status,'Solved')
            continue;
        end

        % Store the numerical solution returned by CVX
        S1_sol = full(S1);
        S2_sol = full(S2);
        Q1_sol = full(Q1);
        Q2_sol = full(Q2);

        % Check the LMI residuals numerically before proceeding
        M1 = A{1}'*S1_sol + S1_sol*A{1} + Q1_sol + g12*(S1_sol-S2_sol);

        M2 = A{2}'*S2_sol + S2_sol*A{2} + Q2_sol + g21*(S2_sol-S1_sol);

        e1 = eig((M1+M1')/2);
        e2 = eig((M2+M2')/2);

        if ~(max(e1) <= -0.5*epsL &&  max(e2) <= -0.5*epsL)
            continue;
        end

        % Check the no-sliding conditions in Eqs. (44)-(45)
        % for both pairs (1,2) and (2,1).
        I = eye(n);

        ok12 = check_sliding_pair(A{1},A{2},S1_sol-S2_sol,I,LAMBDA_MAX,KAPPA_MAX);

        ok21 = check_sliding_pair(A{2},A{1},S2_sol-S1_sol,I,LAMBDA_MAX,KAPPA_MAX);

        % Feasible solution satisfying the LMI tests associated with
        % Eqs. (43)-(45). These matrices are the ones reported in Eq. (70)
        % and used in the simulation script Sim_IV_B.m.
        if ok12 && ok21
            fprintf('\nFeasible solution found.\n');
            fprintf('g12 = %g, g21 = %g\n',g12,g21);

            disp('S1 ='); disp(S1_sol);
            disp('S2 ='); disp(S2_sol);
            disp('Q1 ='); disp(Q1_sol);
            disp('Q2 ='); disp(Q2_sol);

            found = true;
            break;
        end

    end

    if found
        break;
    end
end

if ~found
    disp('No feasible solution found on the selected gamma grid.');
end


function ok = check_sliding_pair(Ai,Aj,S,I,LAMBDA_MAX,KAPPA_MAX)

% LMI check associated with the no-sliding conditions (44)-(45).
%
% Here S = Si-Sj and the equality set is
%
%       Gamma_ij = {x : x'*S*x = 0}.
%
% The first feasibility problem provides a sufficient condition for
% Eq. (44), and the second provides a sufficient condition for
% Eq. (45).
%
% The function returns true if either alternative is feasible.

    ok = false;

    % First alternative: sufficient LMI for Eq. (44)
    cvx_begin sdp quiet
        cvx_precision high

        variables lam kap

        lam >= 0;
        lam <= LAMBDA_MAX;
        abs(kap) <= KAPPA_MAX;

        M = (Ai' - lam*Aj' + kap*I)*S + S*(Ai - lam*Aj + kap*I);

        M >= 0;
    cvx_end

    if contains(cvx_status,'Solved')
        ok = true;
        return;
    end

    % Second alternative: sufficient LMI for Eq. (45)
    cvx_begin sdp quiet
        cvx_precision high

        variables lam kap

        lam >= 0;
        lam <= LAMBDA_MAX;
        abs(kap) <= KAPPA_MAX;

        M = (Ai' - lam*Aj' + kap*I)*S + S*(Ai - lam*Aj + kap*I);

        M <= 0;
    cvx_end

    if contains(cvx_status,'Solved')
        ok = true;
    end
end