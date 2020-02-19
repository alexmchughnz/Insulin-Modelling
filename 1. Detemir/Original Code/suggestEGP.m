function [ EGP, EGPArray ] = suggestEGP( sys, tol, maxIter )
%SUGGESTEGP Suggests an EGP for a pt by setting t very large and G_bolus very small, solving for equilibrium
%   Detailed explanation goes here
    
    if nargin == 1
        tol = 0.01; % This is the acceptable difference between G(1000) and G(Fast)
        maxIter = 10; % Unless really, really, really out this seems to work well
    end
    
    % Create new ptStr with big time, and essentially no G_bolus
    sys2 = sys;
    sys2.max_t = sys.max_t;
    sys2.GI.D = 0.00001;
    sys2.SC.Ibolus = 0.00001;
    
%     % Turns off simulation-time warnings and runs forward simulation
%     warning('off', 'pt_fs:insulin_overdose_warning')
%     warning('off', 'pt_fs:hypo_warning')
    sys2 = eq_solve(sys2);
    
    % Gets the current error - difference between 'fasting' and end of
    % simulation glucoses
    currErr = sys2.GC.Gfast - sys2.results.G(end);
    i = 1;
    EGPArray = zeros(maxIter, 1);
    while (abs(currErr) > tol && i <= maxIter)
        % Estimates the change in EGP required to get to G_fast, makes the
        % change and resimulate.
        d_EGP = (sys2.GC.SI * sys2.GC.Gfast * sys2.results.Q(end) ...
            - sys2.GC.SI * sys2.results.G(end) *  sys2.results.Q(end) + ...
            sys2.GC.pg * currErr) * sys2.GC.VG;
        
        sys2.GC.EGP = sys2.GC.EGP + d_EGP;
        EGPArray(i) = sys2.GC.EGP;
        sys2 = eq_solve(sys2);
        currErr = sys2.GC.Gfast - sys2.results.G(end);
        i = i+1;
    end 
    
    if i == maxIter
        % Rough convergence check
        warning('Fn sugggestEGP did not converge sufficiently!')
    end
    
%     % Check to see if simulated final insulin levels are different from
%     % initial values
%     if abs(sys2.results.Q(end) - sys2.results.I_int_0) > 0.5 || ...
%             abs(sys2.D.I_plasma(end) - sys2.I.I_plas_0) > 0.5
%         warstr = ['Solving EGP shows Int and/or Plas Insulin bad', newline, ...
%             'Suggest I_int_0 = ', num2str(sys2.D.I_int(end)), ' and I_plas_0 = ', num2str(sys.D.I_plasma(end)), newline,...
%             'Based on ', num2str(i-1), ' iterations'];
%         warning(warstr)
%     end
    
    EGP = sys2.GC.EGP;
    
%     % Turn warnings back on
%     warning('on', 'pt_fs:insulin_overdose_warning')
%     warning('on', 'pt_fs:hypo_warning')
end

