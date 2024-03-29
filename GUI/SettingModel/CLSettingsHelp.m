classdef CLSettingsHelp
    properties(Constant)
        HELP = struct(...
            'InitStepsize', 'the initial stepsize', ...
            'MinStepsize', 'the minimum stepsize to compute the next point on the curve', ...
            'MaxStepsize', 'the maximum stepsize', ...
            'MaxNewtonIters','maximum number of Newton-Raphson iterations before switching to Newton-Chords in the corrector iterations', ...
            'MaxCorrIters', 'maximum number of correction iterations', ...
            'MaxTestIters', 'maximum number of iterations to locate a zero of a test function', ...
            'VarTolerance', 'tolerance of coordinates: ||δx|| <= VarTolerance  is the second convergencecriterium of the Newton iteration', ...
            'FunTolerance', 'tolerance of function values: ||F(x)|| <= FunTolerance  is the first convergencecriterium of the Newton iteration', ...
            'TestTolerance', 'tolerance of test functions', ...
            'Adapt', 'number of points indicating when to adapt the problem while computing the curve', ...
            'MaxNumPoints', 'maximum number of points on the curve', ...
            'CheckClosed',   'number of points indicating when to start to check if the curve is closed (0=do not check)', ...
            'eigenvalues', 'enable the computation of the eigenvalues', ...
            'multipliers', 'enable the computation of the multipliers', ...
            'option_tsearchorder', 'this value indicates if unit vectors are cycled in increasing order of index (true) or decreasing (false)', ...
            'Period', 'period of limit cycle', ...
            'ntst', 'number of test intervals', ...
            'ncol', 'number of collocation points', ...
            'bt_amplitude', 'Initial amplitude of the homoclinic orbit.  If left empty the amplitude will automatically be determined.', ...
            'bt_ttolerance', 'Distance between saddle point and homoclinic orbit. If left empty the distance will automatically be determined.', ...
            'whichNS', 'select a NS-curve',...
            'option_moorepenrose', 'enable the use of the Moore-Penrose continuation as the Newton-like corrector procedure (default: true)', ...
            'RelTolerance', 'odeset: this tolerance measures the error relative to the magnitude of each solution component',...
            'AbsTolerance', 'odeset: this tolerance is a threshold below which the value of the solution becomes unimportant', ...
            'Normcontrol', 'odeset: when NormControl is on, the solvers control the error e at each step using the norm of the solution rather than its absolute value', ...
            'Refine', 'odeset: if the refinement factor is n>1, then the solver subdivides each step into n smaller intervals and returns solutions at each point', ...
            'BDF', 'odeset: toggle to use backward differentiation formulas (BDFs)', ... 
            'MaxOrder', 'odeset: specify the maximum order used in the numerical differentiation formulas (NDFs) or backward differentiation formulas (BDFs) that are used by the variable-order solvers', ...
            'InitStepSize_sim', 'odeset: suggested initial step size, it sets an upper bound on the magnitude of the first step size that the ODE solver tries',...
            'MaxStepSize_sim', 'odeset: sets an upper bound on the size of any step taken by the solver',...
            'Interval', 'length of the interval of integration', ...
            'option_increment', 'the increment to compute first order derivatives numerically', ...
            'eventfunction', 'odeset: name of the event function (doc odeset). The name entered must be the name of a function accessible in the current Command Window' ...
            );
            %'option_pause', 'option_archive', 'option_output',
            %'amplitude', 'eps'
    end
    
    methods(Static)
        function helpstr = getHelp(keyword)
            if isfield(CLSettingsHelp.HELP, keyword)
                helpstr = CLSettingsHelp.HELP.(keyword);
            else
                helpstr = '';
            end
            
            
        end
        
        
    end
end
