function options = BT_Hom_set_options()

%% Function used to set options needed for init_BT_Hom
% The option structure contains the following fields
% - `ntst` Number of mesh intervals with a default value of 40.
% - `ncol` Number of collocation points used on each interval with a default of 4.
% - `extravec` Three dimensional boolean row vector indicating which _homoclinic
%   parameters_ are selected to be free. The first component refers to the
%   half-return time, while the second and third refer to the distances from the
%   saddle point to the first, respectively, the last point on the homoclinic
%   orbit. The default value is set to `[0 1 1]`. Thus the half-return time `T`
%   is fixed.
% - `order` The order of the homoclinic approximation used with a default value of 3
% - `amplitude` Desired amplitude of the homoclinic solution. If left empty then
%   an conservative estimate is made.
% - `TTolerance` Desired distance between the last point on the homoclinic
%   solution and the saddle point. Thus should at least be smaller than the
%   amplitude. If left empty it is defined by `amplitude*1.0e-03`;
% - `HigherOrderTimeReparametrization` Boolean to indicate if a higher order
%   approximation to the nonlinear time transformation in the
%   Lindstedt-Poincaré method should used. This should always set to `1`.
%   It is only implemented for demonstration purposes.
% - `method` Selects the method to be used to approximate the homoclinic
%   solution. The different methods available are:
%   - orbital (the default),
%   - orbitalv2,
%   - LP (Lindstedt-Poincaré),
%   - LPHypernormalForm,
%   - RegularPerturbation,
%   - RegularPerturbationL2,
%   We refer to {cite}`Bosschaert@2021` for the interpretations.
% - `messages` Boolean to indicate if information about selected parameter should
%   be printed the console. The default value is set to `true`.
% - `correct` Boolean to indicate if the homoclinic solutions should be correct
%   with Newton. The default value is set to `true`.

options.ntst = 40;
options.ncol = 4;
options.extravec = [0 1 1];
options.order = 3;
options.TTolerance = [];
options.amplitude = [];
options.HigherOrderTimeReparametrization = 1;
options.method = 'orbital';
options.messages = true;
options.perturbationparameter = 0.1;
options.correct = true;
options.MaxNumCorrections = 10;
options.amplitudeTToleranceRatio = 1.0e-03;

end
