function suc = generate_directional_derivatives_octave(equations, coordinates, ...
    parameters, system_name, outputdir, order)

  %% open file to write to and print header to file
  fileID = fopen([outputdir, system_name, '.m'],'w');
  fprintf(fileID, strcat('function out = %s\n', ...
                         'out{1} = @init;\n', ...
                         'out{2} = @fun_eval;\n', ...
                         'out{3} = @jacobian;\n', ...
                         'out{4} = @jacobianp;\n', ...
                         'out{5} = @hessians;\n', ...
                         'out{6} = @hessiansp;\n', ...
                         'out{7} = @der3;\n', ...
                         'out{8} = @der4;\n', ...
                         'out{9} = @der5;\n', ...
                         'out{10} = [];\n\n'), system_name);

  coordinates = children(coordinates);
  parameters = children(parameters);

  # print wrapper functions, these are needed since octave currently doesn't
  # support index variables in the generation of functions with
  # matlabFucntion. However, we are forced to use this by the construction
  # of MatCont.
  parameterArguments = strjoin(cellfun(@(p) char(p), parameters, ...
                                       'UniformOutput', false), ',');
  fprintf(fileID, strcat('function dydt = fun_eval(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'dydt = fun_eval_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function jac = jacobian(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'jac = jacobian_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function jacp = jacobianp(t, x, %s)\n', ...
                         'x = num2cell(x);\n', ...
                         'jacp = jacobianp_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function hess = hessians(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'hess = hessians_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function hessp = hessiansp(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'hessp = hessiansp_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function tens3 = der3(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'tens3 = der3_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function tens4 = der4(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'tens4 = der4_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);
  fprintf(fileID, strcat('function tens5 = der5(t, x, %s)\n',...
                         'x = num2cell(x);\n', ...
                         'tens5 = der5_octave(t, x{:}, %s);\n\n'), ...
          parameterArguments, parameterArguments);

  %% dy_dt
  syms t;
  fh = matlabFunction(equations, 'vars', [t, coordinates{:}, parameters{:}]);
         % convert function handle fh to string and process for printing to file
  dydtString = regexprep(func2str(fh), '^@', 'function dydt = fun_eval_octave');
  dydtString = regexprep(dydtString, '\)', ')\ndydt = ', 'once');
  fprintf(fileID, '\n%s;\n\n', dydtString);

  %% init function used for simulation
  y0 = zeros(size(coordinates));
  fprintf(fileID, strcat('function [tspan,y0,options] = init\n', ...
                         'handles = feval(%s);\n', ...
                         'y0=%s;\n', ...
                         "options = odeset('Jacobian',handles(3),'JacobianP',", ...
                         "handles(4),'Hessians',handles(5),'HessiansP',handles(6));\n", ...
                         'tspan = [0 10];\n\n'), system_name, mat2str(y0));

  %% jacobian
  A = jacobian(equations, coordinates);
  fh = matlabFunction(A, 'vars', {t, coordinates{:}, parameters{:}});
  argumentlist = regexprep(func2str(fh),'\).+',')'); % used later on
  argumentlist = regexprep(argumentlist, '@', '');
  % convert function handle fh to string and process for printing to file
  jacobianString = regexprep(func2str(fh), '^@', 'function jac = jacobian_octave');
  jacobianString = regexprep(jacobianString,'\)',')\njac = ', 'once');
  fprintf(fileID, '%s;\n\n', jacobianString);

  %% jacobianp (D0D1 f(x))
  fh = matlabFunction(jacobian(equations, parameters), 'vars', ...
                      {t, coordinates{:}, parameters{:}});
  % convert function handle fh to string and process for printing to file
  jacobianpString = regexprep(func2str(fh), '^@', 'function jacp = jacobianp_octave');
  jacobianpString = regexprep(jacobianpString, '\)', ')\njacp = ', 'once');
  fprintf(fileID, sprintf('%s;\n\n', jacobianpString));

  %% hessian
  n = length(coordinates);
  % print hessian header
  fprintf(fileID, strcat('function hess = hessians_octave', argumentlist,'\n'));
  % create empty matrix to store derivatives
  % note that GNU Octave currently doesn't support higher dimensional symbolic matrices
  d2f = sym('a', [n, n*n]);
  for i = 1:n
    d2f(1:n,(i-1)*n+1:i*n) = diff(A, coordinates(i));
    fh = matlabFunction(d2f(1:n,(i-1)*n+1:i*n), 'vars', {t, coordinates{:}, parameters{:}});
    hessiString = sprintf('hess(:,:,%d) = ', i);
    fprintf(fileID, ...
            strcat(hessiString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
  end

  %% hessianp (D1D1 f(x))
  num_pars = length(parameters);
  % print hessianp header
  fprintf(fileID, strcat('\nfunction hessp = hessiansp_octave', argumentlist,'\n'));
  for i = 1:num_pars
    fh = matlabFunction(diff(A, parameters(i)), ...
                        'vars', {t, coordinates{:}, parameters{:}});
    hesspiString = sprintf('hessp(:,:,%d) = ', i);
    fprintf(fileID, ...
            strcat(hesspiString, regexprep(func2str(fh),'@\(.*?\)',''), ';\n'));
  end

  %% tens3 (D3 f(x))
  % print hessian header
  fprintf(fileID, strcat('\nfunction tens3 = der3_octave', argumentlist,'\n'));
  % create empty matrix to store derivatives
  d3f = sym('a', [n, n*n*n]);
  for i = 1:n
    for j = 1:n
      d3f(1:n,(i+j-1)*n+1:(i+j)*n) = diff(d2f(1:n,(i-1)*n+1:i*n), coordinates(j));
      fh = matlabFunction(d3f(1:n,(i+j-1)*n+1:(i+j)*n), ...
                          'vars', {t, coordinates{:}, parameters{:}});
      tens3ijString = sprintf('tens3(:,:,%d,%d) = ', i, j);
      fprintf(fileID, ...
              strcat(tens3ijString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
    end
  end

  %% tens4 (D4 f(x))
  % print tens4 header
  fprintf(fileID, strcat('\nfunction tens4 = der4_octave', argumentlist,'\n'));
  % create empty matrix to store derivatives
  d4f = sym('a'* [n, n*n*n*n]);
  for i = 1:n
    for j = 1:n
      for k = 1:n
        d4f(1:n,(i+j+k-1)*n+1:(i+j+k)*n) = diff(d3f(1:n,(i+j-1)*n+1:(i+j)*n), coordinates(k));
        fh = matlabFunction(d4f(1:n,(i+j+k-1)*n+1:(i+j+k)*n), ...
                            'vars', {t, coordinates{:}, parameters{:}});
        tens4ijkString = sprintf('tens4(:,:,%d,%d,%d) = ', i, j, k);
        fprintf(fileID, ...
                strcat(tens4ijkString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
      end
    end
  end

  %% tens5 (D5 f(x))
  % print tens5 header
  fprintf(fileID, strcat('\nfunction tens5 = der5_octave', argumentlist,'\n'));
  for i = 1:n
    for j = 1:n
      for k = 1:n
        for l = 1:n
          d5fijkl = diff(d4f(1:n,(i+j+k-1)*n+1:(i+j+k)*n), coordinates(l));
          fh = matlabFunction(d5fijkl, 'vars', {t, coordinates{:}, parameters{:}});
          tens5ijklString = sprintf('tens5(:,:,%d,%d,%d,%d) = ', i, j, k, l);
          fprintf(fileID, ...
                  strcat(tens5ijklString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
        end
      end
    end
  end

  fclose(fileID);
  suc = 1;
end
