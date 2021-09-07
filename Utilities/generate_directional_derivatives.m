function suc = generate_directional_derivatives(equations, coordinates, ...
    parameters, system_name, outputdir, order)

  % The argument order isn't implemented here yet

  %% first create the mat file needed if the GUI is to be used
  gds.system = system_name;
  gds.equations = arrayfun(@string, equations);
  gds.coordinates = coordinates;
  gds.parameters = parameters;
  gds.dim = length(coordinates);
  gds.time = {'t', [0]};
  der = zeros(4,5);
  der(end,:) = ones(1,5);
  gds.der = der;
  save([outputdir, system_name, '.mat'],'gds');


  % take transpose of coordinates without introducing conjugates
  coordinates = reshape(coordinates, length(coordinates), 1);

  % If we are using GNU Octave call octave version of this file,
  % otherwise continue
  if isOctave
    suc = generate_directional_derivatives_octave(equations, coordinates, ...
        parameters, system_name, outputdir, order);
    return
  end

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
                   'out{10} = [];\n'), system_name);

  %% dy_dt
  syms t;
  parnames = children(parameters,1);
  fh = matlabFunction(equations, 'vars', ...
                      {t, coordinates, parnames{:}});
  % convert function handle fh to string and process for printing to file
  dydtString = regexprep(func2str(fh), ...
                         '^@', 'function dydt = fun_eval');
  dydtString = regexprep(dydtString,')',')\ndydt = ', 'once');
  fprintf(fileID, '\n%s;\n\n', dydtString);

  %% init function used for simulation
  y0 = zeros(size(coordinates));
  fprintf(fileID, strcat('function [tspan,y0,options] = init\n', ...
                    'handles = feval(%s);\n', ...
                    'y0=%s;\n', ...
                    "options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));\n", ...
                    'tspan = [0 10];\n'), system_name, mat2str(y0));

  %% jacobian
  A = jacobian(equations, coordinates);
  fh = matlabFunction(A, 'vars', {t, coordinates, parnames{:}});
  argumentlist = regexprep(func2str(fh),'\).+',')'); % used later on
  argumentlist = regexprep(argumentlist, '@', '');
  % convert function handle fh to string and process for printing to file
  jacobianString = regexprep(func2str(fh), ...
                         '^@', 'function jac = jacobian');
  jacobianString = regexprep(jacobianString,')',')\njac = ', 'once');
  fprintf(fileID, sprintf('\n%s;\n', jacobianString));

  %% jacobianp (D0D1 f(x))
  fh = matlabFunction(jacobian(equations, parameters), 'vars', ...
                      {t, coordinates, parnames{:}});
  % convert function handle fh to string and process for printing to file
  jacobianpString = regexprep(func2str(fh), ...
                             '^@', 'function jacp = jacobianp');
  jacobianpString = regexprep(jacobianpString,')',')\njacp = ', 'once');
  fprintf(fileID, sprintf('\n%s;\n\n', jacobianpString));

  %% hessian
  n = length(coordinates);
  % print hessian header
  fprintf(fileID, strcat('function hess = hessians', argumentlist,'\n'));
  % create empty matrix to store derivatives
  d2f = sym('a', [n, n, n]);
  for i = 1:n
    d2f(:,:,i) = diff(A, coordinates(i));
    fh = matlabFunction(d2f(:,:,i), 'vars', {t, coordinates, parnames{:}});
    hessiString = strcat('hess(:,:,', string(i), ') = ');
    fprintf(fileID, ...
            strcat(hessiString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
  end

  %% hessianp (D1D1 f(x))
  num_pars = length(parameters);
  % print hessian header
  fprintf(fileID, strcat('\nfunction hessp = hessiansp', argumentlist,'\n'));
  for i = 1:num_pars
    fh = matlabFunction(diff(A, parameters(i)), ...
                        'vars', {t, coordinates, parnames{:}});
    hesspiString = strcat('hessp(:,:,', string(i), ') = ');
    fprintf(fileID, ...
            strcat(hesspiString, regexprep(func2str(fh),'@\(.*?\)',''), ';\n'));
  end

  %% tens3 (D3 f(x))
  % print hessian header
  fprintf(fileID, strcat('\nfunction tens3 = der3', argumentlist,'\n'));
  % create empty matrix to store derivatives
  d3f = sym('a', [n, n, n, n]);
  for i = 1:n
    for j = 1:n
      d3f(:,:,i,j) = diff(d2f(:,:,i), coordinates(j));
      fh = matlabFunction(d3f(:,:,i,j), 'vars', {t, coordinates, parnames{:}});
      tens3ijString = strcat('tens3(:,:,', string(i), ',', string(j), ') = ');
      fprintf(fileID, ...
              strcat(tens3ijString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
    end
  end

  %% tens4 (D4 f(x))
  % print tens4 header
  fprintf(fileID, strcat('\nfunction tens4 = der4', argumentlist,'\n'));
  % create empty matrix to store derivatives
  d4f = sym('a', [n, n, n, n, n]);
  for i = 1:n
    for j = 1:n
      for k = 1:n
        d4f(:,:,i,j,k) = diff(d3f(:,:,i,j), coordinates(k));
        fh = matlabFunction(d4f(:,:,i,j,k), 'vars', {t, coordinates, parnames{:}});
        tens4ijkString = strcat('tens4(:,:,', ...
                                string(i), ',', string(j), ',', string(k), ') = ');
        fprintf(fileID, ...
                strcat(tens4ijkString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
      end
    end
  end

  %% tens5 (D5 f(x))
  % print tens5 header
  fprintf(fileID, strcat('\nfunction tens5 = der5', argumentlist,'\n'));
  for i = 1:n
    for j = 1:n
      for k = 1:n
        for l = 1:n
        d5fijkl = diff(d4f(:,:,i,j,k), coordinates(l));
        fh = matlabFunction(d5fijkl, 'vars', {t, coordinates, parnames{:}});
        tens5ijklString = strcat('tens5(:,:,', ...
                                 string(i),',', string(j), ',', string(k), ...
                                 ',', string(l), ') = ');
        fprintf(fileID, ...
                strcat(tens5ijklString, regexprep(func2str(fh), '@\(.*?\)',''), ';\n'));
        end
      end
    end
  end

  fclose(fileID);
  suc = 1;
end
