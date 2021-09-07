function suc = generate_multilinear_forms(system_name, system, xx, pars, order, outputdir)

% If we are using GNU Octave call octave version of this file,
% otherwise continue
if isOctave
  suc = generate_multilinear_forms_octave(system_name, system, xx, pars, order, outputdir);
  return
end

dim_phase = length(xx); % phase dimension
dim_par = length(pars); % parameter dimension

%% vectors
v1 = sym('v1', [dim_phase, 1]);
v2 = sym('v2', [dim_phase, 1]);
v3 = sym('v3', [dim_phase, 1]);
p1 = sym('p1', [dim_par, 1]);
p2 = sym('p2', [dim_par, 1]);
p3 = sym('p3', [dim_par, 1]);

%% calculuate multilinear forms
A  = jacobian(system, xx);
J1 = jacobian(system, pars);
B  = jacobian(A*v1, xx)*v2;
J2 = jacobian(J1*p1, pars)*p2;
A1 = jacobian(A*v1, pars)*p1;
if order > 2
    C  = jacobian(B, xx)*v3;
    J3 = jacobian(J2, pars)*p3;
    B1 = jacobian(B, pars)*p1;
    A2 = jacobian(A1, pars)*p2;
end

%% create system directory if not exist
system_dir = [outputdir, system_name,'/'];
if ~exist(system_dir,'dir')
   mkdir(system_dir);
end

%% temporally write multilinear forms to files in the directory of the defined system
matlabFunction(A,  'vars', {xx, pars}, 'file', [system_dir,'A.m']);
matlabFunction(J1, 'vars', {xx, pars}, 'file', [system_dir,'J1.m']);
matlabFunction(B,  'vars', {xx, pars, v1, v2}, 'file', [system_dir,'B.m']);
matlabFunction(A1, 'vars', {xx, pars, v1, p1}, 'file', [system_dir,'A1.m']);
matlabFunction(J2, 'vars', {xx, pars, p1, p2}, 'file', [system_dir,'J2.m']);
if order > 2
    matlabFunction(C,  'vars', {xx, pars, v1, v2, v3}, 'file', [system_dir,'C.m']);
    matlabFunction(J3, 'vars', {xx, pars, p1, p2, p3}, 'file', [system_dir,'J3.m']);
    matlabFunction(B1, 'vars', {xx, pars, v1, v2, p1}, 'file', [system_dir,'B1.m']);
    matlabFunction(A2, 'vars', {xx, pars, v1, p1, p2}, 'file', [system_dir,'A2.m']);
end

%% export all files to a single file
% open file to write to and print header to file
fidOut=fopen([outputdir, system_name, '_multilinearforms.m'], 'w');
fprintf(fidOut, strcat('function F = %s_multilinearforms\n', ...
                   'F.A = @A;\n', ...
                   'F.J1 = @J1;\n', ...
                   'F.B = @B;\n', ...
                   'F.A1 = @A1;\n', ...
                   'F.J2 = @J2;\n'), system_name);
if order > 2
    fprintf(fidOut, strcat('F.C = @C;\n', ...
                           'F.J3 = @J3;\n', ...
                           'F.B1 = @B1;\n', ...
                           'F.A2 = @A2;\n\n'));
else
    fprintf(fidOut, '\n');
end

files = {'A.m', 'J1.m', 'B.m', 'J2.m', 'A1.m'};
if order > 2
    files = {'A.m', 'J1.m', 'B.m', 'J2.m', 'A1.m', 'C.m', 'J3.m', 'B1.m', 'A2.m'};
end
numberOfFiles = length(files);
for i=1:numberOfFiles
    fid=fopen([system_dir, files{i}],'r');      % open each file in turn
    fwrite(fidOut,fread(fid,'*char'),'*char');  % read remainder as char* image and echo back out
    fwrite(fidOut,newline);
    fclose(fid);                                % done with that one ...
    delete([system_dir, files{i}]);             % remove temporally created file
end
fclose(fidOut);                                 % and close the output

suc = 1;
end
