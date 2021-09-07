---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.8.0
  kernelspec:
    display_name: Matlab
    language: matlab
    name: matlab
---

# Generate system files for predator-prey system

This script generates the __system files__ for the following predator-prey
ecosystem

\begin{cases}
\begin{aligned}
\dot x_1 &= x_1 - \frac{x_1 x_2}{1+\alpha x_1} - \epsilon x_1^2, \\
\dot x_2 &= -\gamma x_2 + \frac{x_1 x_2}{1+\alpha x_1} - \delta x_2^2, \\
\end{aligned}
\end{cases}

These are used in the [Predator-prey dome](Bazykin.ipynb).


## Add MatCont path and load sym package if GNU Octave is used

```matlab
matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, 'Utilities']);
cd(matcontpath) % this file should be executed form the main directory of MatCont
if isOctave
  pkg load symbolic % for GNU Octave
end
```

## Define the system name

```matlab
system_name = 'Bazykin';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'x1', 'x2'};
parnames = {'alpha', 'delta'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `alpha` or `par(1)`. There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define fixed parameters

```matlab
gamma = 1;
epsilon = 0.01;
```

## Define the system

```matlab
dx1_dt = x1 - x1*x2/(1+alpha*x1) - epsilon*x1^2;
dx2_dt = -gamma*x2 + x1*x2/(1+alpha*x1) - delta*x2^2;
system = [dx1_dt; dx2_dt];
```

There are no modifications needed after this line.

## Differentiate and generate code, exporting it to system_name.m

```matlab
suc = generate_directional_derivatives(...
  system,...   % n x 1 array of derivative symbolic expressions
  coords,... % 1 x n array of symbols for states
  par,...      % 1 x np (or np x 1) array of symbols used for parameters
  system_name... % argument specifying the system name
);
```

## Higher-order parameter-dependent multi-linear form.

Exporting it to system_name_multilinearforms.m. These multi-linear forms are
currently only used in the computation of the parameter-dependent center
manifold for the codimension two Bogdanov-Takens bifurcation.

```matlab
order = 3;
suc = generate_multilinear_forms(system_name, system, coords, par, order);
```
