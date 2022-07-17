---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: Matlab
    language: matlab
    name: matlab
---

# Generate system files for predator-prey system with constant rate harvesting

In this script the __system files__ for the predator-prey system given by

$$
\begin{cases}
\begin{aligned}
\dot x &= rx\left(1-\frac x k\right) - y \frac x {e+x}, \\
\dot y &= y\left(-d+\frac x {e+x}\right) - h, \\
\end{aligned}
\end{cases}
$$

are generated. These are used in the [predator prey](./PredatorPrey.ipynb) demo.


## Add MatCont path and load sym package if GNU Octave is used

```matlab
matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, '/Utilities']);
if isOctave
  pkg load symbolic % for GNU Octave
end
```

## Set the system name

```matlab
system_name = 'PredatorPrey';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'x', 'y'};
parnames = {'d', 'h'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `alpha` or
`par(1)`. There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define fixed parameters

```matlab
r = 1;
e = 1;
k = 2;
```

## Define the system

```matlab
dx_dt = r*x*(1-x/k) - y*x/(e+x);
dy_dt = y*(-d+x/(e+x)) - h;
system = [dx_dt; dy_dt];
```

In general there are no modifications needed after this line.

## Differentiate and generate code (directional derivatives)

Exporting it to `<system_name>.m`. This method uses directional derivatives.
Then using polarization identities derivatives can be calculated in arbitrary
direction.

```matlab
suc = generate_directional_derivatives(...
  system,...   % n x 1 array of derivative symbolic expressions
  coords,... % 1 x n array of symbols for states
  par,...      % 1 x np array of symbols used for parameters
  system_name,... % argument specifying the system name
  [matcontpath, 'Systems/']... % directory to save to file
);
```

## Higher-order parameter-dependent multi-linear form.

Exporting it to `<system_name>_multilinearforms.m`. These multi-linear forms are
currently only used in the computation of the parameter-dependent center
manifold for the codimension two Bogdanov-Takens bifurcation.

```matlab
order = 3;
suc = generate_multilinear_forms(system_name, system, coords, par, order, ...
        [matcontpath, 'Systems/']);
```
