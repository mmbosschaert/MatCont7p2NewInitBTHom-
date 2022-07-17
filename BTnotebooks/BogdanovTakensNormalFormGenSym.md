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

In this script the __system files__ for critical Bogdanov-Takens codimension
two normal form given by

$$
\begin{cases}
\begin{aligned}
\dot w_0 &= w_1, \\
\dot w_1 &= \beta_1  + \beta_2 w_1 - w_0^2 + w_0 w_1. \\
\end{aligned}
\end{cases}
$$

are generated. These are used in the [BogdanovTakens.md](BogdanovTakens.md) demo.


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
system_name = 'BogdanovTakensNormalForm';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'w0', 'w1'};
parnames = {'beta1', 'beta2'};
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

## Define the system

```matlab
dw0_dt = w1; 
dw1_dt = beta1 + beta2*w1 - w0^2 + w0*w1;
system = [dw0_dt; dw1_dt];
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
