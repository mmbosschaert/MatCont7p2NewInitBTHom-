---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Matlab
    language: matlab
    name: matlab
---

# Generate system files for CO-oxidation in a platinum model

In this script the __system files__ for the CO-oxidation model

$$
\begin{cases}
\begin{aligned}
z &= 1 - x - y - s, \\
\dot x &= 2k_1z^2 - 2k_{-1}x^2 - k_3xy, \\
\dot y &= k_2z - k_{-2}y - k_3xy, \\
\dot s &= k_4(z - \lambda s).
\end{aligned}
\end{cases}
$$

are generated. These are used in the [CO-oxidation.ipynb](./CO-oxidation.ipynb)
demo.


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
system_name = 'CO_oxidation';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'x', 'y', 's'};
parnames =  {'k1', 'km1', 'k3', 'k2', 'km2', 'k4', 'lambda'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `k1` or
`par(1)`. There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define the system

```matlab
z = 1-x-y-s;
dx_dt = 2*k1*z^2 - 2*km1*x^2 - k3*x*y;
dy_dt = k2*z - km2*y - k3*x*y;
ds_dt = k4*(z - lambda*s);
system = [dx_dt; dy_dt; ds_dt];
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
