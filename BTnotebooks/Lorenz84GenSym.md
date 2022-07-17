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

# Generate system files for the extended Lorenz-84 model

This Jupyter Notebook generates the __system files__ for the extended Lorenz-84 model
given by

$$
\begin{cases}
\begin{aligned}
\dot{X} &=-Y^{2}-Z^{2}-\alpha X+\alpha F-\xi U^{2} \\
\dot{Y} &=X Y-\beta X Z-Y+G \\
\dot{Z} &=\beta X Y+X Z-Z \\
\dot{U} &=-\delta U+\xi U X+S
\end{aligned}
\end{cases}
$$

In this system, $X$ models the intensity of a baroclinic wave, $Y$ and $Z$ the
sin and cos coefficients of the wave respectively, the variable $U$ is added to
study the influence of external parameters such as temperature. We fix the
parameters as follows $\alpha=0.25, \beta=1, G=0.25, \delta=1.04$ and
$\xi=0.987$. The continuation parameters are $F$ and $S$.


These are used in the [extended Lorenz-84 model](./extendedLorenz84model.ipynb) demo.

## Add MatCont path and load sym package if GNU Octave is used


```matlab
matcontpath = '../';
addpath(matcontpath);
addpath([matcontpath, 'Utilities'])
if isOctave
  pkg load symbolic % for GNU Octave
end
```

## Set the system name

```matlab
system_name = 'extendedLorenz84';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'X', 'Y', 'Z', 'U'};
parnames = {'F', 'S'};
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
alpha = 0.25;
beta = 1;
G = 0.25;
delta = 1.04;
xi = 0.987;
```

## Define the system

```matlab
dX_dt = -Y^2-Z^2-alpha*X+alpha*F-xi*U^2;
dY_dt = X*Y-beta*X*Z-Y+G;
dZ_dt = beta*X*Y+X*Z-Z;
dU_dt = -delta*U+xi*U*X+S;
system = [dX_dt; dY_dt; dZ_dt; dU_dt];
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
