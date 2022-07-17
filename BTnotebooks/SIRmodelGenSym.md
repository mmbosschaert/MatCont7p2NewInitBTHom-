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

# Generate system files for SIR model with the standard incidence rate

This Jupyter Notebook generates the __system file__ for the following SIR model
{cite}`Chunhua@2014` with the standard incidence rate

$$
\begin{cases}
\begin{aligned}
    \frac{d S}{d t} &=A-d S-\frac{\beta S I}{S+I+R}, \\ 
    \frac{d I}{d t} &=-(d+\nu) I-\mu(b,I) I+\frac{\beta S I}{S+I+R}, \\ 
    \frac{d R}{d t} &=\mu(b,I) I-d R,
\end{aligned}
\end{cases}
$$

where $A>0$ is the recruitment rate of susceptible population;
$d>0$ is the per capita natural death rate of the population; $\nu>0$ is
the per capita disease-induced death rate; $\mu>0$ is the per capita
recovery rate of infectious individual. The funtion $\mu$ is set to

$$
\mu(b,I) = \mu_0 + (\mu_1 - \mu_0) \frac{b}{I+b}.
$$

These are used in the [SIR model](SIRmodel.ipynb) demo.

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
system_name = 'SIRmodel';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'S', 'I', 'R'};
parnames = {'mu1', 'b'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `mu1` or
`par(1)`. There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define fixed parameters

```matlab
A = 20;
mu0 = 10;
d = 1/10;
nu = 1;
beta = 11.5;
```

## Define the system

```matlab
mu = @(b,I) mu0 + (mu1-mu0)*b/(I+b);
dS_dt = A-d*S-beta*S*I/(S+I+R);
dI_dt = -(d+nu)*I-mu(b,I)*I+beta*S*I/(S+I+R);
dR_dt = mu(b,I)*I-d*R;
system = [dS_dt; dI_dt; dR_dt];
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
