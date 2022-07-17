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

# Generate system files for Morris-Lecar in a platinum model

In this script the __system files__ for the Morris-Lecar model {cite}`Morris1981`

$$
\begin{cases}
\begin{aligned}
c\dot V &= I_{app} - I_{ion}, \\
\dot w &= \phi \frac{w_\infty -w}{\tau},
\end{aligned}
\end{cases}
$$

where

$$
\begin{aligned}
m_\infty &= 0.5 \left(1+\tanh\left(\frac{v-v_1}{v_2}\right)\right), \\
w_\infty &= 0.5 \left(1+\tanh\left(\frac{v-v_3}{v_4}\right)\right), \\
i_{ion} &= g_{ca} m_{\infty} (v-v_{ca}) + g_k w (v-v_k) + g_l (v-v_l), \\
\tau &= \text{sech}\left(\frac{v-v_3}{2v_4}\right).
\end{aligned}
$$

are generated. These are used in the [Morris-Lecar demo](./Morris-Lecar.ipynb).


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
system_name = 'Morris_Lecar';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'V', 'w'};
parnames = {'Iapp', 'v3'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.  Due to
the following two lines we may, for example, use either `alpha` or `par(1)`.
There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define fixed parameters

```matlab
C=20;
VL=-60;
VCa=120;
VK=-84;
gL=2;
gCa=4.4;
gK=8;
v1=-1.2;
v2=18;
v4=30;
phi=1/25;
```

## Define the system

```matlab
minf = 0.5*(1+tanh((V-v1)/v2));
winf = 0.5*(1+tanh((V-v3)/v4));
Iion = gCa*minf*(V-VCa) + gK*w*(V-VK) + gL*(V-VL);
tau = sech((V-v3)/2/v4);
dV_dt = (Iapp - Iion)/C;
dw_dt = phi/tau*(winf-w);
system = [dV_dt; dw_dt];
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
