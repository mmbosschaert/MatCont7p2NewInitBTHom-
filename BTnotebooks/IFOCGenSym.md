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

# Generate system files for the indirect field oriented control system

This Jupyter Notebook generates the __system files__ for the indirect field oriented control system given by

$$
\begin{cases}
\begin{aligned}
\dot x_1 &= -c_1 x_1 + c_2 x_4 - \frac{k c_1}{u_2^0} x_2 x_4, \\
\dot x_2 &= -c_1 x_2 + c_2 u_2^0 + \frac{k c_1}{u_2^0} x_1 x_4, \\
\dot x_3 &= -c_3 x_3 - c_4 c_5 (x_2x_4 - u_2^0 x_1) + (c_4 T_m + c_3 w_{ref}), \\
\dot x_4 &= -(k_i-k_p)x_3 - k_p c_4 c_5 ( x_2 x_4 - u_2^0 x_1) + k_p (c_4 T_m + c_3 w_{ref}). \\
\end{aligned}
\end{cases}
$$

Here $x_1 , x_2 , x_3$ and $x_4$ are the state variables, where $x_1$ and $x_2$
represent, respectively, direct and quadrature components of the rotor ﬂux;
$x_3$ is the rotor speed error; and $x_4$ denotes the quadrature axis component
of the stator current, respectively.  We also deﬁne the following constants and
parameters: $u_2^0$ is a constant reference for the rotor flux magnitude; $c_1$
to $c_5$ are machine parameters; $k_p$ and $k_i$ are the proportional (P) and
the integral (I) control gains, respectively; $w_{ref}$ is the speed reference;
$T_m$ the load torque; $k$ the measure of rotor time constant mismatches.


These are used in the [IFOC](IFOC.ipynb) demo.


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
system_name = 'IFOC';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'x1', 'x2', 'x3', 'x4'};
parnames={'k', 'Tm'};
```

## Create symbols for coordinates and parameters
The array `par` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `k` or
`par(1)`. There should no changes be need of this code.

```matlab
syms(parnames{:});       % create symbol for alpha and delta
par=cell2sym(parnames);  % now alpha1 is par(1) etc
syms(coordsnames{:});    % create symbol for alpha and delta
coords=cell2sym(coordsnames); % create 1 x n vector for coordinates
```

## Define fixed parameters

```matlab
c1 = 4.4868;
c2 = 0.3567;
c3 = 0;
c4 = 9.743;
c5 = 1.911;
u20 = 11.3;
kp = 4.5;
ki = 500;
wref = 0;
```

## Define the system

```matlab
dx1_dt = -c1*x1 + c2*x4 - k*c1/u20*x2*x4;
dx2_dt = -c1*x2 + c2*u20 + k*c1/u20*x1*x4;
dx3_dt = -c3*x3 - c4*c5*(x2*x4 - u20*x1) + (c4*Tm + c3*wref);
dx4_dt = -(ki-kp*c3)*x3 - kp*c4*c5*(x2*x4 - u20*x1) + kp*(c4*Tm + c3*wref);
system = [dx1_dt; dx2_dt; dx3_dt; dx4_dt];
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
