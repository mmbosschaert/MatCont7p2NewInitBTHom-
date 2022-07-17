---
jupyter:
  jupytext:
    cell_metadata_filter: -all
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

# Generate system files for Hodgkin-Huxley equations

The Hodgkin-Huxley equations {cite}`HodgkinHuxley@1952` relate the
difference in electric potential across the cell membrane $(V)$ and gating
variables $(m, n$ and $h$ ) for ion channels to the stimulus intensity $(I)$
and temperature $(T)$, as follows:
 
$$
\begin{cases}
\dot{V} ={}& -G(V, m, n, h)+I \\
\dot{m} ={}& \Phi(T)\left[(1-m) \alpha_{m}(V)-m \beta_{m}(V)\right] \\
\dot{n} ={}& \Phi(T)\left[(1-n) \alpha_{n}(V)-n \beta_{n}(V)\right] \\
\dot{h} ={}& \Phi(T)\left[(1-h) \alpha_{h}(V)-h \beta_{h}(V)\right]
\end{cases}
$$

where $\dot{x}$ stands for $\mathrm{d} x / \mathrm{d} t$ and $\Phi$ is
given by $\Phi(T)=3^{(\mathrm{T}-6.3) / 10}$. The other functions involved
are:

$$
G(V, m, n, h)=\bar{g}_{\mathrm{Na}} m^{3}
h\left(V-\bar{V}_{\mathrm{Na}}\right)+\bar{g}_{\mathrm{K}}
n^{4}\left(V-\bar{V}_{\mathrm{K}}\right)+\bar{g}_{\mathrm{L}}\left(V-\bar{V}_{\mathrm{L}}\right)
$$

and the equations modeling the variation of membrane permeability are:

$$
\begin{array}{ll}
\alpha_{m}(V)=\Psi\left(\frac{V+25}{10}\right) & \beta_{m}(V)=4 e^{V / 18}
\\
\alpha_{n}(V)=0.1 \Psi\left(\frac{V+10}{10}\right) & \beta_{n}(V)=0.125
e^{V / 80} \\
\alpha_{h}(V)=0.07 e^{V / 20} & \beta_{h}(V)=\left(1+e^{(V+30) /
10}\right)^{-1}
\end{array}
$$

with

$$
\Psi(x)=\left\{\begin{array}{ll}
x /\left(e^{x}-1\right) & \text { if } x \neq 0 \\
1 & \text { if } x=0
\end{array}\right.
$$

The parameters $\bar{g}_{\text {ion }}$ and $\bar{V}_{\text {ion}}$
representing maximum conductance and equilibrium potential for the ion were
obtained from experimental data by Hodgkin and Huxley, with the values given
below:

$$
\begin{array}{lll}
\bar{g}_{\mathrm{Na}}=120 \mathrm{mS} / \mathrm{cm}^{2}, &
\bar{g}_{\mathrm{K}}=36 \mathrm{mS} / \mathrm{cm}^{2}, &
\bar{g}_{\mathrm{L}}=0.3 \mathrm{mS} / \mathrm{cm}^{2} \\
\bar{V}_{\mathrm{Na}}=-115 \mathrm{mV}, & \bar{V}_{\mathrm{K}}=12
\mathrm{mV}, & \bar{V}_{\mathrm{L}}=10.599 \mathrm{mV}
\end{array}
$$

The values of $\bar{V}_{\mathrm{Na}}$ and $\bar{V}_{\mathrm{K}}$ can be
controlled experimentally {cite}`HodgkinHuxley@1952a`.
The temperature is set to $T=6.3^{\circ}$.

These are used in the [HodgkinHuxley](HodgkinHuxley.ipynb).

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
system_name = 'HodgkinHuxley';
```

## Create coordinates and parameter names as strings 

```matlab
coordsnames = {'V', 'm', 'n', 'h'};
parnames={'VbarK', 'I'};
```

## Create symbols for parameters
The array `|par|` is the array of symbols in the same order as parnames.
Due to the following two lines we may, for example, use either `mu_1` or
`par(1)`

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
gbarNa = 120;
gbarK  = 36;
gbarL  = 0.3;
VbarNa = -115;
VbarL  = 10.599;
T = 6.3;
```

## Define the system

```matlab
Psi = @(x) x/(exp(x)-1);
alpha_m = @(V) Psi( (V+25)/10 );
alpha_n = @(V) 0.1*Psi( (V+10)/10);
alpha_h = @(V) 0.07*exp(V/20);

beta_m = @(V) 4*exp(V/18);
beta_n = @(V) 0.125*exp(V/80);
beta_h = @(V) 1/(1+exp((V+30)/10));

G = @(V, m, n, h) gbarNa*m^3*h*(V-VbarNa) + gbarK*n^4*(V-VbarK) + gbarL*(V-VbarL);
Phi = @(T) 3^(T-6.3)/10;
dV_dt = -G(V, m, n, h)+I;
dm_dt = Phi(T)*((1-m)*alpha_m(V)-m*beta_m(V));
dn_dt = Phi(T)*((1-n)*alpha_n(V)-n*beta_n(V));
dh_dt = Phi(T)*((1-h)*alpha_h(V)-h*beta_h(V));
system = [dV_dt; dm_dt; dn_dt; dh_dt];
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
