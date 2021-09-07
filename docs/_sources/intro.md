# Homoclinic Predictors for Bogdanov-Takens Bifurcation points in ODEs

In this Jupyter Book we will demonstrate the effectiveness of the homoclinic
predictors derived in {cite}`Bosschaert@2021` to start continuations of
homoclinic orbits in two parameters in MatCont near generic codimension two
Bogdanov-Takens bifurcation points in $n$-dimensional ordinary differential
equations (ODEs) of the form
```{math}
:label: ODE
\dot x(t) = f(x(t), \alpha),
```
where $x: \mathbb R \to \mathbb R^n$ and $\alpha \in \mathbb R^m$. The
function $f : \mathbb R^n \times \mathbb R^m \to \mathbb R^n$ is assumed to
be as smooth as necessary.

In total there are eight different models considered. Each the models consists
of two notebooks:
- one to generate the necessary _derivative files_, see for example
  [](./Morris-LecarGenSym.ipynb),
- and one _model notebook_ to analyse the model.

In the first two models [](./Morris-Lecar.ipynb),  and [](./CO-oxidation.ipynb)
we partly follow the presentation in {cite}`Bashir@2015`, that is, we define an
equilibrium of {eq}`ODE`, continue the equilibrium point one parameter, encounter
a codimension one Hopf or fold bifurcation, continue the codimension one
bifurcation point in two parameters and encounter one or more Bogdanov-Takens
points. The homoclinic orbits emanating from these Bogdanov-Takens point can
then be continued using one of the homoclinic predictors form
{cite}`Bosschaert@2021` implemented into MatCont. In these first two models we
will describe in detail how a homoclinic approximation near the generic
codimension two Bogdanov-Takens bifurcation is obtained in MatCont.

Although the procedure described above is a common situation for discovering
Bogdanov-Takens points in ODEs, in many situations it is possible to derive the
Bogdanov-Takens point directly, either analytically or numerically, form the
ODE. This will be the approach for the remaining models were we will focus
solely on the homoclinic orbits emanating from the Bogdanov-Takens point and
not the Hopf and fold curves.

At the end of each _model notebook_ we create a convergence plot of the
different methods to approximation the homoclinic solution derived in
{cite}`Bosschaert@2021`.

In the notebook [BogdanovTakens.ipynb](BogdanovTakens.ipynb) the importance of
using a higher order approximation of the non-linear time transformation used
in {cite}`Bosschaert@2021` is demonstrated.


```{note}
The _derivative files_ can also be generated with the graphical user interface of
MatCont. In fact, if preferred, the continuation of the various bifurcation
branches can be done using only the graphical user interface.
```
