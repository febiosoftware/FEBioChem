## Introduction
Consider $N$ chemical species that undergo diffusion and can react with each other, and are inside a flow that moves with velocity $\mathbf{v}$. The governing equation takes on the following form.

\[
\begin{equation}
\label{eq:strong-form}
{{\partial }_{t}}{{u}^{\left( i \right)}}= - \nabla \cdot \left( \mathbf{J}^{(i)}\left(\mathbf{u},\nabla\mathbf{u} \right) \right)-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right),i=1\ldots N
\end{equation}
\]

Here, ${{u}^{\left( i \right)}}$ represents the concentration of chemical species i, $\mathbf{J}^{(i)}$ is the concentration flux, and ${{R}^{\left( i \right)}}$ accounts for all local reactions that contribute to the concentration of species $i$. The vector ${{\mathbf{u}}^{T}}=\left[ {{u}^{\left( 1 \right)}},...,{{u}^{\left( N \right)}} \right]$ represents all the concentrations.

Note that we assume that the concentration flux in general can be a function of both the concentrations $\mathbf{u}$ and their gradients $\nabla \mathbf{u}$. As a simple example, consider the case of Fickian diffusion, where the concentration flux takes on the familiar form.

\[
\mathbf{J}^{(i)} = -{{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}}
\]

Here, ${{\mathbf{D}}^{\left( i \right)}}$ the diffusion tensor of species $i$. 

## Weak Formulation

First, an integral version of this equation must be established. Let’s reorganize the terms as follows,

\[
-{{\partial }_{t}}{{u}^{\left( i \right)}}-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)-\nabla \cdot \left( {{\mathbf{J}}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right)=0
\]

Then, consider $w$ a trial function and integrate over the domain.

\[\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}\left[ -{{\partial }_{t}}{{u}^{\left( i \right)}}-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)-\nabla \cdot \left( {{\mathbf{J}}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right) \right]dV}}=0\]

This can be rewritten as,

\[\begin{align}
  & -\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}{{\partial }_{t}}{{u}^{\left( i \right)}}dV}}-\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)}}-\sum\limits_{i=1}^{N}{\int\limits_{S}{{{w}^{\left( i \right)}}{{{J_n}}^{\left( i \right)}}dS}} \\ 
 & +\sum\limits_{i=1}^{N}{\int\limits_{V}{\nabla {{w}^{\left( i \right)}}\cdot {{\mathbf{J}}^{\left( i \right)}}dV}}+\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}{{R}^{\left( i \right)}}dV}}=0 \\ 
\end{align}\]

Note that the third term is integrated over the surface for which $J_n = \mathbf{J}\cdot \mathbf{n}$ is prescribed.

## Discretization

Next, the finite element approximations are introduced.

\[{{u}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}u_{a}^{\left( i \right)},{{w}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}w_{a}^{\left( i \right)}},{{{\dot{u}}}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}\dot{u}_{a}^{\left( i \right)}}}\]

This results in the following,

\[\begin{align}
  & -\sum\limits_{i=1}^{N}{\sum\limits_{b}{\dot{u}_{b}^{\left( i \right)}\left( \int\limits_{V}{{{N}_{a}}{{N}_{b}}dV} \right)}}-\sum\limits_{i=1}^{N}{\sum\limits_{b}^{{}}{u_{b}^{\left( i \right)}\left( \int\limits_{V}{{{N}_{a}}\nabla \cdot \left( \mathbf{v}{{N}_{b}} \right)}dV \right)}}-\sum\limits_{i=1}^{N}{\int\limits_{S}{{{N}_{a}}{{\mathbf{q}}^{\left( i \right)}}\cdot \mathbf{n}dS}} \\ 
 & +\sum\limits_{i=1}^{N}{\sum\limits_{b}{u_{b}^{\left( i \right)}\left( \int\limits_{V}{\nabla {{N}_{a}}\cdot {{\mathbf{J}}^{\left( i \right)}}dV} \right)}}+\sum\limits_{i=1}^{N}{\int\limits_{V}{{{N}_{a}}{{R}^{\left( i \right)}}dV}}=0 \\ 
\end{align}\]

This can be written as a matrix equation,

\[
\begin{equation}
\label{eq:semi-discrete}
\mathbf{M\dot{U}}+\mathbf{C}\mathbf{U}-\mathbf{D\left(U\right)}=\mathbf{F}\left( \mathbf{U} \right)
\end{equation}
\]

\[\begin{align}
  & {{M}_{pq}}=\int\limits_{V}{{{N}_{a}}{{N}_{b}}dV} \\ 
  &{{C}_{pq}} = \int\limits_{V}{{{N}_{a}}\nabla \cdot \left( \mathbf{v}{{N}_{b}} \right)dV}\\
 & {{D}_{pq}}=\int\limits_{V}{\nabla {{N}_{a}}\cdot {{\mathbf{J}}^{\left( i \right)}}dV} \\ 
 & {{F}_{p}}=\int\limits_{V}{{{N}_{a}}{{R}^{\left( i \right)}}dV} - \sum\limits_{i=1}^{N}{\int\limits_{S}{{{N}_{a}}{{{J_n}}^{\left( i \right)}}dS}}\\ 
\end{align}\]

Here, degrees of freedom are denoted by $p$, where $p$ is the degree of freedom $i$ of node $a$. Similarly, $q$ refers to the degree of freedom $j$ of node $b$.

## Time Integration
Equation (\ref{eq:semi-discrete}) is the semi-discrete equation. Next, a time integration method needs to be introduced. First, we enforce equation (\ref{eq:semi-discrete}) at time $t_{n+1}$. 

\[
\begin{equation}
\mathbf{M\dot{U}}_{n+1}+\mathbf{C}\mathbf{U}_{n+1}-\mathbf{D}\left(\mathbf{U}_{n+1}\right)=\mathbf{F}\left( \mathbf{U}_{n+1} \right)
\end{equation}
\]

Then, we use the following update rules. 

\[
  \begin{align}
  & \mathbf{U}_{n+1}=\mathbf{U}_{n} + \Delta{t}\dot{\mathbf{U}}_{n+\alpha} \\
  & \dot{\mathbf{U}}_{n+\alpha} = \alpha\dot{\mathbf{U}}_{n+1}+\left(1 - \alpha \right)\dot{\mathbf{U}}_{n}
  \end{align}
\]

From this we can find ${{\mathbf{\dot{U}}}_{n+1}}$,

\[\,{{\mathbf{\dot{U}}}_{n+1}}=\frac{1}{\alpha }\left[ \frac{1}{\Delta t}\left( {{\mathbf{U}}_{n+1}}-{{\mathbf{U}}_{n}} \right)-\left( 1-\alpha  \right){{{\mathbf{\dot{U}}}}_{n}} \right]\]

Substitution back into the semi-discrete equation results in the *generalized trapezoidal rule*.

\[\left( \mathbf{M}+\alpha \Delta t\mathbf{C} \right){{\mathbf{U}}_{n+1}}-\alpha \Delta t\mathbf{D}\left( {{\mathbf{U}}_{n+1}} \right)=\alpha \Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)+\mathbf{M}{{\mathbf{U}}_{n}}+\left( 1-\alpha  \right)\Delta t\mathbf{M}{{\mathbf{\dot{U}}}_{n}}\]

which is a nonlinear system of equations. This can be solved using Newton's method. 


### Backward Euler

For $\alpha=1$ , we recover the backward Euler method.

\[\left( \mathbf{M}+\Delta t\mathbf{C} \right){{\mathbf{U}}_{n+1}}-\Delta t\mathbf{D}\left( {{\mathbf{U}}_{n+1}} \right)=\Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)+\mathbf{M}{{\mathbf{U}}_{n}}\]

We solve it using Newton’s method. First, we define the residual.

\[\mathbf{R}\left( {{\mathbf{U}}_{n+1}} \right)=\left( \mathbf{M}+\Delta t\mathbf{C} \right){{\mathbf{U}}_{n+1}}-\Delta t\mathbf{D}\left( {{\mathbf{U}}_{n+1}} \right)-\Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)-\mathbf{M}{{\mathbf{U}}_{n}}\]

A Taylor expansion around the current guess gives,

\[\mathbf{R}\left( \mathbf{U}_{n+1}^{k}+\Delta \mathbf{U} \right)=\mathbf{R}\left( \mathbf{U}_{n+1}^{k} \right)+D\mathbf{R}\left( \mathbf{U}_{n+1}^{k} \right)\Delta \mathbf{U}=\mathbf{0}\]

Here,

\[\begin{align}
   D\mathbf{R}\left( {{\mathbf{U}}_{n+1}} \right)\left[ \Delta \mathbf{u} \right] &=\left( \mathbf{M}+\Delta t\mathbf{C} \right)\Delta \mathbf{U}-\Delta t\left[ D\mathbf{D}\left( {{\mathbf{U}}_{n+1}} \right)\left[ \Delta \mathbf{U} \right]-\frac{d\mathbf{F}}{d\mathbf{U}}\Delta \mathbf{U} \right] \\ 
 & =\mathbf{K}\Delta \mathbf{u}  
\end{align}\]


This can be solved for $\Delta \mathbf{U}$.

\[\mathbf{K}\,\Delta \mathbf{U}=-\mathbf{R}\]

Then, the solution can be updated,

\[\mathbf{U}_{n+1}^{k+1}=\mathbf{U}_{n+1}^{k}+\Delta \mathbf{U}\]

We still need to evaluate $D\mathbf{D}\left( {{\mathbf{U}}_{n+1}} \right)\left[ \Delta \mathbf{U} \right]$. First note that,

\[\delta \mathbf{J}=\frac{\partial \mathbf{J}}{\partial u}\delta u+\frac{\partial \mathbf{J}}{\partial \nabla u}\nabla \left( \delta u \right)\]

Thus,

\[\Delta {{D}_{pq}}=\int\limits_{V}{\nabla {{N}_{a}}\cdot \Delta {{\mathbf{J}}^{\left( i \right)}}dV}=\int\limits_{V}{\nabla {{N}_{a}}\cdot \left( \frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial {{u}^{\left( j \right)}}}\Delta {{u}^{\left( j \right)}}+\frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial \nabla {{u}^{\left( j \right)}}}\nabla \left( \Delta {{u}^{\left( j \right)}} \right) \right)dV}\]

Or finally,

\[\Delta {{D}_{pq}}=\left\{ \int\limits_{V}{\nabla {{N}_{a}}\cdot \left( \frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial {{u}^{\left( j \right)}}}{{N}_{b}}+\frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial \nabla {{u}^{\left( j \right)}}}\nabla {{N}_{b}} \right)dV} \right\}\Delta {{u}^{\left( j \right)}}=\mathbf{D}\,\Delta \mathbf{u}\]

And thus,

\[\mathbf{K}=\left( \mathbf{M}+\Delta t\mathbf{C} \right)-\Delta t\left( \mathbf{D}\,-\frac{d\mathbf{F}}{d\mathbf{U}} \right)\]

## Examples of Diffusion Laws
### Fick's Law
As an example of how to evaluate the necessary derivatives of the concentration flux $\mathbf{J}$, consider again the case of Fickian diffusion. The equation for the flux is repeated here. 

\[{{\mathbf{J}}^{\left( i \right)}}=-\mathbf{D}\cdot \nabla {{u}^{\left( i \right)}}\]

From this, 

\[\frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial {{u}^{\left( i \right)}}}=\mathbf{0},\quad \frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial \nabla {{u}^{\left( i \right)}}}=-\mathbf{D}\]

### Isotropic Diffusivity
Consider now a special case of Fick's law, where the diffusivity is isotropic, but may depend on the concentration.

\[{{\mathbf{J}}^{\left( i \right)}}=-D\left({{u}^{\left( i \right)}} \right)\nabla {{u}^{\left( i \right)}}\]

Here, $D\left( {{u}^{\left( i \right)}} \right)$ is a scalar function of the concentration. Then,

\[\frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial {{u}^{\left( i \right)}}}=-\left( \frac{dD}{d{{u}^{\left( i \right)}}} \right)\nabla {{u}^{\left( i \right)}},\quad \frac{\partial {{\mathbf{J}}^{\left( i \right)}}}{\partial \nabla {{u}^{\left( i \right)}}}=-D\,\mathbf{1}\]

## Chemical Reactions

Thus far, F remained unspecified, or better R. We look now at the specific form this function takes when considering chemical reactions.

A general chemical reaction is specified as,

\[\sum\limits_{i}{{{{{v}'}}_{ij}}}{{M}_{i}}\overset{{{k}_{j}}}{\mathop{\to }}\,\sum\limits_{i}{{{{{v}''}}_{ij}}{{M}_{i}}}\]

Here, ${{M}_{i}}$ denotes the chemical species.

For chemical reactions, $R^{(i)}$ takes on the following form,

\[{{R}^{\left( i \right)}}=\sum\limits_{j=1}^{J}{{{v}_{ij}}{{r}_{j}}}\]

Here, the sum is over $J$ reactions and ${{v}_{ij}}={{{v}''}_{ij}}-{{{v}'}_{ij}}$ and,

\[{{r}_{j}}={{k}_{j}}\prod\limits_{i=1}^{I}{{{\left[ {{u}^{\left( i \right)}} \right]}^{{{\mu }_{ij}}}}}\]

the reaction rate of reaction $j$. Usually, the exponent is taken to be ${{\mu }_{ij}}={{{\nu }'}_{ij}}$, as according to the law of mass action.

Let’s evaluate the stiffness matrix for $\mathbf{F}$.

\[{{k}_{pq}}=\int\limits_{{{V}^{e}}}{{{N}_{a}}\frac{d}{d{{u}_{q}}}{{R}^{\left( i \right)}}\left( \mathbf{u} \right)}dV\]

Evaluating the derivative of $R$ with respect to degree of freedom $q$, we find,

\[\frac{d{{R}^{\left( i \right)}}}{d{{u}_{q}}}={{N}_{b}}{{\gamma }_{ij}}\]

where,

\[{{\gamma }_{ij}}=\sum\limits_{k=1}^{J}{{{v}_{ik}}\frac{d{{r}_{k}}}{d{{u}^{\left( j \right)}}}}\]

Thus, the element stiffness matrix is given by,

\[{{k}_{pq}}=\int\limits_{{{V}^{e}}}{{{\gamma }_{ij}}{{N}_{a}}{{N}_{b}}}dV\]

## Mixtures

Thus far, we considered the case where the solute is free, i.e. not inside a porous material. In the case if the solutes are contained inside a porous material, then the equations above remain valid, except now that the concentrations are with respect to the mixture volume. The 'true' concentration ${{c}^{\left( i \right)}}$, i.e. the concentration w.r.t. to the solvent is related to the 'bulk' concentration via,

\[{{u}^{\left( i \right)}}=\varphi {{c}^{\left( i \right)}}\]

where $\varphi$  is the fluid volume fraction, i.e. the fraction of the volume that is accessible to the fluid (and therefore to the solvents).

In the absence of chemical reactions and convection, equation (\ref{eq:strong-form}) becomes,

\[{{\partial }_{t}}\left( \varphi {{c}^{\left( i \right)}} \right)=\nabla \cdot \left( \varphi {{\mathbf{D}}^{\left( i \right)}}\nabla {{c}^{\left( i \right)}} \right),i=1\ldots N\]

Note that $\varphi$ can be dependent on the position, so in general does not drop out of the equation.

Note that now $\mathbf{D}$ is the diffusivity in the porous mixture and can be different (i.e. lower) than the free diffusivity.

Performing the usual steps,

\[-\sum\limits_{i}{\int\limits_{V}{{{w}^{\left( i \right)}}{{\partial }_{t}}\left( \varphi {{u}^{\left( i \right)}} \right)dV}}-\sum\limits_{i}{\int\limits_{V}{\nabla {{w}^{\left( i \right)}}\cdot \left( \varphi {{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}} \right)dV}}-\sum\limits_{i}{\int\limits_{S}{{{w}^{\left( i \right)}}\mathbf{q}_{0}^{\left( i \right)}dS}}=0\]

where ${{\mathbf{q}}^{\left( i \right)}}=-\varphi {{\mathbf{D}}^{\left( i \right)}}\cdot \nabla {{\mathbf{c}}^{\left( i \right)}}$.  Introducing discretization,

\[-\sum\limits_{i}{\sum\limits_{b}{\dot{c}_{b}^{\left( i \right)}\int\limits_{{{V}^{e}}}{\varphi {{N}_{a}}{{N}_{b}}dV}}}-\sum\limits_{i}{\sum\limits_{b}{c_{b}^{\left( i \right)}\int\limits_{V}{\nabla {{N}_{a}}\cdot \varphi {{\mathbf{D}}^{\left( i \right)}}\nabla {{N}_{b}}dV}}}-\sum\limits_{i}{\int\limits_{S}{{{w}^{\left( i \right)}}\mathbf{q}_{0}^{\left( i \right)}dS}}=0\]

### Solid-bound Molecules

In the above we assumed that the domain was composed of a mixture, defined by two phases, a solid and a fluid phase. The fluid was assumed to be a solvent in which several species were dissolved. We emphasize again that the volume fractions of the dissolved solutions are negligible and thus that the solid and fluid volume fractions add up to unity.

\[{{\varphi }^{s}}+{{\varphi }^{f}}=1\]

We now consider the case that a chemical species exists as part of the solid phase. Such a species is called a _solid-bound molecule_ (or _sbm_). Since it is bound to the solid (immovable) phase, it is assumed that a solid-bound molecule does not undergo diffusion. Let us denote the apparent density of the sbm by ${{\rho }^{\sigma }}$. It then follows that the governing equation for a solid-bound molecule is given by,

\[\frac{d{{\rho }^{\sigma }}}{dt}={{\hat{\rho }}^{\sigma }}\]

where ${{\hat{\rho }}^{\sigma }}$ is the mass supply for the sbm. In this context, the mass supply will be calculated from the chemical reactions.

Time integration is done via the backward Euler method.

\[\rho _{t+\Delta t}^{\sigma ,k+1}=\rho _{t}^{\sigma }+\Delta t\hat{\rho }_{t+\Delta t}^{\sigma ,k}\]

Alternatively, the trapezoidal rule could be used,

\[\rho _{t+\Delta t}^{\sigma ,k+1}=\rho _{t}^{\sigma }+\frac{\Delta t}{2}\left( \hat{\rho }_{t+\Delta t}^{\sigma ,k}+\hat{\rho }_{t}^{\sigma } \right)\]

A change in the apparent density of an SBM has an effect of the solid volume fraction, which is now defined via,

\[{{\varphi }^{s}}=\varphi _{0}^{s}+\sum\limits_{\sigma }{{{\varphi }^{\sigma }}}\]

where the volume fraction of the sbm is defined via,

\[{{\varphi }^{\sigma }}=\frac{{{\rho }^{\sigma }}}{\rho _{T}^{\sigma }}\]

That is, it is the ratio of the apparent density over the true density of the sbm. The fluid volume fraction, or porosity, is then defined via $\varphi =1-{{\varphi }^{s}}$.

