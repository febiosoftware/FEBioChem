## Introduction
Consider $N$ chemical species that undergo diffusion and can react with each other, and are inside a flow that moves with velocity $\mathbf{v}$. The governing equation takes on the following form.

\[
\begin{equation}
\label{eq:strong-form}
{{\partial }_{t}}{{u}^{\left( i \right)}}=\nabla \cdot \left( {{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}} \right)-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right),i=1\ldots N
\end{equation}
\]

Here, ${{u}^{\left( i \right)}}$ represents the concentration of chemical species i, ${{\mathbf{D}}^{\left( i \right)}}$ the diffusion tensor of species $i$, and ${{R}^{\left( i \right)}}$ accounts for all local reactions that contribute to the concentration of species $i$. The vector ${{\mathbf{u}}^{T}}=\left[ {{u}^{\left( 1 \right)}},...,{{u}^{\left( N \right)}} \right]$ represents all the concentrations.

## Weak Formulation

First, an integral version of this equation must be established. Let’s reorganize the terms as follows,

\[
-{{\partial }_{t}}{{u}^{\left( i \right)}}-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)+\nabla \cdot \left( {{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right)=0
\]

Then, consider $w$ a trial function and integrate over the domain.

\[\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}\left[ -{{\partial }_{t}}{{u}^{\left( i \right)}}-\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)+\nabla \cdot \left( {{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}} \right)+{{R}^{\left( i \right)}}\left( \mathbf{u} \right) \right]dV}}=0\]

This can be rewritten as,

\[\begin{align}
  & -\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}{{\partial }_{t}}{{u}^{\left( i \right)}}dV}}-\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}\nabla \cdot \left( \mathbf{v}{{u}^{\left( i \right)}} \right)}}+\sum\limits_{i=1}^{N}{\int\limits_{S}{{{w}^{\left( i \right)}}{{\mathbf{q}}^{\left( i \right)}}\cdot \mathbf{n}dS}} \\ 
 & -\sum\limits_{i=1}^{N}{\int\limits_{V}{\nabla {{w}^{\left( i \right)}}\cdot {{\mathbf{D}}^{\left( i \right)}}\nabla {{u}^{\left( i \right)}}dV}}+\sum\limits_{i=1}^{N}{\int\limits_{V}{{{w}^{\left( i \right)}}{{R}^{\left( i \right)}}dV}}=0 \\ 
\end{align}\]

Here, $\mathbf{q}=-\mathbf{D}\nabla u$, is the flux vector and the second term is integrated over the surface for which $\mathbf{q}$ is prescribed.

## Discretization

Next, the finite element approximations are introduced.

\[{{u}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}u_{a}^{\left( i \right)},{{w}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}w_{a}^{\left( i \right)}},{{{\dot{u}}}^{\left( i \right)}}=\sum\limits_{a}{{{N}_{a}}\dot{u}_{a}^{\left( i \right)}}}\]

This results in the following,

\[\begin{align}
  & -\sum\limits_{i=1}^{N}{\sum\limits_{b}{\dot{u}_{b}^{\left( i \right)}\left( \int\limits_{V}{{{N}_{a}}{{N}_{b}}dV} \right)}}-\sum\limits_{i=1}^{N}{\sum\limits_{b}^{{}}{u_{b}^{\left( i \right)}\left( \int\limits_{V}{{{N}_{a}}\nabla \cdot \left( \mathbf{v}{{N}_{b}} \right)}dV \right)}}+\sum\limits_{i=1}^{N}{\int\limits_{S}{{{N}_{a}}{{\mathbf{q}}^{\left( i \right)}}\cdot \mathbf{n}dS}} \\ 
 & -\sum\limits_{i=1}^{N}{\sum\limits_{b}{u_{b}^{\left( i \right)}\left( \int\limits_{V}{\nabla {{N}_{a}}\cdot {{\mathbf{D}}^{\left( i \right)}}\nabla {{N}_{b}}dV} \right)}}+\sum\limits_{i=1}^{N}{\int\limits_{V}{{{N}_{a}}{{R}^{\left( i \right)}}dV}}=0 \\ 
\end{align}\]

This can be written as a matrix equation,

\[
\begin{equation}
\label{eq:semi-discrete}
\mathbf{M\dot{U}}+\mathbf{DU}=\mathbf{F}\left( \mathbf{U} \right)
\end{equation}
\]

\[\begin{align}
  & {{M}_{pq}}=\int\limits_{V}{{{N}_{a}}{{N}_{b}}dV} \\ 
 & {{D}_{pq}}=\int\limits_{V}{\nabla {{N}_{a}}\cdot {{\mathbf{D}}^{\left( i \right)}}\nabla {{N}_{b}}dV}+\int\limits_{V}{{{N}_{a}}\nabla \cdot \left( \mathbf{v}{{N}_{b}} \right)dV} \\ 
 & {{F}_{p}}=\sum\limits_{i=1}^{N}{\int\limits_{S}{{{N}_{a}}{{\mathbf{q}}^{\left( i \right)}}\cdot \mathbf{n}dS}}+\int\limits_{V}{{{N}_{a}}{{R}^{\left( i \right)}}dV} \\ 
\end{align}\]

Here, degrees of freedom are denoted by $p$, where $p$ is the degree of freedom $i$ of node $a$. Similarly, $q$ refers to the degree of freedom $j$ of node $b$.

Equation (\ref{eq:semi-discrete}) is the semi-discrete equation. Next, a time integration method needs to be introduced. Using the generalized trapezoidal rule this can be written as,

\[\left( \mathbf{M}+\alpha \Delta t\mathbf{D} \right){{\mathbf{U}}_{n+1}}=\left( \mathbf{M}-\left( 1-\alpha  \right)\Delta t\mathbf{D} \right){{\mathbf{U}}_{n}}+\Delta t\left( \alpha {{\mathbf{F}}_{n+1}}+\left( 1-\alpha  \right){{\mathbf{F}}_{n}} \right)\]

which is a nonlinear system of equations.

## Time Integration

### Backward Euler

For $\alpha=1$ , we recover the backward Euler method.

\[\left( \mathbf{M}+\Delta t\mathbf{D} \right){{\mathbf{U}}_{n+1}}=\mathbf{M}{{\mathbf{U}}_{n}}+\Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)\]

We solve it using Newton’s method. First, we define the residual.

\[\begin{align}
  & \mathbf{R}\left( {{\mathbf{U}}_{n+1}} \right)=\mathbf{\tilde{K}}{{\mathbf{U}}_{n+1}}-\Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)-\mathbf{M}{{\mathbf{U}}_{n}} \\ 
 & =\mathbf{M}\left( {{\mathbf{U}}_{n+1}}-{{\mathbf{U}}_{n}} \right)+\Delta t\mathbf{D}{{\mathbf{U}}_{n+1}}-\Delta t\mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)  
\end{align}\]

where $\mathbf{\tilde{K}}=\mathbf{M}+\Delta t\mathbf{D}$. A Taylor expansion around the current guess gives,

\[\mathbf{R}\left( \mathbf{U}_{n+1}^{k}+\Delta \mathbf{U} \right)=\mathbf{R}\left( \mathbf{U}_{n+1}^{k} \right)+\mathbf{K}\left( \mathbf{U}_{n+1}^{k} \right)\Delta \mathbf{U}=\mathbf{0}\]

Here,

\[\begin{align}
  & \mathbf{K}=\frac{d\mathbf{R}}{d\mathbf{U}}=\mathbf{\tilde{K}}-\Delta t\frac{d\mathbf{F}}{d\mathbf{U}} \\ 
 & =\mathbf{M}+\Delta t\left( \mathbf{D}-\frac{d\mathbf{F}}{d\mathbf{U}} \right)  
\end{align}\]

This can be solved for $\Delta \mathbf{U}$.

\[\mathbf{K}\,\Delta \mathbf{U}=-\mathbf{R}\]

Then, the solution can be updated,

\[\mathbf{U}_{n+1}^{k+1}=\mathbf{U}_{n+1}^{k}+\Delta \mathbf{U}\]

### Generalized Trapezoidal Rule

Considering the general case,

\[\left( \mathbf{M}+\alpha \Delta t\mathbf{D} \right){{\mathbf{U}}_{n+1}}=\left( \mathbf{M}-\left( 1-\alpha  \right)\Delta t\mathbf{D} \right){{\mathbf{U}}_{n}}+\Delta t\left( \alpha \mathbf{F}\left( {{\mathbf{U}}_{n+1}} \right)+\left( 1-\alpha  \right)\mathbf{F}\left( {{\mathbf{U}}_{n}} \right) \right)\]

The residual in this case can be written as,

\[\mathbf{R}\left( {{\mathbf{U}}_{n+1}} \right)=\mathbf{M}\left( {{\mathbf{U}}_{n+1}}-{{\mathbf{U}}_{n}} \right)+\Delta t\mathbf{D}\left( \alpha {{\mathbf{U}}_{n+1}}+\left( 1-\alpha  \right){{\mathbf{U}}_{n}} \right)-\Delta t\left( \alpha {{\mathbf{F}}_{n+1}}+\left( 1-\alpha  \right){{\mathbf{F}}_{n}} \right)\]

It follows that the stiffness matrix is given by,

\[\mathbf{K}=\mathbf{M}+\alpha \Delta t\left( \mathbf{D}-\frac{d\mathbf{F}}{d\mathbf{U}} \right)\]

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

