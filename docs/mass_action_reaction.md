A mass-action reaction implements the classical model of a chemical reaction, converting reactants into products according to the law of mass action. 

For this reaction type, the reaction rate is given by:

\[
    r=k \prod_{i=1}^n(c_i)^{\nu_i}
\]

Here, $c_i$ are the concentrations of the reactants and $\nu_i$ is the stoichiometric coefficient of reactant $i$.

This type of reaction requires two parameters:

* `equation` : The reaction equation. 
* `rate_constant` : the rate constant $k$ in the equation above. 

The reaction equation is defined by writing the reaction formula as follows.

* Use the name of the species as defined in the Globals section to reference it.
* A stoichiometric coefficient of 1 does not need to be defined, but otherwise it must an integer number that appears before the species’ name and separated from it via a multiplication symbol `*`. 
* A plus symbol `+` is used to enumerate all reactants and products.
* The symbols `-->` or `->` can be used to separate the reactants from the products. 
* The reactant or product lists can be empty. 

Here are some examples. Assume that all the chemical species are properly defined in the `Globals` section. 

```
A-->B
A+B-->C
2*A+2*B->C+B
A-->
-->B
```
