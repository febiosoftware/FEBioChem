The `negative-order reaction` reaction is a reaction whose rate depends inversely on the concentration of the reactant. For such reactions, the reaction rate will increase as the concentration of the reactant decreases. The reaction is represented as a conversion of a reactant $R$ into a product $P$,

\[
    R \to P
\]

The reaction rate is defined by,

\[
    r = \frac{k}{\left( c_R + \epsilon \right)^{n}}
\]

Here, $k$ is the rate coefficient, $\epsilon$ is the concentration offset and $n$ is the order of the reaction. The concentration offset $\epsilon$ prevents the reaction rate from becoming singular as the reactant concentration approaches zero. The reaction order $n$ may be any non-negative real number."

The following table lists the input parameters for this type of reaction. 

| Name     | Description                      | Default | Range   |
|----------|----------------------------------|---------|---------|
| `k`      | rate coefficient $k$             |       0 | $\ge 0$ |
| `offset` | concentration offset ($\epsilon$)|       0 | $\ge 0$ |
| `n`      | reaction order $n$.              |       1 | $\ge 0$ |

