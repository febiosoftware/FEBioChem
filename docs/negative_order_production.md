The `negative-order production` reaction is a production reaction whose rate depends inversely on the concentration of the product species. Consequently, the production rate will decrease as the concentration of the product species increases. Negative-order production kinetics can be used to model self-limiting production processes in which accumulation of a species reduces its own production rate.

The reaction equation is,

\[
    \emptyset \to P
\]

The reaction rate is defined by,

\[
    r = r_0 + \frac{k}{\left( c_P + \epsilon \right)^{n}}
\]

Here, $k$ is the rate constant. The rate offset $r_0$ provides a constant contribution to the production rate and may be positive or negative. The concentration offset $\epsilon$ prevents the reaction rate from becoming singular as the concentration approaches zero. The order $n$ may be any non-negative real number. Setting $n=0$ reduces the reaction to a constant-rate production process with rate $r=r_0+k$.

The following table lists the input parameters for this type of reaction. 

| Name                   | Description                      | Default | Range            |
|------------------------|----------------------------------|---------|------------------|
| `rate_offset`          | rate offset $r_0$                |       0 | $\in \mathbb{R}$ |
| `rate_constant`        | rate constant $k$                |       0 | $\ge 0$          |
| `concentration_offset` | concentration offset ($\epsilon$)|       0 | $\ge 0$          |
| `order`                | reaction order $n$.              |       1 | $\ge 0$          |
| `product`              | The name of the product species  |         |                  |
