# Black-Scholes
Numerical solution to the Black-Scholes PDE using finite difference methods FDM.


The original formulation of the Black-Scholes option pricing model was first published in the article: Robert C. Merton (1973) The Pricing of Options and Corporate Liabilities. The Journal of Political Economy.\\
The model consists of a partial differential equation (PDE) describing the dynamics of the efficient(no arbitrage) price $P(S(t), t)$ of a financial option on an underlying risky asset (ie. stock) with price $S(t)$, where $t$ is the time until the option contract expires.
The derivation of the PDE is based on a series of assumptions including:

* The price of the underlying asset is a Geometric Brownian motion process with constant drift (which we shall assume to be zero) and volatility (standard deviation of log-returns) $ \sigma $.
* In cases where the underlying asset is a stock, it does not pay dividends.
* No arbitrage: It is not possible to device a trading strategy based on combined trading of options and the underlying asset in a way that results in risk-free profits.
*  Market participants can borrow and lend any amount of cash and/or shares(even fractional) of stock(ie. short selling) at a risk-free rate denoted by $r$.


The PDE solved in this repo is the folloring 1D linear black scholes:

$$ u_t - \frac{1}{2} \sigma^2 x^2 u_{xx} - r x u_x + x u = 0 $$

With initial conditions:
In order to obtain solutions to the PDEs we will consider 3 different initial conditions:

European put:
$$u(x, 0) = (K-x)^+$$

Butterfly spread:
$$u(x, 0) = (x-K)^+ - 2(x-(K+H))^+ +(x+2H))^+$$

Binary call:
$$u(x, 0) = sgn^+(x - K)$$



Some results are:

Stable solution with butterfly initial condition and a forward euler numerical method:
![](https://github.com/erlendlokna/Numerical-solution-to-the-Black-Scholes-PDE/blob/main/figures/StableButterflySpreadFE.png)

Unstable forward euler solution:
![](https://github.com/erlendlokna/Numerical-solution-to-the-Black-Scholes-PDE/blob/main/figures/unstableFE.png)


Backwards euler and crank nicholson is also implemented, but results are not shown here.

