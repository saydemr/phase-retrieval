# ICT Phase Retrieval Implementation

Implementation of the functions in [*"Phase retrieval in propagation-based X-ray imaging beyond the limits of transport of intensity and contrast transfer function approaches"*](https://doi.org/10.1364/OL.530330) by Faragó et. al. (2024). The code simulates a sphere and retrieves its thickness from the simulated intensity using the ICT method.  


# Explanation of Steps
1.  **Simulation (`S(x)` & `T(x)`)**: Modeled a sphere by calculating the chord length 
    at each pixel. The complex transmission `T(x)` combines attenuation (beta) and 
    phase shift (delta).

2.  **Propagation (Fresnel)**: Used the Fourier Transform approach. 
    The propagator `exp(-i*pi*lambda*z*u^2)` describes how spatial frequencies 
    evolve over distance `z`. Padding was applied to prevent edge artifacts caused 
    by the periodic nature of the FFT.

3.  **Phase Retrieval (ICT)**:
    Implemented the Single-Distance ICT formula from the referenced paper.
    - Computed the Contrast Transfer Function (CTF) `H(u)`.
    - Inverted the formation model in Fourier space: `F{I_0} = F{I_z} / H(u)`.
    - Conditional regularization (`alpha_map`) was added to the denominator. **Derivation is provided below**.
    - Finally, converted the recovered contact intensity `I_0` back to phase `phi`
      and then to thickness `S(x)` using the physical constants.


# Regularization of the CTF Function
In the referenced paper, authors state that they used a conditional regularization approach. 
They use $\alpha = 0$ for frequencies up to first maximum and $\alpha = 1$ for frequencies above the first maximum, of $H(u)$.
In the first version of the code, I implemented a simple regularization with a constant $\alpha$ value of 1e-2, assuming that it is a simple constant to avoid division by zero. Turns out that this is not the case, and the paper explicitly states that they used a conditional regularization approach later in the simulation results paragraph.


We have:

$$
H(u) = \cos(\pi\lambda z|u|^{2}) + \frac{\delta}{\beta}\sin(\pi\lambda z|u|^{2})
$$

Then,

$$
\frac{dH}{dA} = -\sin(A) + \frac{\delta}{\beta}\cos(A) = 0, \quad \text{where } A = \pi\lambda z|u|^{2}
$$

Rearranging the terms gives the condition for the maxima (first extremum is maxima, see [https://www.desmos.com/calculator/9uw9yhy3fm](https://www.desmos.com/calculator/9uw9yhy3fm)):

$$
\sin(A) = \frac{\delta}{\beta}\cos(A) \implies \tan(A) = \frac{\delta}{\beta}
$$

Then, we can solve for the spatial frequency at the first peak:

$$
\pi\lambda z |u|_{peak}^2 = \arctan\left(\frac{\delta}{\beta}\right)
$$

$$
|u|_{peak} = \sqrt{\frac{1}{\pi\lambda z} \arctan\left(\frac{\delta}{\beta}\right)}
$$


# Phase to Thickness Conversion
Given the transmission function in the assignment:

$$
T(x) = \exp\left\lbrace-\frac{2\pi}{\lambda} S(x) [\beta + i\delta]\right\rbrace
$$

and equation (1) from the paper:

$$
T(x) = \exp\left\lbrace-B(x) + i\varphi(x)\right\rbrace
$$


The imaginary parts must be equal.

$$
\exp\left\lbrace i\varphi(x)\right\rbrace = \exp\left\lbrace -i \frac{2\pi}{\lambda}\delta S(x)\right\rbrace
$$

Equating the arguments of the exponential functions:

$$
i\varphi(x) = -i \frac{2\pi}{\lambda}\delta S(x)
$$

$$
\varphi(x) = -\frac{2\pi}{\lambda}\delta S(x)
$$

Letting the wavenumber $k = \frac{2\pi}{\lambda}$ and isolating $S(x)$ gives:

$$
S(x) = -\frac{\varphi(x)}{k \delta}
$$

with the minus sign in the beginning.


# Running the Code
pyproject.toml contains the dependencies. You can install them using uv:
```bash
uv install
```
To run the code, execute:
```bash
uv run main.py
```
This will generate the plots for the simulated intensity and the retrieved thickness.
If you do not have uv, you can install it from [https://docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/).
