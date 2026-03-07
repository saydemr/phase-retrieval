# Explanation of Steps
1.  **Simulation (`S(x)` & `T(x)`)**: We modeled a sphere by calculating the chord length 
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
    - Regularization (`alpha`) was added to the denominator to handle frequencies 
      where `H(u)` approaches zero (zero-crossings), acting like a Wiener filter.
    - Finally, converted the recovered contact intensity `I_0` back to phase `phi`
      and then to thickness `S(x)` using the physical constants.

# Complications
- The propagation distance was given as 1cm in the initial description but later specified as 20cm. 
  I chose to implement the 20cm distance as explicitly requested in step 2.2. The paper explicitly talks about earlier methods requiring small propagation distances and ICT being able to work at larger distances, so I wanted to look at the longer propagation distance. Also, using 1cm gave very small variations in the retrieved thickness which was not interesting to me.

- The formula for converting phase to thickness was not explicitly given, so I found the formula (without the minus sign) in the literature. The negative sign was necessary to ensure the thickness came out positive. I could not find an explicit formula in Farago et al., it might be trivial for someone familiar with the topic.