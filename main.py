import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2, fftfreq

## Simulation Parameters
# Physical constants
WAVELENGTH = 1e-10          # wavelength (m)
WAVENUMBER = 2 * np.pi / WAVELENGTH
DELTA = 1e-7                # refractive index decrement
BETA = 1e-9                 # refractive index absorption part

# Geometry parameters
PIXEL_SIZE = 1e-6           # 1 micrometer (m)
GRID_SIZE = 512             # N x N pixels
PROP_DIST = 0.01


# 1. Simulate Transmission Function of a Sphere
def create_sphere_projection(grid_size, pixel_size, radius_pixels):
    """
    Creates the true projected thickness S(x) of a sphere on a 2D grid.
    S(x) = 2 * sqrt(R^2 - r^2)
    """
    x = np.linspace(-grid_size/2, grid_size/2, grid_size) * pixel_size
    y = np.linspace(-grid_size/2, grid_size/2, grid_size) * pixel_size
    X, Y = np.meshgrid(x, y)
    
    R = radius_pixels * pixel_size
    r_squared = X**2 + Y**2
    
    # Calculate thickness: 2 * sqrt(R^2 - r^2) inside the sphere, 0 outside
    thickness = np.zeros_like(r_squared)
    mask = r_squared < R**2
    thickness[mask] = 2 * np.sqrt(R**2 - r_squared[mask])
    
    return thickness

# Create the sphere
# Chose R = 100 pixels to be clearly visible within the 512x512 frame.
radius_px = 100
S_x = create_sphere_projection(GRID_SIZE, PIXEL_SIZE, radius_px)

# Compute Complex Transmission Function T(x)
# Formula from Assignment: T(x) = exp( -k * S(x) * [beta + i * delta])
# k = 2*pi/lambda.
# Argument is -k * beta * S(x) - i * k * delta * S(x)
T_x = np.exp(-WAVENUMBER * S_x * (BETA + 1j * DELTA))

# 2. Fresnel Propagation (Forward Model)
def fresnel_propagate(wavefront, dist, pixel_size, wavelength):
    """
    Propagates a wavefront using the Transfer Function method i.e., Fresnel Propagator.
    Includes padding to mitigate edge artifacts (circular convolution).
    """
    N = wavefront.shape[0]
    
    # Pad the wavefront to avoid wrapping artifacts
    # We pad to double the size
    pad_width = N // 2
    wavefront_padded = np.pad(wavefront, pad_width, mode='edge')
    N_pad = wavefront_padded.shape[0]
    
    # Spatial Frequencies
    u = fftfreq(N_pad, d=pixel_size)
    v = fftfreq(N_pad, d=pixel_size)
    U, V = np.meshgrid(u, v)
    
    # Squared spatial frequency magnitude |u|^2
    U_sq = U**2 + V**2
    
    # Fresnel Propagator in Fourier Domain: P(u) = exp(-i * pi * lambda * z * |u|^2)
    # The assignment asks for P_tilde(u) = exp(-i * pi * lambda * z * |u|^2)
    propagator = np.exp(-1j * np.pi * wavelength * dist * U_sq)
    
    # Perform Propagation: IFFT( FFT(T) * P )
    ft_wavefront = fft2(wavefront_padded)
    ft_prop = ft_wavefront * propagator
    wavefront_z = ifft2(ft_prop)
    
    # Crop back to original size
    start = pad_width
    end = start + N
    return wavefront_z[start:end, start:end]

# Propagate wavefield
U_z = fresnel_propagate(T_x, PROP_DIST, PIXEL_SIZE, WAVELENGTH)

# Compute Intensity I_z(x) = |U_z(x)|^2
I_z = np.abs(U_z)**2

# 3. Phase Retrieval (Inverse Model)
# Intensity Contrast Transfer (ICT) in Farago et al., Optics Letters 49, 18, 5159, Eq. 9 (adapted for single distance, i.e. M=1)
def ict_phase_retrieval(intensity_z, dist, pixel_size, wavelength, delta, beta):
    """
    Retrieves the object phase/thickness using the ICT formula.
    """
    N = intensity_z.shape[0]
    
    # Spatial frequencies
    u = fftfreq(N, d=pixel_size)
    v = fftfreq(N, d=pixel_size)
    U, V = np.meshgrid(u, v)
    U_sq = U**2 + V**2
    
    # 1. Calculate the Contrast Transfer Function H(u)
    # H(u) = cos(pi * lambda * z * |u|^2) + (delta / beta) * sin(pi * lambda * z * |u|^2)
    phase_factor = np.pi * wavelength * dist * U_sq
    H_u = np.cos(phase_factor) + (delta / beta) * np.sin(phase_factor)
    
    # 2. Fourier Transform of the measured Intensity
    I_z_hat = fft2(intensity_z)

    # 3. Frequency-Dependent Regularization
    # The first maximum of H(u) occurs when tan(A) = delta/beta.
    # We find the spatial frequency squared |u|^2 corresponding to this peak.
    # see Readme for detailed derivation.
    # np.arctan returns values in (-pi/2, pi/2). Since delta / beta > 0, result is positive.
    u_sq_peak = np.arctan(delta / beta) / (np.pi * wavelength * dist)
    
    # Create the alpha map
    # alpha = 0 for frequencies below the first peak 
    # alpha = 1 for frequencies above the first peak 
    # no square root needed since we are comparing U_sq directly to u_sq_peak
    alpha_map = np.where(U_sq <= u_sq_peak, 0.0, 1.0)

    # 4. Inversion (Regularized Division)
    # We recover FFT{I_0} using Eq (9) with M=1
    # F{I_0} = (H(u) * F{I_z}) / (H(u)^2 + alpha)
    numerator = H_u * I_z_hat
    denominator = H_u**2 + alpha_map
    I_0_hat = numerator / denominator
    
    # 5. Inverse FFT to get retrieved contact intensity I_0
    I_0_rec = np.real(ifft2(I_0_hat))
    
    # Clip to avoid log of negative numbers
    I_0_rec = np.clip(I_0_rec, 1e-10, None)
    
    # 6. Retrieve Phase phi(x)
    # Formula Eq (8)-(9): phi(x) = (delta / 2*beta) * ln(I_0_rec)
    phi_rec = (delta / (2 * beta)) * np.log(I_0_rec)
    
    # 7. Convert Phase to Thickness S(x)
    # phi(x) = -k * delta * S(x)
    # So: S(x) = -phi(x) / (k * delta)
    k = 2 * np.pi / wavelength
    S_rec = -phi_rec / (k * delta)
    
    return S_rec

# Perform retrieval
S_reconstructed = ict_phase_retrieval(I_z, PROP_DIST, PIXEL_SIZE, WAVELENGTH, DELTA, BETA)


## 4. Visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: True Thickness S(x)
im1 = axes[0].imshow(S_x * 1e6, cmap='gray', extent=[-256, 256, -256, 256])
axes[0].set_title(r"True Projected Thickness S(x) [$\mu$m]")
axes[0].set_xlabel(r"x [$\mu$m]")
plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)

# Plot 2: Simulated Intensity I_z(x)
im2 = axes[1].imshow(I_z, cmap='gray', extent=[-256, 256, -256, 256])
axes[1].set_title(f"Simulated Intensity @ z={PROP_DIST*100:.0f}cm")
axes[1].set_xlabel(r"x [$\mu$m]")
plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)

# Plot 3: Reconstructed Thickness
im3 = axes[2].imshow(S_reconstructed * 1e6, cmap='gray', extent=[-256, 256, -256, 256])
axes[2].set_title(r"Reconstructed Thickness [$\mu$m] (ICT)")
axes[2].set_xlabel(r"x [$\mu$m]")
plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()

# put cross-section plots below to a figure with 2 rows and 1 column

fig, axes = plt.subplots(2, 1, figsize=(10, 10))


# Cross-section comparison
mid = GRID_SIZE // 2
axes[0].plot(S_x[mid, :] * 1e6, 'k--', label='True Thickness', linewidth=2)
axes[0].plot(S_reconstructed[mid, :] * 1e6, 'r-', label='Reconstructed (ICT)', alpha=0.7)
axes[0].set_title(f"Cross-section Profile Comparison, z={PROP_DIST*100:.0f}cm")
axes[0].set_xlabel("Pixel Index")
axes[0].set_ylabel(r"Thickness [$\mu$m]")

# cross section comparison, around the center of the sphere, with 16 pixel width
axes[1].plot(range(240,272), S_x[mid, 240:272] * 1e6, 'k--', label='True Thickness', linewidth=2)
axes[1].plot(range(240,272), S_reconstructed[mid, 240:272] * 1e6, 'r-', label='Reconstructed (ICT)', alpha=0.7)
axes[1].set_title(f"Zoomed Cross-section Profile Comparison, z={PROP_DIST*100:.0f}cm")
axes[1].set_xlabel("Pixel Index")
axes[1].set_ylabel(r"Thickness [$\mu$m]")

plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

# Phase shift vs position plot, through y = 0 line
phi_rec = (2 * np.pi / WAVELENGTH) * DELTA * S_reconstructed
plt.figure(figsize=(10, 5))
plt.plot(phi_rec[mid, :], 'b-', label='Retrieved Phase Shift (rad)', alpha=0.7)
plt.title(f"Retrieved Phase Shift Profile, z={PROP_DIST*100:.0f}cm, y=0 line")
plt.xlabel("Pixel Index")
plt.ylabel("Phase Shift (radians)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()


print(f"Max True Thickness: {np.max(S_x)*1e6:.4f} um")
print(f"Max Rec Thickness:  {np.max(S_reconstructed)*1e6:.4f} um")