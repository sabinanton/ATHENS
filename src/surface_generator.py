import numpy as np
import sys
import scipy as sp
import scipy.special as spc
import matplotlib.pyplot as plt
import scipy.ndimage as nd
import scipy.signal as sn
from stl import mesh

class Surface:
    """
    A class used to represent a 3D surface with stochastic characteristics.

    The surface is generated based on mu and sigma arrays, which control the shape and
    roughness of the surface. This class also handles saving the surface geometry to a file.
    """

    def __init__(self, mu_array, sigma_array, sample_length, local_parameters):
        """
        Initializes the Surface object.

        Parameters:
        - mu_array: Array of coefficients for the mu parameter of the surface.
        - sigma_array: Array of coefficients for the sigma parameter of the surface.
        - sample_length: Length of the surface to be generated.
        - local_parameters: Additional local parameters associated with the surface.
        """
        self.mu_array = mu_array
        self.sigma_array = sigma_array
        self.T = 1  # Autocorrelation length
        self.sample_length = sample_length
        self.param_local = local_parameters
        self.samples = None
        self.X = None
        self.Y = None
        self.Z = None
        self.mesh = None
        self.sigma_T = None

    def find_coefficients(self, func, n_coeff, lim, N):
        """
        Computes the Hermite polynomial coefficients for a given function.

        Parameters:
        - func: The function to fit.
        - n_coeff: The number of Hermite coefficients to compute.
        - lim: The range over which to compute the coefficients.
        - N: The number of points to use for the fitting.

        Returns:
        - Hermite coefficients as a numpy array.
        """
        x = np.linspace(-lim, lim, N)
        return np.polynomial.hermite.hermfit(x, func(x), n_coeff)

    def compute_Hermite_expansion(self, coeff_array, tau):
        """
        Evaluates the Hermite expansion for the given coefficients and input values.

        Parameters:
        - coeff_array: Array of Hermite coefficients.
        - tau: Input value or array for evaluating the expansion.

        Returns:
        - The evaluated Hermite expansion.
        """
        return np.polynomial.Hermite(coeff_array)(tau)

    def compute_psd(self, Z):
        """
        Computes the power spectral density (PSD) of the surface Z.

        Parameters:
        - Z: The surface height matrix.

        Returns:
        - The power spectral density as a numpy array.
        """
        return np.fft.fftshift(np.abs(np.fft.fft2(np.fft.ifftshift(Z))) ** 2)

    def generate_height_pdf(self, Z):
        """
        Generates the height distribution function for the surface based on mu and sigma arrays.

        Parameters:
        - Z: The surface height matrix.

        Returns:
        - The modified surface height matrix with a Gaussian height distribution.
        """
        indices = np.random.randn(Z.shape[0] * Z.shape[1])
        mus = self.compute_Hermite_expansion(self.mu_array, indices)
        sigmas = self.compute_Hermite_expansion(self.sigma_array, indices)
        Z_flat = np.random.randn(Z.shape[0] * Z.shape[1]) * sigmas + mus
        Z = Z_flat.reshape(Z.shape)

        return Z

    def correct_PSD(self, Z, target_psd):
        """
        Adjusts the PSD of the surface to match the target PSD.

        Parameters:
        - Z: The surface height matrix.
        - target_psd: The target power spectral density to match.

        Returns:
        - The modified surface height matrix with the corrected PSD.
        """
        Z_fft = (np.fft.fft2(Z))
        Z_new_fft = target_psd * np.exp(1j * np.angle(Z_fft))
        Z_new = np.real(np.fft.ifft2((Z_new_fft)))

        return Z_new

    def correct_height(self, Z_h, Z_p):
        """
        Corrects the height of the surface based on a reference surface.

        Parameters:
        - Z_h: The target height distribution surface.
        - Z_p: The reference surface.

        Returns:
        - The modified surface with corrected heights.
        """
        index_p = np.argsort(Z_p.flatten())
        Z_p = Z_p.flatten()
        Z_p[index_p] = Z_h.flatten()

        return Z_p.reshape((self.samples, self.samples))

    def generate(self, samples, iterations, verbose):
        """
        Generates the surface by iteratively correcting its PSD and height distribution.

        Parameters:
        - samples: Number of samples along each dimension of the surface.
        - iterations: Number of iterations for the PSD and height correction.
        - verbose: Boolean flag to print details during the surface generation.
        """
        self.samples = samples
        x = np.linspace(-self.sample_length / 2, self.sample_length / 2, samples)
        y = np.linspace(-self.sample_length / 2, self.sample_length / 2, samples)
        self.X, self.Y = np.meshgrid(x, y)
        self.Z = np.zeros_like(self.X)
        gamma = np.random.randn(samples ** 2).reshape(self.Z.shape)
        epsilon = np.random.randn(samples ** 2).reshape(self.Z.shape)
        epsilon_order = np.sort(epsilon.flatten()).reshape(epsilon.shape)
        gamma_order = np.sort(gamma.flatten()).reshape(epsilon.shape)
        psd_target = (np.abs(np.fft.fft2(np.exp(- np.abs(self.X ** 2 + self.Y ** 2) / self.T ** 2))))

        for i in range(iterations):
            epsilon_p = self.correct_PSD(epsilon, psd_target)
            epsilon = self.correct_height(epsilon_order, epsilon_p)
            gamma_p = self.correct_PSD(gamma, psd_target)
            gamma = self.correct_height(gamma_order, gamma_p)

        self.epsilon = epsilon
        self.gamma = gamma
        self.Z = np.polynomial.Hermite(self.sigma_array)(gamma) * epsilon + np.polynomial.Hermite(self.mu_array)(gamma)
        self.sigma_T = np.var(self.Z) ** 0.5 / self.compute_autocorr_length(self.Z)
        if verbose:
            print("Generated surface with roughness parameters", np.var(self.epsilon) ** 0.5, self.compute_autocorr_length(self.epsilon))

        self.mesh = np.empty((self.samples - 1, self.samples - 1, 2, 3), dtype=int)

        r = np.arange(self.samples * self.samples).reshape(self.samples, self.samples)

        self.mesh[:, :, 0, 0] = r[:-1, :-1]
        self.mesh[:, :, 1, 0] = r[:-1, 1:]
        self.mesh[:, :, 0, 1] = r[:-1, 1:]

        self.mesh[:, :, 1, 1] = r[1:, 1:]
        self.mesh[:, :, :, 2] = r[1:, :-1, None]

        self.mesh.shape = (-1, 3)

        self.T = self.compute_autocorr_length(self.epsilon)

    def compute_autocorr_length(self, Z):
        """
        Computes the autocorrelation length of the surface.

        Parameters:
        - Z: The surface height matrix.

        Returns:
        - The autocorrelation length as a float.
        """
        Z_mean = np.mean(Z)
        centered_Z = Z - Z_mean
        fft_surface = np.fft.fft2(centered_Z)
        power_spectrum = np.abs(fft_surface) ** 2
        autocorr = np.fft.ifft2(power_spectrum)
        autocorr = np.fft.fftshift(autocorr)
        autocorr = autocorr / np.max(autocorr)
        tau = (self.X ** 2 + self.Y ** 2) ** 0.5
        T = np.max(tau[autocorr >= np.exp(-1)])

        return T

    def save(self, project_name, file_name):
        """
        Saves the generated surface to SRF files.

        Parameters:
        - project_name: The name of the project directory.
        - file_name: The base name of the files to save.
        """

        f = open(project_name + "/surface_files/" + file_name + ".srf", "w+")
        f.write(str(len(self.mu_array)) + "\n")
        for i in range(len(self.mu_array)):
            f.write(str(self.mu_array[i]) + " ")
        f.write("\n")
        for i in range(len(self.sigma_array)):
            f.write(str(self.sigma_array[i]) + " ")
        f.write("\n")
        f.write(str(self.T) + "\n")
        f.write(str(len(self.param_local)) + "\n")
        for i in range(len(self.param_local)):
            f.write(str(self.param_local[i]) + " ")
        f.write("\n")
        f.close()


def compute_average_autocorrelation_function(mu_k, sigma_k, num_surfaces, num_samples, sample_length, iterations):
    """
    Computes the average autocorrelation function across multiple surfaces.

    Parameters:
    - mu_k: Hermite coefficients for the mu parameter.
    - sigma_k: Hermite coefficients for the sigma parameter.
    - num_surfaces: The number of surfaces to generate for averaging.
    - num_samples: The number of samples along each axis for each surface.
    - sample_length: The length of the surface to be generated.
    - iterations: The number of iterations for PSD and height correction.

    Returns:
    - The average autocorrelation length across the generated surfaces.
    """
    T = 0.0

    for i in range(num_surfaces):
        # Create a Surface object and generate the surface
        surf = Surface(mu_k, sigma_k, sample_length, [0, 0, 0])
        surf.generate(num_samples, iterations, verbose=False)
        T_s = surf.compute_autocorr_length(surf.epsilon)
        T += T_s

    return T / num_surfaces


def mu(gamma):
    """
    Defines the mu function as a function of gamma.

    Parameters:
    - gamma: Input values for the mu function.

    Returns:
    - The mu value based on the input gamma values.
    """
    return spc.erf(gamma * 2.0 - 1.0) * 0.0


def sigma(gamma):
    """
    Defines the sigma function as a function of gamma.

    Parameters:
    - gamma: Input values for the sigma function.

    Returns:
    - The sigma value based on the input gamma values.
    """
    return 0.8 + (1 + spc.erf(gamma * 2.0 - 1.0)) * 0.0


# Margin for gamma values in the Hermite fit
margin = 8.0
# Number of Hermite coefficients
num_order = 40
# Number of samples along each axis for the surface
num_samples = 100
# Length of the surface
sample_length = 200

# Generate gamma values for the Hermite fitting
gamma_fit = np.linspace(-margin, margin, 1000)

# Compute the Hermite coefficients for mu and sigma using a polynomial Hermite fit
mu_k = np.polynomial.hermite.hermfit(gamma_fit, mu(gamma_fit), num_order)
sigma_k = np.polynomial.hermite.hermfit(gamma_fit, sigma(gamma_fit), num_order)

# Create a Surface object with the fitted Hermite coefficients
surf = Surface(mu_k, sigma_k, sample_length, [0.0, 0.0, 0.0])

# Generate the surface using the specified number of samples and iterations
surf.generate(num_samples, iterations=400, verbose=False)

# Compute the average autocorrelation length across 50 surfaces
surf.T = compute_average_autocorrelation_function(mu_k, sigma_k, 50, 50, 10, 100)

# Save the generated surface to a file
surf.save("Plane_test_rough", "plane")
