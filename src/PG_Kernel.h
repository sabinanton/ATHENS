
/**
 * @brief Represents a class for the proposed rough surface kernel applied to Poly-Gaussian surfaces.
 * 
 * This class implements the rough GSI kernel algorithm for trajectory sampling.
 * It returns the trajectories of free-molecular-flow gas particles on a rough surface.
 * 
 * @tparam T The type of the elements in the PG_Kernel.
 */
template <typename T>
class PG_Kernel {

    private:

        Surface<T> surface; // Surface object for storing surface properties
        Gas<T> gas; // Gas object for storing gas properties

        Matrix<T> roots; // Matrix of roots for Hermite polynomials
        Matrix<T> weights; // Matrix of weights for Hermite polynomials
        Matrix<T> roots_2d_r; // Matrix of roots for 2D Hermite polynomials (row)
        Matrix<T> roots_2d_c; // Matrix of roots for 2D Hermite polynomials (column)
        Matrix<T> weights_2d; // Matrix of weights for 2D Hermite polynomials
        Matrix<T> gamma; // Matrix of gamma values for PG_Kernel
        Matrix<T> gamma_dot; // Matrix of gamma_dot values for PG_Kernel
        Matrix<T> gamma_2d; // Matrix of gamma values for 2D PG_Kernel
        Matrix<T> gamma_dot_2d; // Matrix of gamma_dot values for 2D PG_Kernel
        Matrix<T> p_gamma; // Matrix of p_gamma values for PG_Kernel
        Matrix<T> p_gamma_dot; // Matrix of p_gamma_dot values for PG_Kernel
        Matrix<T> p_gamma_2d; // Matrix of p_gamma values for 2D PG_Kernel
        Matrix<T> p_gamma_dot_2d; // Matrix of p_gamma_dot values for 2D PG_Kernel
        Matrix<T> mu; // Matrix of mu values for PG_Kernel
        Matrix<T> sigma; // Matrix of sigma values for PG_Kernel
        Matrix<T> mu_2d; // Matrix of mu values for 2D PG_Kernel
        Matrix<T> sigma_2d; // Matrix of sigma values for 2D PG_Kernel
        Matrix<T> mu_gamma; // Matrix of mu_gamma values for PG_Kernel
        Matrix<T> sigma_gamma; // Matrix of sigma_gamma values for PG_Kernel
        Matrix<T> mu_gamma_2d; // Matrix of mu_gamma values for 2D PG_Kernel
        Matrix<T> sigma_gamma_2d; // Matrix of sigma_gamma values for 2D PG_Kernel

        unsigned int int_order; // Order of integration for Hermite polynomials

        std::string sim_name; // Name of the simulation

        T margin, ac_length, mu_gas, T_gas; // Parameters for PG_Kernel

        void (*local_kernel)(T*, T*, T*, const Gas<T>&, Surface<T>&); // Pointer to local kernel function

        Hermite<T> herm_tools; // Hermite tools object


    public:

        /**
         * @brief Default constructor for PG_Kernel.
         * 
         * Initializes the PG_Kernel object with default values.
         */
        PG_Kernel() {

            Hermite<T> herm_tools;

            int_order = 2;
            margin = 5.0;
            ac_length = 1.0;
            mu_gas = 4.0;
            T_gas = 400.0;
            
        }

        /**
         * @brief Constructor for PG_Kernel.
         * 
         * Initializes the PG_Kernel object with the given parameters.
         * 
         * @param surface The Surface object representing the surface properties.
         * @param gas The Gas object representing the gas properties.
         * @param local_kernel Pointer to the local kernel function.
         * @param incident_velocity Matrix of incident velocities.
         * @param num_particles Number of particles.
         * @param sim_name Name of the simulation.
         */
        PG_Kernel(Surface<T>& surface, const Gas<T>& gas, void (*local_kernel)(T*, T*, T*, const Gas<T>& , Surface<T>&)) {

            int_order = 2;
            margin = 5.0;

            this->surface = surface;
            this->gas = gas;
            this->ac_length = surface.autocorrelation_length();

            roots = herm_tools.get_roots(int_order);
            weights = herm_tools.get_weights(int_order);

            std::tuple<Matrix<T>, Matrix<T>> roots_tuple = roots.meshgrid(roots);
            roots_2d_r = std::get<1>(roots_tuple);
            roots_2d_c = std::get<0>(roots_tuple);
            weights_2d = weights.outer(weights);

            gamma = roots;
            gamma_dot = roots;
            gamma_2d = roots_2d_r;
            gamma_dot_2d = roots_2d_c;
            mu = herm_tools.evaluate(surface.get_mu_coefficients(), gamma);
            sigma = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma);
            mu_2d = herm_tools.evaluate(surface.get_mu_coefficients(), gamma_2d);
            sigma_2d = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma_2d);
            mu_gamma = herm_tools.evaluate_deriv(surface.get_mu_coefficients(), gamma);
            sigma_gamma = herm_tools.evaluate_deriv(surface.get_sigma_coefficients(), gamma);
            mu_gamma_2d = herm_tools.evaluate_deriv(surface.get_mu_coefficients(), gamma_2d);
            sigma_gamma_2d = herm_tools.evaluate_deriv(surface.get_sigma_coefficients(), gamma_2d);

            p_gamma = exp(gamma * gamma * (-1.0) / 2.0) * 1.0 / std::sqrt(2.0 * M_PI);
            p_gamma_2d = exp(gamma_2d * gamma_2d * (-1.0) / 2.0) * 1.0 / std::sqrt(2.0 * M_PI);
            p_gamma_dot = exp(gamma_dot * gamma_dot * (-1.0) * ac_length * ac_length / 4.0) * ac_length / std::sqrt(4.0 * M_PI);
            p_gamma_dot_2d = exp(gamma_dot_2d * gamma_dot_2d * (-1.0) * ac_length * ac_length / 4.0) * ac_length / std::sqrt(4.0 * M_PI);

            mu_gas = gas.get_molar_mass();
            T_gas = gas.get_temperature();

            this->local_kernel = local_kernel;
            
        }

        /**
         * @brief Destructor for PG_Kernel.
         * 
         * Cleans up the PG_Kernel object.
         */
        ~PG_Kernel() {

            int_order = 0;
            margin = 0.0;
            ac_length = 0.0;
            mu_gas = 0.0;
            T_gas = 0.0;
            this->local_kernel = nullptr;

        }  

        /**
         * @brief Get the Gas object.
         * 
         * @return Gas<T>& The Gas object.
         */
        Gas<T>& get_gas() {return this->gas;}

        /**
         * @brief Get the Surface object.
         * 
         * @return Surface<T>& The Surface object.
         */
        Surface<T>& get_surface() {return this->surface;}

        /**
         * @brief Perform random sampling in 1D using the Metropolis-Hastings algorithm.
         * 
         * @param pdf Pointer to the PDF function.
         * @param parameters Vector of parameters for the PDF function.
         * @param x0 Initial position.
         * @param size Size of the sampling range.
         * @param domain_x Array of domain limits.
         * @return T The sampled value.
         */
        T random_sampler_1D(T (PG_Kernel::*pdf)(T, std::vector<T>), std::vector<T> parameters, T x0, T size, T domain_x[]);

        /**
         * @brief Perform random sampling in 2D using the Metropolis-Hastings algorithm.
         * 
         * @param pdf Pointer to the PDF function.
         * @param parameters Vector of parameters for the PDF function.
         * @param x0 Initial x position.
         * @param y0 Initial y position.
         * @param size_x Size of the x sampling range.
         * @param size_y Size of the y sampling range.
         * @param domain_x Array of x domain limits.
         * @param domain_y Array of y domain limits.
         * @return std::tuple<T, T> The sampled x and y values.
         */
        std::tuple<T, T> random_sampler_2D(T (PG_Kernel::*pdf)(T, T, std::vector<T>), std::vector<T> parameters, T x0, T y0, T size_x, T size_y, T domain_x[], T domain_y[]);

        /**
         * @brief Perform random sampling of angles in 2D using the Metropolis-Hastings algorithm.
         * 
         * @param pdf Pointer to the PDF function.
         * @param parameters Vector of parameters for the PDF function.
         * @param x0 Initial x position.
         * @param y0 Initial y position.
         * @param size_x Size of the x sampling range.
         * @param size_y Size of the y sampling range.
         * @param domain_x Array of x domain limits.
         * @param domain_y Array of y domain limits.
         * @return std::tuple<T, T> The sampled theta_r1 and theta_r2 values.
         */
        std::tuple<T, T> random_sampler_angles(T (PG_Kernel::*pdf)(T, T, std::vector<T>), std::vector<T> parameters, T x0, T y0, T size_x, T size_y, T domain_x[], T domain_y[]);

        /**
         * @brief Calculate the PDF of the control process gamma and its derivative gamma_dot, for Poly-Gaussian surfaces.
         * 
         * @param gamma The gamma value.
         * @param gamma_dot The gamma_dot value.
         * @param parameters Vector of parameters for the PDF function.
         * @return T The PDF value.
         */
        T mixture_pdf(T gamma, T gamma_dot, std::vector<T> parameters); 

        /**
         * @brief Calculate the reflected particle angular PDF according to the Kirchhoff model.
         * 
         * @param theta_r1 The theta_r1 angle value.
         * @param theta_r2 The theta_r2 angle value.
         * @param parameters Auxiliary parameters required: incidence angle and control process values.
         * @return T The PDF value.
         */
        T kirchhoff_pdf(T theta_r1, T theta_r2, std::vector<T> parameters);

        /**
         * @brief Calculate the probability of a reflected particle to escape the surface.
         * 
         * @param theta_r1 The theta_r1 angle.
         * @param height The non-dimensional height of the particle.
         * @param direction The direction of the particle trajectory.
         * @return T The value of the probability.
         */
        T shadow(T theta_r1, T height, int direction);

        /**
         * @brief Calculate the probability of a reflected particle to escape the surface, given a certain travelled horizontal distance.
         * 
         * @param theta_r1 The theta_r1 angle.
         * @param height The non-dimensional height of the particle.
         * @param distance The horizontal distance value.
         * @param direction The direction of the particle trajectory.
         * @return T The value of the probability.
         */
        T occlusion(T theta_r1, T height, T distance, int direction);

        /**
         * @brief Calculate the PDF for the height_0 parameter.
         * 
         * @param height_0 The height_0 value.
         * @param parameters Vector of parameters for the PDF function.
         * @return T The PDF value.
         */
        T height_0_pdf(T height_0, std::vector<T> parameters);

        /**
         * @brief Calculate the PDF for the distance parameter.
         * 
         * @param distance The distance value.
         * @param parameters Vector of parameters for the PDF function.
         * @return T The PDF value.
         */
        T distance_pdf(T distance, std::vector<T> parameters);

        /**
         * @brief Convert global coordinates to local coordinates.
         * 
         * @param global_vector The global vector.
         * @param theta_i The theta_i value.
         * @param theta_r1 The theta_r1 value.
         * @param theta_r2 The theta_r2 value.
         * @return Matrix<T> The local vector.
         */
        Matrix<T> global_to_local(Matrix<T> global_vector, T theta_i, T theta_r1, T theta_r2);

        /**
         * @brief Convert local coordinates to global coordinates.
         * 
         * @param local_vector The local vector.
         * @param theta_i The theta_i value.
         * @param theta_r1 The theta_r1 value.
         * @param theta_r2 The theta_r2 value.
         * @return Matrix<T> The global vector.
         */
        Matrix<T> local_to_global(Matrix<T> local_vector, T theta_i, T theta_r1, T theta_r2);

        /**
         * @brief Convert angles to slopes.
         * 
         * @param theta_i The theta_i value.
         * @param theta_r1 The theta_r1 value.
         * @param theta_r2 The theta_r2 value.
         * @return Matrix<T> The slopes.
         */
        Matrix<T> angles_to_slopes(T theta_i, T theta_r1, T theta_r2);

        /**
         * @brief Sample a single particle trajectory.
         * 
         * @return trajectory<T> The sampled trajectory.
         */
        Matrix<T> sample_particle(Matrix<T> incident_velocity, T* gsi_parameters);
        

};

template <typename T>
T PG_Kernel<T>::random_sampler_1D(T (PG_Kernel::*pdf)(T, std::vector<T>), std::vector<T> parameters, T x0, T size, T domain_x[]) {

    unsigned long steps = 50;

    T x = x0, x_p = 0.0, acceptance = 0.0;

    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 rng(fullSeed);
    std::uniform_real_distribution<T> uniformDist(0.0f, 1.0f);
    std::normal_distribution<T> normDist(0.0f, 1.0f);

    for(auto i = 0; i < steps; i ++) {

        T step = normDist(rng) * size;
        T u_sample = uniformDist(rng);
        x_p = x + step;
        acceptance = 0.0;
        if(x_p >= domain_x[0] && x_p <= domain_x[1]) {
            acceptance = std::min(1.0, (this->*pdf)(x_p, parameters) / (this->*pdf)(x, parameters));
        }
        if(acceptance > u_sample) {
            x = x_p;
        }
    }

    return x;
}

template <typename T>
std::tuple<T, T> PG_Kernel<T>::random_sampler_2D(T (PG_Kernel::*pdf)(T, T, std::vector<T>), std::vector<T> parameters, T x0, T y0, T size_x, T size_y, T domain_x[], T domain_y[]) {

    unsigned long steps = 50;

    T x = x0, x_p, y_p;
    T y = y0;
    T acceptance = 0.0;

    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 rng(fullSeed);
    std::uniform_real_distribution<T> uniformDist(0.0f, 1.0f);
    std::normal_distribution<T> normDist(0.0f, 1.0f);

    for(auto i = 0; i < steps; i ++) {

        T step_x = normDist(rng) * size_x;
        T step_y = normDist(rng) * size_y;
        T u_sample = uniformDist(rng);
        x_p = x + step_x;
        y_p = y + step_y;
        acceptance = 0.0;
        if(x_p >= domain_x[0] && x_p <= domain_x[1] && y_p >= domain_y[0] && y_p <= domain_y[1]) {
            acceptance = std::min(1.0, (this->*pdf)(x_p, y_p, parameters) / (this->*pdf)(x, y, parameters));
        }
        if(acceptance > u_sample) {
            x = x_p;
            y = y_p;
        }
    }

    return std::make_tuple(x, y);
}

template <typename T>
std::tuple<T, T> PG_Kernel<T>::random_sampler_angles(T (PG_Kernel::*pdf)(T, T, std::vector<T>), std::vector<T> parameters, T theta_r1, T theta_r2, T size_x, T size_y, T domain_x[], T domain_y[]) {

    unsigned long steps = 50;

    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 rng(fullSeed);
    std::uniform_real_distribution<T> uniformDist(0.0f, 1.0f);
    std::normal_distribution<T> normDist(0.0f, 1.0f);

    T acceptance = 0.0;

    for(auto i = 0; i < steps; i ++) {

        T step_x = normDist(rng) * size_x;
        T step_y = normDist(rng) * size_y;
        T u_sample = uniformDist(rng);
        T theta_r1_p = theta_r1 + step_x;
        T theta_r2_p = theta_r2 + step_y;
        acceptance = 0.0;
        if(theta_r1_p < 0.0) {
            theta_r1_p = std::abs(theta_r1_p);
            theta_r2_p += M_PI;
        }
        if(theta_r2_p < 0.0) {
            theta_r2_p += 2.0 * M_PI;
        }
        if(theta_r2_p > 2.0 * M_PI) {
            theta_r2_p -= 2.0 * M_PI;
        }

        if(theta_r1_p < domain_x[1]) {
            acceptance = std::min(1.0, (this->*pdf)(theta_r1_p, theta_r2_p, parameters) / (this->*pdf)(theta_r1, theta_r2, parameters));
        }

        if(acceptance > u_sample) {
            theta_r1 = theta_r1_p;
            theta_r2 = theta_r2_p;
        }
    }

    return std::make_tuple(theta_r1, theta_r2);
}

template <typename T>
T PG_Kernel<T>::mixture_pdf(T gamma_local, T gamma_dot_local, std::vector<T> parameters) {

    bool first = (bool)parameters[0];
    T theta_i = parameters[1];
    T distance = parameters[2];

    T mu_local = herm_tools.evaluate(surface.get_mu_coefficients(), gamma_local);
    T mu_dot_local = herm_tools.evaluate_deriv(surface.get_mu_coefficients(), gamma_local) * gamma_dot_local;
    T sigma_local = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma_local);
    T sigma_gamma_local = herm_tools.evaluate_deriv(surface.get_sigma_coefficients(), gamma_local);

    T w = std::sqrt(2.0 * (sigma_local * sigma_local / ac_length / ac_length) + sigma_gamma_local * sigma_gamma_local * gamma_dot_local * gamma_dot_local);
    T eta = 1.0 / std::tan(theta_i);
    T slope_term = 0.5 * (1.0 + std::erf((eta + mu_dot_local) / (w * std::sqrt(2.0)))) * (std::cos(theta_i) + mu_dot_local * std::sin(theta_i)) + 
                    w * std::sin(theta_i) / std::sqrt(2.0 * M_PI) * std::exp(- (eta + mu_dot_local) * (eta + mu_dot_local) / 2.0 / w / w);
    T p_gamma_local = 1.0 / std::sqrt(2.0 * M_PI) * std::exp( - gamma_local * gamma_local / 2.0);
    T p_gamma_dot_local = ac_length / std::sqrt(4.0 * M_PI) * std::exp(- gamma_dot_local * gamma_dot_local * ac_length * ac_length / 4.0);

    if(first) {
        return p_gamma_local * p_gamma_dot_local * slope_term * shadow(theta_i, mu_local, -1);
    }
    else {
        return p_gamma_local * p_gamma_dot_local * slope_term * occlusion(theta_i, mu_local, distance, -1);
    }

}

template <typename T>
T PG_Kernel<T>::kirchhoff_pdf(T theta_r1, T theta_r2, std::vector<T> parameters) {

    T theta_i = parameters[0];
    T gamma = parameters[1];
    T gamma_old = parameters[2];
    T gamma_dot_x = parameters[3];
    T gamma_dot_x_old = parameters[4];
    T gamma_dot_y = parameters[5];
    T gamma_dot_y_old = parameters[6];
    T slope_x = parameters[7];
    T slope_y = parameters[8];
    T tau = parameters[9];

    T C = std::exp(- tau * tau / ac_length / ac_length);

    T mu_gamma = herm_tools.evaluate_deriv(surface.get_mu_coefficients(), gamma);
    T mu_gamma_old = herm_tools.evaluate_deriv(surface.get_mu_coefficients(), gamma_old);
    T sigma = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma);
    T sigma_old = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma_old);
    T sigma_gamma = herm_tools.evaluate_deriv(surface.get_sigma_coefficients(), gamma);
    T sigma_gamma_old = herm_tools.evaluate_deriv(surface.get_sigma_coefficients(), gamma_old);

    T mu_kirchhoff_x = mu_gamma * gamma_dot_x;
    T mu_kirchhoff_y = mu_gamma * gamma_dot_y;
    T mu_kirchhoff_x_old = mu_gamma_old * gamma_dot_x_old;
    T mu_kirchhoff_y_old = mu_gamma_old * gamma_dot_y_old;

    T vx = std::sin(theta_i) - std::sin(theta_r1) * std::cos(theta_r2);
    T vy = - std::sin(theta_r1) * std::sin(theta_r2);
    T vz = - std::cos(theta_i) - std::cos(theta_r1);

    T vx_old = - slope_x * vz;
    T vy_old = - slope_y * vz;
    T vz_old = std::sqrt(1.0 + slope_x * slope_x + slope_y * slope_y);

    T F = (1.0 + std::cos(theta_i) * std::cos(theta_r1) - std::sin(theta_i) * std::sin(theta_r1) * std::cos(theta_r2)) / 
            (std::cos(theta_i) * std::cos(theta_i) + std::cos(theta_r1) * std::cos(theta_i));

    T sigma_kirchhoff_x = std::sqrt((sigma * sigma / ac_length / ac_length) + 0.5 * sigma_gamma * sigma_gamma * gamma_dot_x * gamma_dot_x);
    T sigma_kirchhoff_y = std::sqrt((sigma * sigma / ac_length / ac_length) + 0.5 * sigma_gamma * sigma_gamma * gamma_dot_y * gamma_dot_y);
    T sigma_kirchhoff_x_old = std::sqrt((sigma_old * sigma_old / ac_length / ac_length) + 0.5 * sigma_gamma_old * sigma_gamma_old * gamma_dot_x_old * gamma_dot_x_old);
    T sigma_kirchhoff_y_old = std::sqrt((sigma_old * sigma_old / ac_length / ac_length) + 0.5 * sigma_gamma_old * sigma_gamma_old * gamma_dot_y_old * gamma_dot_y_old);


    T vx_k = (vx + mu_kirchhoff_x * vz);
    T vy_k = (vy + mu_kirchhoff_y * vz);
    T vx_k_old = (vx_old + mu_kirchhoff_x_old * vz);
    T vy_k_old = (vy_old + mu_kirchhoff_y_old * vz);

    // T p_kirchhoff_x = std::sqrt(2.0 * M_PI) * F / vz / std::sqrt(1.0 - C * C) * 
    //                     std::exp(- (vx_k * vx_k - 2.0 * C * vx_k * (slope_x - mu_kirchhoff_x) * vz + vz * vz * (slope_x - mu_kirchhoff_x) * (slope_x - mu_kirchhoff_x))
    //                             / (4.0 * vz * vz * sigma_kirchhoff_x * sigma_kirchhoff_x * (1.0 - C * C)));
    // T p_kirchhoff_y = std::sqrt(2.0 * M_PI) * F / vz / std::sqrt(1.0 - C * C) * 
    //                     std::exp(- (vy_k * vy_k - 2.0 * C * vy_k * (slope_y - mu_kirchhoff_y) * vz + vz * vz * (slope_y - mu_kirchhoff_y) * (slope_y - mu_kirchhoff_y))
    //                             / (4.0 * vz * vz * sigma_kirchhoff_y * sigma_kirchhoff_y * (1.0 - C * C)));
    
    T p_kirchhoff_x = std::sqrt(2.0 * M_PI) * F / vz / std::sqrt(1.0 - C * C) * std::exp(- 1.0 / 4.0 / (1.0 - C * C) * (vx_k * vx_k / sigma_kirchhoff_x / sigma_kirchhoff_x / vz / vz
                                + vx_k_old * vx_k_old / sigma_kirchhoff_x_old / sigma_kirchhoff_x_old / vz_old / vz_old - 
                                2.0 * C * vx_k * vx_k_old / sigma_kirchhoff_x / sigma_kirchhoff_x_old / vz / vz_old));
    T p_kirchhoff_y = std::sqrt(2.0 * M_PI) * F / vz / std::sqrt(1.0 - C * C) * std::exp(- 1.0 / 4.0 / (1.0 - C * C) * (vy_k * vy_k / sigma_kirchhoff_y / sigma_kirchhoff_y / vz / vz
                                + vy_k_old * vy_k_old / sigma_kirchhoff_y_old / sigma_kirchhoff_y_old / vz_old / vz_old - 
                                2.0 * C * vy_k * vy_k_old / sigma_kirchhoff_y / sigma_kirchhoff_y_old / vz / vz_old));

    return p_kirchhoff_x * p_kirchhoff_y * std::sin(theta_r1);
}

template <typename T>
T PG_Kernel<T>::shadow(T theta_r1, T height, int direction) {

    // Check if theta_r1 is greater than 90 degrees
    if(theta_r1 > M_PI / 2.0) {
        return 0.0;
    }
    T EPSILON = 1e-15;
    theta_r1 = std::max(EPSILON, theta_r1);
    T eta = 1.0 / std::tan(theta_r1);
    Matrix<T> w_2d = (sigma_2d * sigma_2d / ac_length / ac_length * 2.0 + sigma_gamma_2d * sigma_gamma_2d * gamma_dot_2d * gamma_dot_2d) ^ 0.5;
    Matrix<T> Delta = (w_2d / std::sqrt(2.0 * M_PI) * exp(((mu_gamma_2d * gamma_dot_2d * (-1.0) * direction + eta) ^ 2.0) / ((w_2d * w_2d) * 2.0) * (-1.0)) - 
                        (mu_gamma_2d * gamma_dot_2d * (-1.0) * direction + eta) * 0.5 * (erfc((mu_gamma_2d * gamma_dot_2d * (-1.0) * direction + eta) / (w_2d * std::sqrt(2.0)))));
    T exponent = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * Delta) / eta;
    T base = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * 0.5 * (erfc((mu_2d * (-1.0) + height) / (sigma_2d * std::sqrt(2.0))) * (-1.0) + 2.0));
    base = std::min(1.0, base);
    return std::pow(base, exponent);
}

template <typename T>
T PG_Kernel<T>::occlusion(T theta_r1, T height, T distance, int direction) {

    theta_r1 = std::max(1e-12, theta_r1);
    T eta = 1.0 / std::tan(theta_r1);
    T height_new = height + distance * eta;
    Matrix<T> w_2d = (sigma_2d * sigma_2d / ac_length / ac_length * 2.0 + sigma_gamma_2d * sigma_gamma_2d * gamma_dot_2d * gamma_dot_2d) ^ 0.5;
    Matrix<T> Delta = (w_2d / std::sqrt(2.0 * M_PI) * exp(((mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) ^ 2.0) / ((w_2d * w_2d) * 2.0) * (-1.0)) - 
                        (mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) * 0.5 * (erfc((mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) / (w_2d * std::sqrt(2.0)))));
    T exponent = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * Delta) / eta;
    T base_old = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * 0.5 * (erfc((mu_2d * (-1.0) + height) / (sigma_2d * std::sqrt(2.0))) * (-1.0) + 2.0));
    T base_new = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * 0.5 * (erfc((mu_2d * (-1.0) + height_new) / (sigma_2d * std::sqrt(2.0))) * (-1.0) + 2.0));

    return std::pow(base_old / base_new, exponent);
}

template <typename T>
T PG_Kernel<T>::distance_pdf(T distance, std::vector<T> parameters) {

    T theta_r1 = std::max(parameters[0], 1e-12);
    T height_0 = parameters[1];
    T gamma_0 = parameters[2];
    T gamma_dot_0 = parameters[3];
    
    T eta = 1.0 / std::tan(theta_r1);
    T height_new = height_0 + eta * distance;
    T mu_0 = herm_tools.evaluate(surface.get_mu_coefficients(), gamma_0);
    T sigma_0 = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma_0);
    Matrix<T> w_2d = (sigma_2d * sigma_2d / ac_length / ac_length * 2.0 + sigma_gamma_2d * sigma_gamma_2d * gamma_dot_2d * gamma_dot_2d) ^ 0.5;
    Matrix<T> Delta = (w_2d / std::sqrt(2.0 * M_PI) * exp(((mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) ^ 2.0) / ((w_2d * w_2d) * 2.0) * (-1.0)) - 
                        (mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) * 0.5 * (erfc((mu_gamma_2d * gamma_dot_2d * (-1.0) + eta) / (w_2d * std::sqrt(2.0)))));
    T exponent = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * Delta) / eta;
    T base_old = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * 0.5 * (erfc((mu_2d * (-1.0) + height_0) / (sigma_2d * std::sqrt(2.0))) * (-1.0) + 2.0));
    T base_new = sum(weights_2d * p_gamma_2d * p_gamma_dot_2d * 0.5 * (erfc((mu_2d * (-1.0) + height_new) / (sigma_2d * std::sqrt(2.0))) * (-1.0) + 2.0));
    T p_h = 1.0 / std::sqrt(2.0 * M_PI) / sigma_0 * std::exp(- (height_new - mu_0) * (height_new - mu_0) / 2.0 / sigma_0 / sigma_0);
    T S = std::pow(base_old, exponent) / std::pow(base_new, exponent + 1.0);

    return S * p_h;
}

template <typename T>
T PG_Kernel<T>::height_0_pdf(T height_0, std::vector<T> parameters) {

    T theta_i = parameters[0];
    T gamma_0 = parameters[1];
    T gamma_dot_0 = parameters[2];

    T mu_0 = herm_tools.evaluate(surface.get_mu_coefficients(), gamma_0);
    T sigma_0 = herm_tools.evaluate(surface.get_sigma_coefficients(), gamma_0);

    T p_h = 1.0 / std::sqrt(2.0 * M_PI) / sigma_0 * std::exp(- std::pow(height_0 - mu_0, 2.0) / 2.0 / sigma_0 / sigma_0);

    return p_h * shadow(theta_i, height_0, -1);
}

template <typename T>
Matrix<T> PG_Kernel<T>::global_to_local(Matrix<T> global_vector, T theta_i, T theta_r1, T theta_r2) {

    Matrix<T> unit_incident = global_vector / global_vector.norm();
    Matrix<T> z_local(3, 1);
    z_local(0, 0) = std::sin(theta_r1) * std::cos(theta_r2) - std::sin(theta_i);
    z_local(1, 0) = std::sin(theta_r1) * std::sin(theta_r2);
    z_local(2, 0) = std::cos(theta_r1) + std::cos(theta_i);
    z_local /= z_local.norm();
    Matrix<T> x_local(3, 1);
    x_local = unit_incident - z_local * unit_incident.tr().dot(z_local)(0, 0);
    x_local /= x_local.norm();
    Matrix<T> y_local(3, 1);
    y_local = x_local.cross(z_local);
    y_local /= y_local.norm();

    T res_x = global_vector.tr().dot(x_local)(0, 0);
    T res_y = global_vector.tr().dot(y_local)(0, 0);
    T res_z = global_vector.tr().dot(z_local)(0, 0);

    Matrix<T> res(3, 1, {res_x, res_y, res_z});
    
    return res;
}

template <typename T>
Matrix<T> PG_Kernel<T>::local_to_global(Matrix<T> local_vector, T theta_i, T theta_r1, T theta_r2) {

    Matrix<T> unit_incident(3, 1, {std::sin(theta_i), 0.0, -std::cos(theta_i)});
    Matrix<T> z_local(3, 1);
    z_local(0, 0) = std::sin(theta_r1) * std::cos(theta_r2) - std::sin(theta_i);
    z_local(1, 0) = std::sin(theta_r1) * std::sin(theta_r2);
    z_local(2, 0) = std::cos(theta_r1) + std::cos(theta_i);
    z_local /= z_local.norm();
    Matrix<T> x_local(3, 1);
    x_local = unit_incident - z_local * unit_incident.tr().dot(z_local)(0, 0);
    x_local /= x_local.norm();
    Matrix<T> y_local(3, 1);
    y_local = x_local.cross(z_local);
    y_local /= y_local.norm();
    
    return x_local * local_vector(0, 0) + y_local * local_vector(1, 0) + z_local * local_vector(2, 0);
}

template <typename T>
Matrix<T> PG_Kernel<T>::angles_to_slopes(T theta_i, T theta_r1, T theta_r2) {

    Matrix<T> unit_incident(3, 1, {std::sin(theta_i), 0.0, -std::cos(theta_i)});
    Matrix<T> unit_reflected(3, 1, {std::sin(theta_r1) * std::cos(theta_r2), std::sin(theta_r1) * std::sin(theta_r2), std::cos(theta_r1)});
    Matrix<T> normal(3, 1);
    normal = unit_reflected - unit_incident;
    normal /= std::abs(normal(2, 0));
    
    return normal;
}

template <typename T>
Matrix<T> PG_Kernel<T>::sample_particle(Matrix<T> incident_velocity, T* gsi_parameters) {

    Matrix<T> final_reflected_velocity(3, 1, {0.0, 0.0, 0.0});

    T incident_norm = incident_velocity.norm();
    T theta_i = std::acos(- incident_velocity(2, 0) / incident_norm);
    T theta_r1 = 0.0, theta_r2 = 0.0, theta_r2_total = 0.0, distance = 1e8;
    T gamma_0 = 0.0, gamma_dot_0_x = 0.0, gamma_dot_0_y = 0.0, slope_x = 0.0, slope_y = 0.0, slope_x_new =0.0, slope_y_new = 0.0;

    long num_col = 0;
    bool collided = true;

    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 rng(fullSeed);
    std::normal_distribution<T> normDist(0.0f, 1.0f);

    T gamma, gamma_dot_x, gamma_dot_y, gamma_old = 0.0, gamma_dot_x_old = 1.0, gamma_dot_y_old = 1.0;
    T domain[] = {- margin, margin};
    std::tie(gamma, gamma_dot_x) = random_sampler_2D(&PG_Kernel::mixture_pdf, {1.0, theta_i, 0.0}, normDist(rng), normDist(rng), 0.5, 0.5, domain, domain);
    gamma_dot_y = normDist(rng) * std::sqrt(2.0) / ac_length;

    T height_domain[] = {-margin, margin};
    T height_0 = random_sampler_1D(&PG_Kernel::height_0_pdf, {theta_i, gamma, gamma_dot_x}, normDist(rng), 0.5, height_domain);

    while(collided) {

        T angle_domain_x[] = {0.0, M_PI - theta_i};
        T angle_domain_y[] = {-M_PI, M_PI};
        std::tie(theta_r1, theta_r2) = random_sampler_angles(&PG_Kernel::kirchhoff_pdf, {theta_i, gamma, gamma_old, gamma_dot_x, gamma_dot_x_old, gamma_dot_y, gamma_dot_y_old, slope_x, slope_y, distance}, 1e-5, 1e-5, 
                            (M_PI - theta_i) * 0.25, (M_PI - theta_i) * 0.25, angle_domain_x, angle_domain_y);
        Matrix<T> slopes_matrix = angles_to_slopes(theta_i, theta_r1, theta_r2);
        slope_x_new = -slopes_matrix(0, 0);
        slope_y_new = -slopes_matrix(1, 0);

        Matrix<T> incident_velocity_local = global_to_local(incident_velocity, theta_i, theta_r1, theta_r2);
        T vel_in_array[] = {incident_velocity_local(0, 0), incident_velocity_local(1, 0), incident_velocity_local(2, 0)};
        T vel_refl_array[] = {0.0, 0.0, 0.0};

        local_kernel(vel_in_array, vel_refl_array, gsi_parameters, gas, surface);

        Matrix<T> reflected_velocity_local(3, 1, {vel_refl_array[0], vel_refl_array[1], vel_refl_array[2]});
        Matrix<T> reflected_velocity = local_to_global(reflected_velocity_local, theta_i, theta_r1, theta_r2);

        theta_r1 = std::acos(reflected_velocity(2, 0) / reflected_velocity.norm());
        theta_r2 = std::atan2(reflected_velocity(1, 0), reflected_velocity(0, 0));

        std::bernoulli_distribution shadow_pdf(1.0 - shadow(theta_r1, height_0, 1));
        collided = shadow_pdf(rng);

        T distance_domain[] = {0.0, margin * 5.0};
        distance = random_sampler_1D(&PG_Kernel::distance_pdf, {theta_r1, height_0, gamma, gamma_dot_x}, 0.0, 0.5, distance_domain);
        height_0 += distance / std::tan(theta_r1);
        theta_i = M_PI - theta_r1;
        theta_r2_total += theta_r2;

        incident_norm = reflected_velocity.norm();
        incident_velocity(0, 0) = incident_norm * std::sin(theta_i);
        incident_velocity(1, 0) = 0.0;
        incident_velocity(2, 0) = - incident_norm * std::cos(theta_i);

        Matrix<T> global_reflected_velocity(3, 1, {std::sin(theta_r1) * std::cos(theta_r2_total), 
                                                    std::sin(theta_r1) * std::sin(theta_r2_total), 
                                                    std::cos(theta_r1)});
        global_reflected_velocity *= incident_norm;
        Matrix<T> position_change(3, 1, {std::sin(theta_r1) * std::cos(theta_r2_total), 
                std::sin(theta_r1) * std::sin(theta_r2_total), std::cos(theta_r1)});
        position_change *= distance / std::sin(theta_r1);

        slope_x = - std::cos(theta_r2) * slope_x_new - std::sin(theta_r2) * slope_y_new;
        slope_y = std::sin(theta_r2) * slope_x_new - std::cos(theta_r2) * slope_y_new;
        gamma_old = gamma;
        gamma_dot_x_old = gamma_dot_x;
        gamma_dot_y_old = gamma_dot_y;

        std::tie(gamma, gamma_dot_x) = random_sampler_2D(&PG_Kernel::mixture_pdf, {0.0, theta_i, distance}, normDist(rng), normDist(rng), 0.5, 0.5, domain, domain);
        gamma_dot_y = normDist(rng) * std::sqrt(2.0) / ac_length;
        num_col ++;
        final_reflected_velocity = global_reflected_velocity;
    }
    return final_reflected_velocity;
}
