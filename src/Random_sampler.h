
/**
 * @brief Computes the sign of a value.
 * 
 * This function returns -1 if the input value is negative, 1 if the input value is positive, 
 * and 0 if the input value is zero.
 * 
 * @tparam T Type of the input value (must be comparable).
 * @param val The value whose sign is to be determined.
 * @return int -1 if val < 0, 1 if val > 0, and 0 if val == 0.
 */
template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief Samples the velocity after reflection using the Cercignani-Lampis-Lord (CLL) gas-surface interaction model.
 *
 * This function computes the post-collision reflected velocities of gas molecules based on the input velocity and 
 * the gas-surface interaction parameters. It uses random number generation to model stochastic scattering events.
 *
 * @tparam T Type of the velocity and parameter values (e.g., float, double).
 * @param vel_in Input velocity of the gas particle (3D vector: x, y, z).
 * @param vel_refl Output reflected velocity of the gas particle (3D vector: x, y, z).
 * @param GSI_parameters Gas-Surface Interaction parameters. GSI_parameters[0] = normal accommodation coefficient (alpha_n),
 * GSI_parameters[1] = tangential momentum accommodation coefficient (sigma_t).
 * @param gas Gas object containing properties like temperature and molar mass.
 * @param surface Surface object representing the boundary condition (e.g., temperature of the wall).
 */
template <typename T>
void sample_CLL(T* vel_in, T* vel_refl, T* GSI_parameters, const Gas<T>& gas, Surface<T>& surface) {

    // Extract Gas-Surface Interaction (GSI) parameters
    T alpha_n = GSI_parameters[0];  ///< Normal accommodation coefficient
    T sigma_t = GSI_parameters[1];  ///< Tangential momentum accommodation coefficient

    // Compute the tangential momentum accommodation parameter
    T alpha_t = sigma_t * (2 - sigma_t);

    // Get the wall temperature from the surface object
    T T_W = surface.get_temperature();

    // Set up random number generation
    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};  ///< Seed sequence for RNG
    std::mt19937 rng(fullSeed);  ///< Mersenne Twister RNG engine
    std::uniform_real_distribution<T> uniformDist(0.0f, 1.0f);  ///< Uniform distribution between 0 and 1

    // Calculate the thermal speed of the gas molecules at the wall temperature
    T vel_wall = sqrt(2.0 * R_gas / gas.get_molar_mass() * T_W);

    // Generate random values for calculating the reflected velocities
    T x1 = uniformDist(rng), x2 = uniformDist(rng), x3 = uniformDist(rng), x4 = uniformDist(rng), x5 = uniformDist(rng), x6 = uniformDist(rng);

    // Calculate the reflection parameters based on random numbers
    T r1 = sqrt(-alpha_n * log(x1));  ///< Normal component for reflection
    T r3 = sqrt(-alpha_t * log(x3));  ///< First tangential component
    T r5 = sqrt(-alpha_t * log(x5));  ///< Second tangential component

    // Calculate the angles for the velocity reflection
    T theta2 = 2.0 * M_PI * x2;
    T theta4 = 2.0 * M_PI * x4;
    T theta6 = 2.0 * M_PI * x6;

    // Compute modified velocities
    T vel_norm_m = fabs(vel_in[2] / vel_wall) * sqrt(1 - alpha_n);   ///< Normalized normal velocity component
    T vel_tan1_m = vel_in[0] / vel_wall * sqrt(1 - alpha_t);         ///< First tangential velocity component

    // Calculate the reflected velocity components
    vel_refl[2] = vel_wall * sqrt(r1 * r1 + vel_norm_m * vel_norm_m + 2 * r1 * vel_norm_m * cos(theta2));  ///< Reflected normal velocity
    vel_refl[0] = vel_wall * (vel_tan1_m + r3 * cos(theta4));  ///< Reflected first tangential velocity
    vel_refl[1] = vel_wall * (r5 * cos(theta6));  ///< Reflected second tangential velocity
}


/**
 * @brief Samples the velocity after reflection using the Diffuse Reemission with Incomplete Accommodation (DRIA) model.
 *
 * This function computes the post-collision reflected velocities of gas molecules using the DRIA model. It considers
 * both the temperature of the gas and the temperature of the surface to determine the reflected velocity distribution.
 *
 * @tparam T Type of the velocity and parameter values (e.g., float, double).
 * @param vel_in Input velocity of the gas particle (3D vector: x, y, z).
 * @param vel_refl Output reflected velocity of the gas particle (3D vector: x, y, z).
 * @param GSI_parameters Gas-Surface Interaction parameters. GSI_parameters[0] = accommodation coefficient (alpha).
 * @param gas Gas object containing properties like temperature, molar mass, and speed.
 * @param surface Surface object representing the boundary condition (e.g., temperature of the wall).
 */
template <typename T>
void sample_DRIA(T* vel_in, T* vel_refl, T* GSI_parameters, const Gas<T>& gas, Surface<T>& surface) {

    // Calculate the magnitude of the input velocity
    T vel_in_norm = sqrt(vel_in[0] * vel_in[0] + vel_in[1] * vel_in[1] + vel_in[2] * vel_in[2]);

    // Extract the accommodation coefficient (alpha) from GSI parameters
    T alpha = GSI_parameters[0];

    // Get the temperatures from the surface and gas objects
    T T_W = surface.get_temperature();  ///< Wall temperature (surface)
    T T_G = gas.get_temperature();      ///< Gas temperature

    // Get the gas speed from the gas object
    T gas_speed = gas.get_speed();

    // Compute the characteristic temperature of the gas using its speed and molar mass
    T T_K = gas_speed * gas_speed * gas.get_molar_mass() / 3.0 / R_gas;

    // Set up random number generation
    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};  ///< Seed sequence for RNG
    std::mt19937 rng(fullSeed);  ///< Mersenne Twister RNG engine
    std::uniform_real_distribution<T> uniformDist(0.0f, 1.0f);  ///< Uniform distribution between 0 and 1
    std::normal_distribution<T> normDist(0.0f, 1.0f);  ///< Normal distribution with mean 0 and std dev 1

    // Compute the effective reflection temperature (T_R) as a weighted average of T_K and T_W
    T T_R = (1.0 - alpha) * T_K + alpha * T_W;

    // Compute the thermal velocity of gas molecules at temperature T_R
    T vel_wall = sqrt(2.0 * R_gas / gas.get_molar_mass() * T_R);

    // Generate random values for calculating the reflected velocities
    T x1 = uniformDist(rng), x2 = uniformDist(rng), x3 = uniformDist(rng), x4 = uniformDist(rng), x5 = uniformDist(rng), x6 = uniformDist(rng);

    // Calculate the reflection parameters based on random numbers
    T r1 = sqrt(-log(x1));  ///< Normal component for reflection
    T r3 = sqrt(-log(x3));  ///< First tangential component
    T r5 = sqrt(-log(x5));  ///< Second tangential component

    // Calculate the angles for the velocity reflection
    T theta2 = 2.0 * M_PI * x2;
    T theta4 = 2.0 * M_PI * x4;
    T theta6 = 2.0 * M_PI * x6;

    // Initialize normal and tangential velocity components
    T vel_norm_m = 0.0;
    T vel_tan_m = 0.0;

    // Calculate the reflected normal velocity
    vel_refl[2] = vel_wall * r1;  ///< Reflected normal velocity

    // Calculate the reflected tangential velocities based on random normal distribution
    vel_refl[0] = sqrt(R_gas * T_R / gas.get_molar_mass()) * normDist(rng);  ///< Reflected first tangential velocity
    vel_refl[1] = sqrt(R_gas * T_R / gas.get_molar_mass()) * normDist(rng);  ///< Reflected second tangential velocity
}


/**
 * @brief Samples the velocity after rough surface reflection using the Kirchhoff - Cercignani-Lampis-Lord (CLL) model.
 *
 * This function simulates the reflection of gas molecules from a rough surface. The incident velocity is altered
 * based on surface roughness, with the CLL model applied for particle reflection. The method uses random number generation
 * to model the stochastic nature of the surface roughness and scattering events.
 *
 * @tparam T Type of the velocity and parameter values (e.g., float, double).
 * @param vel_in Input velocity of the gas particle (3D vector: x, y, z).
 * @param vel_refl Output reflected velocity of the gas particle (3D vector: x, y, z).
 * @param GSI_parameters Gas-Surface Interaction parameters used in the CLL model.
 * @param gas Gas object containing properties like temperature and molar mass.
 * @param surface Surface object representing the boundary condition (e.g., temperature of the wall, roughness).
 */
template <typename T>
void sample_rough_CLL(T* vel_in, T* vel_refl, T* GSI_parameters, const Gas<T>& gas, Surface<T>& surface) {

    // Small constant to prevent division by zero or other numerical issues
    T epsilon = 1e-12;

    // Set up random number generation
    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};  ///< Seed sequence for RNG
    std::mt19937 rng(fullSeed);  ///< Mersenne Twister RNG engine

    // Initialize collision status
    bool collided = true;  ///< Variable indicating whether the particle collided with the surface

    // Calculate the initial incident angle (theta_i) of the gas particle
    T theta_i_initial = acos(fabs(vel_in[2]) / sqrt(vel_in[0] * vel_in[0] + vel_in[1] * vel_in[1] + vel_in[2] * vel_in[2]));

    // Variables for the reflected angles (set to zero initially)
    T theta_r1 = 0.0;
    T theta_r2 = 0.0;

    // Track the current incident angle (theta_i), initialized to the initial incident angle
    T theta_i = theta_i_initial;

    // Initialize the reflection counter
    unsigned int reflections = 0;

    // Matrix to hold the input velocity, ensuring the normal component (z) is negative (to simulate reflection)
    Matrix<T> vel_in_m(3, 1, {vel_in[0], vel_in[1], -fabs(vel_in[2])});

    // Matrix to store the output reflected velocity (initialized to zero)
    Matrix<T> vel_out_m(3, 1, {0.0, 0.0, 0.0});

    // Create a PG_Kernel object to handle rough surface interactions
    PG_Kernel<T> rough_kernel(surface, gas, sample_CLL);  ///< PG_Kernel uses the CLL model to sample particle reflection

    // Sample the reflected velocity using the rough surface kernel
    vel_out_m = rough_kernel.sample_particle(vel_in_m, GSI_parameters);

    // Assign the reflected velocity components to the output array
    vel_refl[0] = vel_out_m(0);  ///< Reflected velocity in the x-direction
    vel_refl[1] = vel_out_m(1);  ///< Reflected velocity in the y-direction
    vel_refl[2] = vel_out_m(2);  ///< Reflected velocity in the z-direction
}
