/**
 * @brief Template class representing a Gas with various properties.
 * 
 * @tparam T Type of the gas properties (e.g., float, double).
 */
template <typename T>
class Gas {

    private:
        T density;       ///< The density of the gas.
        T temperature;   ///< The temperature of the gas.
        T molar_mass;    ///< The molar mass of the gas.
        T speed;         ///< The average speed of the gas molecules.

    public:
        /**
         * @brief Default constructor initializing all properties to 0.0.
         */
        Gas() : density(0.0), temperature(0.0), molar_mass(0.0), speed(0.0) {}

        /**
         * @brief Parameterized constructor initializing with specific values.
         * 
         * @param density Initial density of the gas.
         * @param temperature Initial temperature of the gas.
         * @param molar_mass Initial molar mass of the gas.
         * @param speed Initial average speed of the gas molecules.
         */
        Gas(T density, T temperature, T molar_mass, T speed) 
            : density(density), temperature(temperature), molar_mass(molar_mass), speed(speed) {}

        /**
         * @brief Copy constructor.
         * 
         * @param other Another Gas object to copy from.
         */
        Gas(const Gas<T>& other) 
            : density(other.density), temperature(other.temperature), molar_mass(other.molar_mass), speed(other.speed) {}

        /**
         * @brief Move constructor.
         * 
         * @param other Another Gas object to move from.
         * After moving, the original object will be reset to 0.0.
         */
        Gas(Gas<T>&& other) 
            : density(other.density), temperature(other.temperature), molar_mass(other.molar_mass), speed(other.speed) 
        {
            // Reset the source object to a default state after moving.
            other.density = 0.0;
            other.temperature = 0.0;
            other.molar_mass = 0.0;
            other.speed = 0.0;
        }

        /**
         * @brief Destructor resetting all values to 0.0.
         */
        ~Gas() {
            density = 0.0;
            temperature = 0.0;
            molar_mass = 0.0;
            speed = 0.0;
        }

        /**
         * @brief Copy assignment operator.
         * 
         * @param other Another Gas object to copy from.
         * @return Gas& Reference to the updated object.
         */
        Gas& operator=(const Gas& other) {
            density = other.density;
            temperature = other.temperature;
            molar_mass = other.molar_mass;
            speed = other.speed;
            return *this;
        }

        /**
         * @brief Move assignment operator.
         * 
         * @param other Another Gas object to move from.
         * After moving, the original object will be reset to 0.0.
         * @return Gas& Reference to the updated object.
         */
        Gas& operator=(Gas&& other) {
            density = other.density;
            temperature = other.temperature;
            molar_mass = other.molar_mass;
            speed = other.speed;

            // Reset the source object to a default state after moving.
            other.density = 0.0;
            other.temperature = 0.0;
            other.molar_mass = 0.0;
            other.speed = 0.0;
            return *this;
        }

        /**
         * @brief Calculates the thermal velocity of the gas molecules.
         * 
         * @return Matrix<T> A 3x1 matrix representing the velocity components (v_x, v_y, v_z).
         */
        Matrix<T> get_thermal_velocity();

        /**
         * @brief Getter for the temperature.
         * 
         * @return T Current temperature of the gas.
         */
        T get_temperature() const { return temperature; }

        /**
         * @brief Getter for the density.
         * 
         * @return T Current density of the gas.
         */
        T get_density() const { return density; }

        /**
         * @brief Getter for the molar mass.
         * 
         * @return T Current molar mass of the gas.
         */
        T get_molar_mass() const { return molar_mass; }

        /**
         * @brief Getter for the speed.
         * 
         * @return T Current speed of the gas molecules.
         */
        T get_speed() const { return speed; }
};

/**
 * @brief Calculates the thermal velocity of the gas molecules based on random distribution.
 * 
 * This function generates the velocity components (v_x, v_y, v_z) using a normal distribution 
 * and the Maxwell-Boltzmann distribution equation.
 * 
 * @tparam T Type of the gas properties (e.g., float, double).
 * @return Matrix<T> A 3x1 matrix representing the velocity components (v_x, v_y, v_z).
 */
template <typename T>
Matrix<T> Gas<T>::get_thermal_velocity() {

    // Seed for random number generation
    std::random_device rd;
    std::seed_seq fullSeed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937 rng(fullSeed);  // Mersenne Twister RNG
    std::normal_distribution<T> normDist(0.0f, 1.0f);  // Normal distribution with mean 0 and std 1

    // Calculate the thermal velocities in x, y, and z directions
    T v_x = sqrt(R_gas * temperature / molar_mass) * normDist(rng);
    T v_y = sqrt(R_gas * temperature / molar_mass) * normDist(rng);
    T v_z = sqrt(R_gas * temperature / molar_mass) * normDist(rng);

    // Store the velocities in a 3x1 matrix
    Matrix<T> vel(3, 1, {v_x, v_y, v_z});

    return vel;  // Return the velocity matrix
}
