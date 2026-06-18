/**
 * @class Surface
 * @brief A class representing a surface with geometry and statistical properties.
 * 
 * The Surface class stores the geometric data and statistical properties (e.g., autocorrelation length,
 * mu and sigma coefficients) required for surface modeling and interaction in simulations.
 * 
 * @tparam T Type parameter for numerical data (typically float or double).
 */
template <typename T>
class Surface {

    private:

        /** 
         * @brief Mu coefficients for the surface properties.
         */
        Matrix<T> mu_coefficients;

        /** 
         * @brief Sigma coefficients for the surface properties.
         */
        Matrix<T> sigma_coefficients;

        /** 
         * @brief Autocorrelation length of the surface.
         */
        T ac_length, temperature;


    public:

        /** 
         * @brief Default constructor for Surface class.
         * 
         * Initializes the autocorrelation length to 0.
         */
        Surface();

        /**
         * @brief Parameterized constructor for Surface class.
         * 
         * @param properties_file The file containing surface properties (e.g., mu, sigma coefficients).
         */
        Surface(std::string properties_file);

        /**
         * @brief Copy constructor. Creates a copy of an existing Surface object.
         * 
         * This constructor creates a new Surface object with the same dimensions 
         * and copies all elements from the original surface to the new one.
         * 
         * @param other The Surface object to copy from.
         */
        Surface(const Surface<T>& other);

        /**
         * @brief Move constructor. Moves the data from a given Surface object to the new one.
         * 
         * This constructor "steals" the internal data from the original object and 
         * leaves it in an empty state. This is more efficient than copying when the 
         * original object is no longer needed.
         * 
         * @param other The Surface object whose data is to be moved.
         */
        Surface(Surface<T>&& other);

        /**
         * @brief Destructor for Surface class.
         * 
         * Resets the autocorrelation length to 0.
         */
        ~Surface();

        /**
         * @brief Move assignment operator. Moves the data from the given Surface object.
         * 
         * This operator "steals" the internal data from the 'other' Surface, leaving it in an empty state.
         * 
         * @param other The Surface object to move data from.
         * @return A reference to the current object, to enable operator chaining.
         */
        Surface& operator=(const Surface<T>& other);

         /**
         * @brief Constant assignment operator. Assigns all Surface values to one constant.
         * 
         * @param constant The constant to be assigned.
         * @return A reference to the current object, to enable operator chaining.
         */
        Surface& operator=(Surface<T>&& other);

        /**
         * @brief Imports surface properties from a file.
         * 
         * Reads the surface properties, such as mu and sigma coefficients, and autocorrelation length
         * from the provided properties file.
         * 
         * @param properties_file The file containing surface properties.
         */
        void import_properties(std::string properties_file);

        /**
         * @brief Sets the temperature of the surface.
         * 
         * @param temperature The temperature to be set.
         */
        void set_temperature(T temperature) { this->temperature = temperature; }

        /**
         * @brief Retrieves the mu coefficients for the surface.
         * 
         * @return The mu coefficients as a Matrix object.
         */
        Matrix<T>& get_mu_coefficients() { return mu_coefficients; }

        /**
         * @brief Retrieves the sigma coefficients for the surface.
         * 
         * @return The sigma coefficients as a Matrix object.
         */
        Matrix<T>& get_sigma_coefficients() { return sigma_coefficients; }

        /**
         * @brief Retrieves the autocorrelation length of the surface.
         * 
         * @return The autocorrelation length of the surface.
         */
        T autocorrelation_length() { return ac_length; }

        /**
         * @brief Retrieves the temperature of the surface.
         * 
         * @return The temperature of the surface.
         */
        T get_temperature() { return temperature; }

};

/**
 * @brief Default constructor for Surface.
 * 
 * Initializes the autocorrelation length to 0.
 */
template <typename T>
Surface<T>::Surface() : ac_length(0.0), temperature(300.0) {}

/**
 * @brief Parameterized constructor for Surface.
 * 
 * Initializes the geometry from a file and imports surface properties from another file.
 * 
 * @param properties_file The file containing surface properties (e.g., mu, sigma coefficients).
 */
template <typename T>
Surface<T>::Surface(std::string properties_file) {
    import_properties(properties_file);
    this->temperature = 300.0;
}

 /**
 * @brief Copy constructor. Creates a copy of an existing Surface object.
 * 
 * This constructor creates a new Surface object with the same dimensions 
 * and copies all elements from the original surface to the new one.
 * 
 * @param other The Surface object to copy from.
 */
template <typename T>
Surface<T>::Surface(const Surface<T>& other) : ac_length(other.ac_length), temperature(other.temperature), mu_coefficients(other.mu_coefficients),
sigma_coefficients(other.sigma_coefficients) {}

 /**
 * @brief Move constructor. Moves the data from a given Surface object to the new one.
 * 
 * This constructor "steals" the internal data from the original object and 
 * leaves it in an empty state. This is more efficient than copying when the 
 * original object is no longer needed.
 * 
 * @param other The Surface object whose data is to be moved.
 */
template <typename T>
Surface<T>::Surface(Surface<T>&& other) : ac_length(other.ac_length), temperature(other.temperature), mu_coefficients(other.mu_coefficients),
sigma_coefficients(other.sigma_coefficients) {
    other.ac_length = 0.0;
    other.temperature = 0.0;
}

/**
 * @brief Destructor for Surface.
 * 
 * Sets the autocorrelation length to 0 upon destruction.
 */
template <typename T>
Surface<T>::~Surface() {
    ac_length = 0.0;
    temperature = 0.0;
}

 /**
 * @brief Move assignment operator. Moves the data from the given Surface object.
 * 
 * This operator "steals" the internal data from the 'other' Surface, leaving it in an empty state.
 * 
 * @param other The Surface object to move data from.
 * @return A reference to the current object, to enable operator chaining.
 */
template <typename T>
Surface<T>& Surface<T>::operator=(const Surface<T>& other) {

    this->ac_length = other.ac_length;
    this->temperature = other.temperature;
    this->mu_coefficients = other.mu_coefficients;
    this->sigma_coefficients = other.sigma_coefficients;

    return *this;
}

 /**
 * @brief Constant assignment operator. Assigns all Surface values to one constant.
 * 
 * @return A reference to the current object, to enable operator chaining.
 */
template <typename T>
Surface<T>& Surface<T>::operator=(Surface<T>&& other) {

    this->ac_length = other.ac_length;
    this->temperature = other.temperature;
    this->mu_coefficients = other.mu_coefficients;
    this->sigma_coefficients = other.sigma_coefficients;
    other.ac_length = 0.0;
    other.temperature = 0.0;

    return *this;
}

/**
 * @brief Imports the surface properties from a file.
 * 
 * Reads the mu and sigma coefficients and autocorrelation length from the properties file
 * and stores them in the appropriate member variables.
 * 
 * @param properties_file The file containing surface properties.
 */
template <typename T>
void Surface<T>::import_properties(std::string properties_file) {

    unsigned int order, local_param_number;

    std::ifstream myfile;
    myfile.open(properties_file);
    myfile >> order;

    Matrix<T> mu_k(order, 1);
    Matrix<T> sigma_k(order, 1);

    // Read mu coefficients from the file.
    for(auto i = 0; i < order; i++) {
        myfile >> mu_k(i, 0);
    }

    // Read sigma coefficients from the file.
    for(auto i = 0; i < order; i++) {
        myfile >> sigma_k(i, 0);
    }

    // Read autocorrelation length from the file.
    myfile >> this->ac_length;

    this->mu_coefficients = mu_k;
    this->sigma_coefficients = sigma_k;

    myfile.close();
}
