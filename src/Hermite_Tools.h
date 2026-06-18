
/**
 * @brief Class for computing Hermite polynomials and related operations.
 * 
 * The Hermite class provides methods for evaluating Hermite polynomials, 
 * computing their derivatives, and importing data for roots and weights.
 * 
 * @tparam T The data type used for computations (e.g., float, double).
 */
template <typename T>
class Hermite {

    private:
        Matrix<T> roots; ///< Matrix storing the roots of Hermite polynomials.
        Matrix<T> weights; ///< Matrix storing the weights of Hermite polynomials.
        unsigned int max_order = 100; ///< Maximum order of Hermite polynomials supported.

    public:

        /**
         * @brief Default constructor for Hermite class.
         * 
         * Initializes the Hermite polynomial roots and weights by importing data.
         */
        Hermite();

        /**
         * @brief Destructor for Hermite class.
         */
        ~Hermite();

        /**
         * @brief Evaluates the Hermite polynomial at specified values.
         * 
         * @param coefficients Coefficients of the Hermite polynomial.
         * @param values Matrix of input values at which to evaluate the polynomial.
         * @return Matrix<T> Resulting values of the evaluated polynomial.
         */
        Matrix<T> evaluate(const Matrix<T>& coefficients, Matrix<T> values);

        /**
         * @brief Evaluates the Hermite polynomial at a single value.
         * 
         * @param coefficients Coefficients of the Hermite polynomial.
         * @param value The input value at which to evaluate the polynomial.
         * @return T The resulting value of the evaluated polynomial.
         */
        T evaluate(const Matrix<T>& coefficients, T value);

        /**
         * @brief Computes the derivative of the Hermite polynomial coefficients.
         * 
         * @param coefficients Coefficients of the Hermite polynomial.
         * @return Matrix<T> Coefficients of the derivative of the Hermite polynomial.
         */
        Matrix<T> coefficients_deriv(const Matrix<T>& coefficients);

        /**
         * @brief Evaluates the derivative of the Hermite polynomial at specified values.
         * 
         * @param coefficients Coefficients of the Hermite polynomial.
         * @param values Matrix of input values at which to evaluate the derivative.
         * @return Matrix<T> Resulting values of the evaluated derivative.
         */
        Matrix<T> evaluate_deriv(const Matrix<T>& coefficients, Matrix<T> values);

        /**
         * @brief Evaluates the derivative of the Hermite polynomial at a single value.
         * 
         * @param coefficients Coefficients of the Hermite polynomial.
         * @param value The input value at which to evaluate the derivative.
         * @return T The resulting value of the evaluated derivative.
         */
        T evaluate_deriv(Matrix<T>& coefficients, T value);

        /**
         * @brief Imports data for roots or weights from a file.
         * 
         * @param file The filename containing the data.
         * @return Matrix<T> The imported data as a matrix.
         */
        Matrix<T> import_data(const char* file);

        /**
         * @brief Retrieves the roots of the Hermite polynomial of specified order.
         * 
         * @param order The order of the Hermite polynomial.
         * @return Matrix<T> The roots of the Hermite polynomial.
         */
        Matrix<T> get_roots(unsigned int order);

        /**
         * @brief Retrieves the weights of the Hermite polynomial of specified order.
         * 
         * @param order The order of the Hermite polynomial.
         * @return Matrix<T> The weights of the Hermite polynomial.
         */
        Matrix<T> get_weights(unsigned int order);

        /**
         * @brief Computes the coefficients for the Hermite polynomial using a given function.
         * 
         * @param order The order of the Hermite polynomial.
         * @param func Pointer to the function used to compute the coefficients.
         * @return Matrix<T> The coefficients of the Hermite polynomial.
         */
        Matrix<T> compute_coefficients(unsigned int order, T (*func)(T));
};

/**
 * @brief Default constructor for Hermite class.
 * 
 * Initializes the Hermite polynomial by importing roots and weights from data files.
 */
template <typename T>
Hermite<T>::Hermite() {
    this->roots = import_data("code_files/roots.txt");
    this->weights = import_data("code_files/weights.txt");
    this->weights /= exp((this->roots ^ 2.0) * (-1.0));
}

/**
 * @brief Destructor for Hermite class.
 */
template <typename T>
Hermite<T>::~Hermite() {}

/**
 * @brief Evaluates the Hermite polynomial at specified matrix values.
 * 
 * @param coefficients Coefficients of the Hermite polynomial.
 * @param values Matrix of input values at which to evaluate the polynomial.
 * @return Matrix<T> Resulting values of the evaluated polynomial.
 */
template <typename T>
Matrix<T> Hermite<T>::evaluate(const Matrix<T>& coefficients, Matrix<T> values) {

    unsigned int order = coefficients.get_rows();
    Matrix<T> Hn_0(values), Hn_1(values), Hn(values), result(values);

    // Initialize the first two Hermite polynomials
    Hn_0 = 1.0;
    Hn_1 = values * 2.0;
    result = Hn_0 * coefficients(0, 0) + Hn_1 * coefficients(1, 0);

    // Compute higher order Hermite polynomials iteratively
    for (auto i = 2; i < order; i++) {
        Hn = Hn_1 * values * 2.0 - Hn_0 * 2.0 * (i - 1.0);
        Hn_0 = Hn_1;
        Hn_1 = Hn;
        result += (Hn * coefficients(i, 0));
    }

    return result;
}

/**
 * @brief Evaluates the Hermite polynomial at a single value.
 * 
 * @param coefficients Coefficients of the Hermite polynomial.
 * @param value The input value at which to evaluate the polynomial.
 * @return T The resulting value of the evaluated polynomial.
 */
template <typename T>
T Hermite<T>::evaluate(const Matrix<T>& coefficients, T value) {

    T Hn_0, Hn_1, Hn, result;
    unsigned int order = coefficients.get_rows();
    Hn_0 = 1.0;
    Hn_1 = value * 2.0;
    result = Hn_0 * coefficients(0, 0) + Hn_1 * coefficients(1, 0);

    // Compute higher order Hermite polynomials iteratively
    for (auto i = 2; i < order; i++) {
        Hn = Hn_1 * value * 2.0 - Hn_0 * 2.0 * (i - 1.0);
        Hn_0 = Hn_1;
        Hn_1 = Hn;
        result += Hn * coefficients(i, 0);
    }

    return result;
}

/**
 * @brief Computes the derivative of Hermite polynomial coefficients.
 * 
 * @param coefficients Coefficients of the Hermite polynomial.
 * @return Matrix<T> Coefficients of the derivative of the Hermite polynomial.
 */
template <typename T>
Matrix<T> Hermite<T>::coefficients_deriv(const Matrix<T>& coefficients) {

    unsigned int order = coefficients.get_rows();
    Matrix<T> new_coefficients(order - 1, 1);

    // Compute the derivative of the polynomial coefficients
    for (auto i = 1; i < order; i++) {
        new_coefficients(i - 1, 0) = 2.0 * i * coefficients(i, 0);
    }

    return new_coefficients;
}

/**
 * @brief Evaluates the derivative of the Hermite polynomial at specified matrix values.
 * 
 * This function first computes the derivative of the Hermite polynomial coefficients
 * and then evaluates the resulting derivative polynomial at the given input matrix values.
 * 
 * @param coefficients Coefficients of the Hermite polynomial.
 * @param values Matrix of input values at which to evaluate the derivative.
 * @return Matrix<T> The resulting values of the evaluated derivative.
 */
template <typename T>
Matrix<T> Hermite<T>::evaluate_deriv(const Matrix<T>& coefficients, Matrix<T> values) {
    Matrix<T> deriv_coefficients = coefficients_deriv(coefficients);
    return evaluate(deriv_coefficients, values);
}

/**
 * @brief Evaluates the derivative of the Hermite polynomial at a single value.
 * 
 * This function first computes the derivative of the Hermite polynomial coefficients
 * and then evaluates the resulting derivative polynomial at the given input value.
 * 
 * @param coefficients Coefficients of the Hermite polynomial.
 * @param value The input value at which to evaluate the derivative.
 * @return T The resulting value of the evaluated derivative.
 */
template <typename T>
T Hermite<T>::evaluate_deriv(Matrix<T>& coefficients, T value) {
    Matrix<T> deriv_coefficients = coefficients_deriv(coefficients);
    return evaluate(deriv_coefficients, value);
}

/**
 * @brief Imports data (e.g., roots or weights) from a file.
 * 
 * The imported data is expected to be formatted with the first two values in the file
 * representing the number of rows and columns, followed by the actual data values.
 * 
 * @param file The filename containing the data.
 * @return Matrix<T> The imported data as a matrix.
 */
template <typename T>
Matrix<T> Hermite<T>::import_data(const char* file) {
    unsigned int rows, cols;
    T buffer;
    std::ifstream myfile;
    myfile.open(file);
    myfile >> rows >> cols;

    Matrix<T> data(rows, cols);

    // Read data into the matrix
    for (auto i = 0; i < rows; i++) {
        for (auto j = 0; j < cols; j++) {
            myfile >> buffer;
            data(i, j) = buffer;
        }
    }
    myfile.close();
    return data;
}

/**
 * @brief Retrieves the roots of the Hermite polynomial of the specified order.
 * 
 * This function retrieves the roots from the precomputed data. If the order exceeds the maximum order, 
 * an error is printed, and an empty matrix is returned.
 * 
 * @param order The order of the Hermite polynomial.
 * @return Matrix<T> The roots of the Hermite polynomial for the specified order.
 */
template <typename T>
Matrix<T> Hermite<T>::get_roots(unsigned int order) {
    Matrix<T> roots_local(order, 1);
    if (order >= max_order) {
        std::cerr << "Maximum order exceeded!\n";
        return roots_local;
    } else {
        for (auto i = 0; i < order; i++) {
            roots_local(i, 0) = this->roots(order - 1, i);
        }
        return roots_local;
    }
}

/**
 * @brief Retrieves the weights of the Hermite polynomial of the specified order.
 * 
 * This function retrieves the weights from the precomputed data. If the order exceeds the maximum order, 
 * an error is printed, and an empty matrix is returned.
 * 
 * @param order The order of the Hermite polynomial.
 * @return Matrix<T> The weights of the Hermite polynomial for the specified order.
 */
template <typename T>
Matrix<T> Hermite<T>::get_weights(unsigned int order) {
    Matrix<T> weights_local(order, 1);
    if (order >= max_order) {
        std::cerr << "Maximum order exceeded!\n";
        return weights_local;
    } else {
        for (auto i = 0; i < order; i++) {
            weights_local(i, 0) = this->weights(order - 1, i);
        }
        return weights_local;
    }
}

/**
 * @brief Computes the coefficients of the Hermite polynomial for a given function.
 * 
 * This function computes the coefficients of the Hermite polynomial of the specified order
 * by evaluating the provided function at the roots of the polynomial and using the weights
 * for integration. The coefficients are calculated based on these values.
 * 
 * @param order The order of the Hermite polynomial.
 * @param func Pointer to the function to be used for computing the coefficients.
 * @return Matrix<T> The computed coefficients of the Hermite polynomial.
 */
template <typename T>
Matrix<T> Hermite<T>::compute_coefficients(unsigned int order, T (*func)(T)) {
    Matrix<T> weights = get_weights(order);
    Matrix<T> roots = get_roots(order);
    Matrix<T> func_values(order, 1);
    Matrix<T> coefficients(order, 1);

    T factorial = 1.0;

    // Evaluate the function at each root
    for (auto i = 0; i < order; i++) {
        func_values(i, 0) = *func(roots(i, 0));
    }

    // Compute the coefficients
    for (auto i = 0; i < order; i++) {
        coefficients(i, 0) = weights.tr().dot(func_values)(0, 0) / std::pow(2, i) / factorial / std::sqrt(M_PI);
        factorial *= i;
    }

    return coefficients;
}
