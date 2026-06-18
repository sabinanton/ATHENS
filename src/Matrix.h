/**
 * @file
 * @brief File containing the definition of the Matrix class.
 */



/**
 * @class Matrix
 * @brief This class implements a matrix object used for linear algebra and vectorized operations.
 *
 * The Matrix class provides a range of operations that can be performed on matrix objects, 
 * including addition, subtraction, multiplication, division, dot product, cross product, 
 * and various other operations common in linear algebra.
 */

template <typename T>
class Matrix {

    private:

        long rows, cols; ///< Number of rows and columns in the matrix
        T* data; ///< Pointer to the array holding the matrix's elements

    public:

        /**
         * @brief Default constructor. Initializes an empty matrix.
         */
        Matrix();

        /**
         * @brief Constructor that initializes the matrix with the specified dimensions.
         * 
         * The created matrix has the specified number of rows and columns, all set to zero.
         * 
         * @param dim_x Number of rows.
         * @param dim_y Number of columns.
         */
        Matrix(long dim_x, long dim_y);

        /**
         * @brief Constructor that initializes the matrix with the specified dimensions and data.
         * 
         * The created matrix has the specified number of rows and columns. The data for the matrix
         * is provided by an initializer list. If the size of the list does not match the expected size
         * of the matrix (rows * columns), an error message is displayed.
         * 
         * @param dim_x Number of rows.
         * @param dim_y Number of columns.
         * @param l Initializer list containing data to populate the matrix.
         */
        Matrix(long dim_x, long dim_y, std::initializer_list<T> l);

        /**
         * @brief Destructor. Cleans up the allocated memory.
         */
        ~Matrix();

        /**
         * @brief Copy constructor. Creates a copy of an existing Matrix object.
         * 
         * This constructor creates a new Matrix object with the same dimensions 
         * and copies all elements from the original matrix to the new one.
         * 
         * @param other The Matrix object to copy from.
         */
        Matrix(const Matrix& other);

        /**
         * @brief Move constructor. Moves the data from a given Matrix object to the new one.
         * 
         * This constructor "steals" the internal data from the original object and 
         * leaves it in an empty state. This is more efficient than copying when the 
         * original object is no longer needed.
         * 
         * @param other The Matrix object whose data is to be moved.
         */
        Matrix(Matrix&& other);

        T* getData() {return this->data;}
        long get_rows() const {return rows;}
        long get_cols() const {return cols;}

        /**
         * @brief Copy assignment operator. Copies the data from the given Matrix object.
         * 
         * This operator replaces the current Matrix's data with a copy of the data from the 'other' Matrix.
         * 
         * @param other The Matrix object to copy data from.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator=(const Matrix& other);

        /**
         * @brief Move assignment operator. Moves the data from the given Matrix object.
         * 
         * This operator "steals" the internal data from the 'other' Matrix, leaving it in an empty state.
         * 
         * @param other The Matrix object to move data from.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator=(Matrix&& other);

        /**
         * @brief Constant assignment operator. Assigns all Matrix values to one constant.
         * 
         * @param constant The constant to be assigned.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator=(T constant);

        /**
         * @brief Compound assignment operator (addition). Adds another Matrix to this Matrix.
         * 
         * This operator performs element-wise addition between this Matrix and the 'other' Matrix. 
         * Both matrices must have the same dimensions, otherwise, an exception is thrown.
         * 
         * @param other The Matrix object to be added to this one.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator+=(const Matrix& other);

        /**
         * @brief Compound assignment operator (subtraction). Subtracts another Matrix from this Matrix.
         * 
         * This operator performs element-wise subtraction of the 'other' Matrix from this Matrix. 
         * Both matrices must have the same dimensions, otherwise, an exception is thrown.
         * 
         * @param other The Matrix object to be subtracted from this one.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator-=(const Matrix& other);

         /**
         * @brief Compound assignment operator (multiplication). Multiplies this Matrix with another Matrix.
         * 
         * This operator performs matrix multiplication between this Matrix and the 'other' Matrix. 
         * The number of columns in this Matrix must be equal to the number of rows in the 'other' matrix, 
         * otherwise, an exception is thrown.
         * 
         * @param other The Matrix object to be multiplied with this one.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator*=(const Matrix& other);

        /**
         * @brief Compound assignment operator (division). Divides this Matrix by another Matrix.
         * 
         * This operator performs element-wise division of this Matrix by the 'other' Matrix. 
         * Both matrices must have the same dimensions, and the 'other' matrix must not contain any zero elements,
         * otherwise, an exception is thrown.
         * 
         * @param other The Matrix object to divide this one by.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator/=(const Matrix& other);

        /**
         * @brief Compound assignment operator (addition). Adds a scalar to this Matrix.
         * 
         * This operator performs element-wise addition of a scalar value to this Matrix. 
         * The scalar value is added to every element of this Matrix.
         * 
         * @param other The scalar value to be added to this Matrix's elements.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator+=(T other);

        /**
         * @brief Compound assignment operator (subtraction). Subtracts a scalar from this Matrix.
         * 
         * This operator performs element-wise subtraction of a scalar value from this Matrix. 
         * The scalar value is subtracted from every element of this Matrix.
         * 
         * @param other The scalar value to be subtracted from this Matrix's elements.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator-=(T other);

        /**
         * @brief Compound assignment operator (multiplication). Multiplies this Matrix by a scalar.
         * 
         * This operator performs element-wise multiplication of this Matrix by a scalar value. 
         * Each element of this Matrix is multiplied by the scalar value.
         * 
         * @param other The scalar value to multiply this Matrix's elements by.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator*=(T other);

        /**
         * @brief Compound assignment operator (division). Divides this Matrix by a scalar.
         * 
         * This operator performs element-wise division of this Matrix by a scalar value. 
         * Each element of this Matrix is divided by the scalar value. The scalar must not be zero,
         * otherwise, an exception is thrown.
         * 
         * @param other The scalar value to divide this Matrix's elements by.
         * @return A reference to the current object, to enable operator chaining.
         */
        Matrix& operator/=(T other);


        /**
         * @brief Subscript operator for element access (read/write).
         * 
         * This operator provides direct access to the matrix data using a two-dimensional index.
         * It performs boundary checks and throws an exception if the indices are out of range.
         * 
         * @param i The row index.
         * @param j The column index.
         * @return A reference to the element at the specified position.
         */
        T& operator()( const long i, const long j ) { return data[i*cols + j]; }

        /**
         * @brief Subscript operator for element access (read-only).
         * 
         * This operator provides direct access to the matrix data using a two-dimensional index.
         * It performs boundary checks and throws an exception if the indices are out of range.
         * 
         * @param i The row index.
         * @param j The column index.
         * @return A constant reference to the element at the specified position.
         */
        const T& operator()( const long i, const long j ) const { return data[i*cols + j]; }

         /**
         * @brief Subscript operator for element access in a single row (read/write).
         * 
         * This operator provides direct access to the matrix data considering it as a single row.
         * It performs boundary checks and throws an exception if the index is out of range.
         * 
         * @param n The index in the matrix considering it as a single row.
         * @return A reference to the element at the specified position.
         */
        T& operator()( const long n ) { return data[n]; }

        /**
         * @brief Subscript operator for element access in a single row (read-only).
         * 
         * This operator provides direct access to the matrix data considering it as a single row.
         * It performs boundary checks and throws an exception if the index is out of range.
         * 
         * @param n The index in the matrix considering it as a single row.
         * @return A constant reference to the element at the specified position.
         */
        const T& operator()( const long n ) const { return data[n]; }

        Matrix operator+(const Matrix& other);
        Matrix operator-(const Matrix& other);
        Matrix operator*(const Matrix& other);
        Matrix operator/(const Matrix& other);

        Matrix operator*(T constant);
        Matrix operator/(T constant);
        Matrix operator^(T constant);
        Matrix operator+(T constant);
        Matrix operator-(T constant);

        /**
         * @brief Calculates the dot product of two matrices.
         * 
         * This function calculates the dot product, or scalar product, of two matrices.
         * Both matrices must have the same dimensions, otherwise, an exception is thrown.
         * 
         * @param other The matrix with which to calculate the dot product.
         * @return A new Matrix object containing the result of the dot product.
         */
        Matrix dot(const Matrix& other);

        /**
         * @brief Calculates the cross product of two matrices.
         * 
         * This function calculates the cross product, or vector product, of two matrices.
         * Both matrices must be 3x1 (three rows, one column), otherwise, an exception is thrown.
         * 
         * @param other The matrix with which to calculate the cross product.
         * @return A new Matrix object containing the result of the cross product.
         */
        Matrix cross(const Matrix& other);

        /**
         * @brief Calculates the outer product of two matrices.
         * 
         * This function calculates the outer product of two 1D vectors.
         * Both vectors must be column vectors, otherwise, an exception is thrown.
         * 
         * @param other The matrix with which to calculate the outer product.
         * @return A new Matrix object containing the result of the outer product.
         */
        Matrix outer(const Matrix& other);

        /**
         * @brief Calculates the meshgrid of 2 1D matrices.
         * 
         * This function calculates the meshgrid matrix of two 1D vectors.
         * Both vectors must be column vectors, otherwise, an exception is thrown.
         * 
         * @param other The matrix with which to calculate the meshgrid.
         * @return A new Matrix object containing the result of the mehsgrid.
         */
        std::tuple<Matrix<T>, Matrix<T>> meshgrid(const Matrix& other);

        /**
         * @brief Transposes the matrix.
         * 
         * This function returns a new matrix which is the transpose of the current matrix.
         * 
         * @return A new Matrix object that is the transpose of the current matrix.
         */
        Matrix tr();

        /**
         * @brief Calculates the Frobenius norm of the matrix.
         * 
         * This function calculates the square root of the sum of the absolute squares of its elements.
         * 
         * @return The Frobenius norm of the matrix.
         */
        T norm();

        /**
         * @brief Finds the maximum element in the matrix.
         * 
         * This function searches for the largest element in the matrix.
         * 
         * @return The value of the largest element in the matrix.
         */
        T max();

        /**
         * @brief Finds the minimum element in the matrix.
         * 
         * This function searches for the smallest element in the matrix.
         * 
         * @return The value of the smallest element in the matrix.
         */
        T min();

        /**
         * @brief Prints the matrix to the console.
         * 
         * This function prints a textual representation of the matrix to the console.
         * Each row is printed on a new line, and the elements are separated by spaces.
         */
        void print();

};

template <typename T>
Matrix<T>::Matrix() : rows(0), cols(0), data(nullptr) {}

template <typename T>
Matrix<T>::Matrix(long dim_x, long dim_y) : rows(dim_x), cols(dim_y), data(new T[rows * cols]) {
    for(long i = 0; i < rows * cols; i ++) {
        data[i] = 0.0;
    }
}

template <typename T>
Matrix<T>::Matrix(long dim_x, long dim_y, std::initializer_list<T> l) : rows(dim_x), cols(dim_y), data(new T[rows * cols]) {
    
    if(rows * cols != l.size()) {
        std::cout<<"The size of the list ("<<l.size()<<") does not match the matrix shape ("<<dim_x<<", "<<dim_y<<")\n";
    }
    else{
        std::uninitialized_copy(l.begin(), l.end(), data);
    }
}

template<typename T>
Matrix<T>::~Matrix() {
    this->rows = 0;
    this->cols = 0;
    delete[] data;
    data = nullptr;
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other) {
    if(this != &other) {
        this->rows = other.rows;
        this->cols = other.cols;
        this->data = new T[rows * cols];

        for(long i = 0; i < rows * cols; i ++){
            this->data[i] = other.data[i];
        }
    }
}

template <typename T>
Matrix<T>::Matrix(Matrix<T>&& other) {
    if(this != &other) { 
        this->rows = other.rows;
        this->cols = other.cols;

        this->data = other.data;
        other.data = nullptr;
    }
    
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    if(this != &other) {
        delete[] data;
        this->rows = other.rows;
        this->cols = other.cols;
        data = new T[rows * cols];

        for(long i = 0; i < rows * cols; i ++) {
            this->data[i] = other.data[i];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& other) {
    if(this != &other) {
        delete[] data;
        this->rows = other.rows;
        this->cols = other.cols;
        this->data = other.data;
        other.rows = 0;
        other.cols = 0;
        other.data = nullptr;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(T constant) {
    for(long i = 0; i < rows * cols; i ++) {
            this->data[i] = constant;
        }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
     if(this->rows != other.rows || this->cols != other.cols) {
        std::cout<<"Matrix dimensions do not match!\n";
    }
    else{
        for(long i = 0; i < rows * cols; i ++){
            this->data[i] += other.data[i];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
     if(this->rows != other.rows || this->cols != other.cols) {
        std::cout<<"Matrix dimensions do not match!\n";
    }
    else{
        for(long i = 0; i < rows * cols; i ++){
            this->data[i] -= other.data[i];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
     if(this->rows != other.rows || this->cols != other.cols) {
        std::cout<<"Matrix dimensions do not match!\n";
    }
    else{
        for(long i = 0; i < rows * cols; i ++){
            this->data[i] *= other.data[i];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& other) {
     if(this->rows != other.rows || this->cols != other.cols) {
        std::cout<<"Matrix dimensions do not match!\n";
    }
    else{
        for(long i = 0; i < rows * cols; i ++){
            this->data[i] /= other.data[i];
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(T other) {
    for(long i = 0; i < rows * cols; i ++){
        data[i] *= other;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(T other) {
    for(long i = 0; i < rows * cols; i ++){
        data[i] /= other;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(T other) {
    for(long i = 0; i < rows * cols; i ++){
        data[i] += other;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(T other) {
    for(long i = 0; i < rows * cols; i ++){
        data[i] -= other;
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) {
   Matrix<T> res(*this);
   res += other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(T other) {
   Matrix<T> res(*this);
   res += other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) {
   Matrix<T> res(*this);
   res -= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(T other) {
   Matrix<T> res(*this);
   res -= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) {
   Matrix<T> res(*this);
   res *= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(T other) {
   Matrix<T> res(*this);
   res *= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(T other) {
   Matrix<T> res(*this);
   res /= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const Matrix<T>& other) {
   Matrix<T> res(*this);
   res /= other;
   return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator^(T constant) {
    Matrix<T> res(*this);
    for(long i = 0; i < rows; i ++){
        for(long j = 0; j < cols; j ++) {
            res(i, j) = pow(this->data[i * cols + j], constant);
        }
    }
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::dot(const Matrix<T>& other) {

    Matrix<T> res(this->rows, other.cols);

    if(this->cols != other.rows) {
        std::cout<<"Matrix dimensions of ("<<this->rows<<", "<<this->cols<<") and ("<<other.rows<<", "<<other.cols<<") are not compatible!\n";
    }
    else{
        for(long i = 0; i < this->rows; i ++) {
            for(long j = 0; j < other.cols; j ++) {
                for(long k = 0; k < this->cols; k ++) {
                    res.data[i * other.cols + j] += this->data[i * this->cols + k] * other.data[k * other.cols + j];
                }
            }
        }
    }
    return res;
}

template <typename T> 
Matrix<T> Matrix<T>::cross(const Matrix<T>& other) {

    Matrix<T> result(3, 1);
    if(this->rows != 3 || this->cols != 1 || other.rows != 3 || other.cols != 1) {
        std::cerr<<"Dimensions incompatible with the cross product operation!\n";
    }
    else{
        result(0, 0) = this->data[1 * cols + 0] * other(2, 0) - this->data[2 * cols + 0] * other(1, 0);
        result(1, 0) = this->data[2 * cols + 0] * other(0, 0) - this->data[0 * cols + 0] * other(2, 0);
        result(2, 0) = this->data[0 * cols + 0] * other(1, 0) - this->data[1 * cols + 0] * other(0, 0);
    }
    return result;
}

template <typename T> 
Matrix<T> Matrix<T>::outer(const Matrix<T>& other) {

    Matrix<T> result(this->rows, other.rows);
    if(this->cols != 1 || other.cols != 1) {
        std::cerr<<"Dimensions incompatible with the outer product operation!\n";
    }
    else{
        for(long i = 0; i < this->rows; i ++) {
            for(long j = 0; j < other.rows; j ++) {
                result(i, j) = this->data[j] * other(i);
            }
        }
    }
    return result;
}

template <typename T> 
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::meshgrid(const Matrix<T>& other) {

    Matrix<T> mat1(this->rows, other.rows), mat2(this->rows, other.rows);
    if(this->cols != 1 || other.cols != 1) {
        std::cerr<<"Dimensions incompatible with the outer product operation!\n";
    }
    else{
        for(long i = 0; i < this->rows; i ++) {
            for(long j = 0; j < other.rows; j ++) {
                mat1(i, j) = this->data[i];
                mat2(i, j) = other(j);
            }
        }
    }
    return std::make_tuple(mat1, mat2);
}

template <typename T>
T Matrix<T>::norm() {
    T res = 0.0;

    for(long i = 0; i < rows * cols; i ++) {
        res += this->data[i] * this->data[i];
    }

    return sqrt(res);
}

template <typename T> 
T Matrix<T>::max() {
    T max_value = -1e13;

    for(auto i = 0; i < rows * cols; i ++) {
        if(max_value < data[i]) max_value = data[i];
    }
    return max_value;
}

template <typename T> 
T Matrix<T>::min() {
    T min_value = 1e13;

    for(auto i = 0; i < rows * cols; i ++) {
        if(min_value > data[i]) min_value = data[i];
    }
    return min_value;
}

template <typename T>
Matrix<T> Matrix<T>::tr() {
    Matrix<T> new_matrix(cols, rows);

    for(auto i = 0; i < rows; i ++) {
        for(auto j = 0; j < cols; j ++) {
            new_matrix(j, i) = data[i * cols + j];
        }
    }

    return new_matrix;
}

template <typename T>
Matrix<T> cos(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::cos(other(i));
    }
    return result;
}

template <typename T>
Matrix<T> sin(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::sin(other(i));
    }
    return result;
}

template <typename T>
Matrix<T> tan(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::tan(other(i));
    }
    return result;
}

template <typename T>
Matrix<T> exp(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::exp(other(i));
    }
    return result;
}

template <typename T>
Matrix<T> erf(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::erf(other(i));
    }
    return result;
}

template <typename T>
Matrix<T> erfc(const Matrix<T>& other) {

    Matrix<T> result(other.get_rows(), other.get_cols());

    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        result(i) = std::erfc(other(i));
    }
    return result;
}

template <typename T>
T sum(const Matrix<T>& other) {
    T res = 0.0;
    for(long i = 0; i < other.get_rows() * other.get_cols(); i ++){
        res += other(i);
    }
    return res;
}

template <typename T>
void Matrix<T>::print() {
    for(auto i = 0; i < rows; i ++) {
        for(auto j = 0; j < cols; j ++) {
            std::cout<<data[i * cols + j] <<" ";
        }
        std::cout<<"\n";
    }
}