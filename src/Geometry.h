/**
 * @brief Class representing a geometric structure composed of vertices, normals, and properties.
 * 
 * This class is used to manage 3D geometries, typically loaded from STL files. It stores vertices, normals, 
 * gas-surface interaction (GSI) properties, and surface properties. The class provides methods for manipulating 
 * the geometry (shift, scale, rotate), adding properties, and calculating bounding boxes and areas.
 * 
 * @tparam T Type of the numeric values used for geometry (e.g., float, double).
 */
template <typename T>
class Geometry {

    private:

        unsigned int num_vertices; ///< Number of vertices (triangles) in the geometry
        char header[81];           ///< Header information from the STL file
        bool record = true;        ///< A flag allowing the recording of geometry transformations
        Matrix<T> p1[4], p2[4], p3[4]; ///< Matrices for storing vertex points p1, p2, p3
        Matrix<T> norm_x, norm_y, norm_z, norm_mod; ///< Matrices for storing the normal vectors and their magnitudes
        T A_ref, c_ref;            ///< Reference area and reference speed for the geometry
        std::vector<Matrix<T>> GSI_prop; ///< Vector of gas-surface interaction (GSI) properties
        std::vector<Matrix<T>> surface_prop; ///< Vector of surface properties
        std::vector<const char*> GSI_prop_names; ///< Names of the GSI properties
        std::vector<const char*> surface_prop_names; ///< Names of the surface properties
        Matrix<T> bounding_box;   ///< Bounding box matrix for the geometry
        std::vector<int> rotation_axis_history;   ///< Array containing the history of rotation axes
        std::vector<T> rotation_history;   ///< Array containing the history of rotations
        std::vector<Matrix<T>> shift_history;   ///< Array containing the history of translations
        std::vector<T> scale_history;   ///< Array containing the history of scalings

    public:

        /**
         * @brief Default constructor for Geometry class.
         * 
         * Initializes the geometry with zero vertices, and default values for reference area and speed.
         */
        Geometry();

        /**
         * @brief Constructor that loads geometry data from a binary STL file.
         * 
         * @param filename The name of the STL file to load.
         */
        Geometry(std::string filename);

        /**
         * @brief Destructor for Geometry class.
         * 
         * Resets the number of vertices to zero.
         */
        ~Geometry();

        /**
         * @brief Get the number of vertices (triangles) in the geometry.
         * 
         * @return The number of vertices.
         */
        unsigned int get_num_vertices() const { return num_vertices; }

        /**
         * @brief Shift the geometry by an offset.
         * 
         * @param offset Array containing x, y, z offset values.
         */
        void shift(T* offset);

        /**
         * @brief Shift the geometry by a Matrix offset.
         * 
         * @param offset Matrix representing x, y, z offset values.
         */
        void shift(Matrix<T> offset);

        /**
         * @brief Scale the geometry by a factor.
         * 
         * @param factor The scaling factor to apply.
         */
        void scale(T factor);

        /**
         * @brief Rotate the geometry around a specific axis.
         * 
         * @param axis Axis of rotation (0 for x, 1 for y, 2 for z).
         * @param angle Angle of rotation in radians.
         */
        void rotate(int axis, T angle);

        /**
         * @brief Reset the orientation, position and scale of the geometry.
         */
        void reset();

        /**
         * @brief Add a scalar GSI property to the geometry.
         * 
         * @param property_name The name of the GSI property.
         * @param value The scalar value to assign to the property.
         */
        void add_GSI_property(const char* property_name, T value);

        /**
         * @brief Add a matrix GSI property to the geometry.
         * 
         * @param property_name The name of the GSI property.
         * @param values Matrix of values for the GSI property.
         */
        void add_GSI_property(const char* property_name, const Matrix<T>& values);

        /**
         * @brief Add a scalar surface property to the geometry.
         * 
         * @param property_name The name of the surface property.
         * @param value The scalar value to assign to the property.
         */
        void add_surface_property(const char* property_name, T value);

        /**
         * @brief Add a matrix surface property to the geometry.
         * 
         * @param property_name The name of the surface property.
         * @param values Matrix of values for the surface property.
         */
        void add_surface_property(const char* property_name, const Matrix<T>& values);

        /**
         * @brief Set the reference area for the geometry.
         * 
         * @param A_ref The reference area to set.
         */
        void set_A_ref(T A_ref) { this->A_ref = A_ref; }

        /**
         * @brief Set the reference speed for the geometry.
         * 
         * @param c_ref The reference speed to set.
         */
        void set_c_ref(T c_ref) { this->c_ref = c_ref; }

        /**
         * @brief Compute the bounding box for the geometry.
         */
        void compute_bounding_box();

        /**
         * @brief Get the reference area.
         * 
         * @return The reference area.
         */
        T get_A_ref() { return A_ref; }

        /**
         * @brief Get the reference speed.
         * 
         * @return The reference speed.
         */
        T get_c_ref() { return c_ref; }

        /**
         * @brief Compute the area of a specific triangle in the geometry.
         * 
         * @param index The index of the triangle.
         * @return The computed area of the triangle.
         */
        T compute_area(unsigned int index);

        /**
         * @brief Print the geometry details, including vertices, normals, and properties.
         */
        void print();

        // Getter methods for the vertex coordinates and normal vectors.
        Matrix<T>& get_p1_x() { return this->p1[0]; }
        Matrix<T>& get_p1_y() { return this->p1[1]; }
        Matrix<T>& get_p1_z() { return this->p1[2]; }
        Matrix<T>& get_p2_x() { return this->p2[0]; }
        Matrix<T>& get_p2_y() { return this->p2[1]; }
        Matrix<T>& get_p2_z() { return this->p2[2]; }
        Matrix<T>& get_p3_x() { return this->p3[0]; }
        Matrix<T>& get_p3_y() { return this->p3[1]; }
        Matrix<T>& get_p3_z() { return this->p3[2]; }
        Matrix<T>& get_norm_x() { return this->norm_x; }
        Matrix<T>& get_norm_y() { return this->norm_y; }
        Matrix<T>& get_norm_z() { return this->norm_z; }

        /**
         * @brief Get a GSI property by index.
         * 
         * @param index Index of the GSI property.
         * @return The GSI property matrix.
         */
        Matrix<T>& get_GSI_property(int index) { return GSI_prop[index]; }

        /**
         * @brief Get the name of a GSI property by index.
         * 
         * @param index Index of the GSI property.
         * @return The name of the GSI property.
         */
        char* get_GSI_property_name(int index) { return GSI_prop_names[index]; }

        /**
         * @brief Get the bounding box of the geometry.
         * 
         * @return The bounding box matrix.
         */
        Matrix<T> get_bounding_box();

        /**
         * @brief Get the center of the bounding box of the geometry.
         * 
         * @return A matrix representing the center of the bounding box.
         */
        Matrix<T> get_center();

        /**
         * @brief Get all GSI properties.
         * 
         * @return A vector of GSI property matrices.
         */
        std::vector<Matrix<T>>& get_GSI_properties() { return GSI_prop; }

        /**
         * @brief Get all GSI property names.
         * 
         * @return A vector of GSI property names.
         */
        std::vector<const char*>& get_GSI_property_names() { return GSI_prop_names; }

        /**
         * @brief Get all surface properties.
         * 
         * @return A vector of surface property matrices.
         */
        std::vector<Matrix<T>>& get_surface_properties() { return surface_prop; }

        /**
         * @brief Get all surface property names.
         * 
         * @return A vector of surface property names.
         */
        std::vector<const char*>& get_surface_property_names() { return surface_prop_names; }
};

/**
 * @brief Opens a binary file for reading.
 * 
 * @param filename The name of the file to open.
 * @return A binary input file stream.
 */
std::ifstream open_binary_file(const char* filename) {
    return std::ifstream{filename, std::ifstream::in | std::ifstream::binary};
}

/**
 * @brief Casts a pointer to its equivalent character pointer.
 * 
 * @tparam T Type of the pointer.
 * @param pointer The pointer to cast.
 * @return The pointer cast to char*.
 */
template <typename T>
char* as_char_ptr(T* pointer) {
    return reinterpret_cast<char*>(pointer);
}

/**
 * @brief Reads a single binary value from an input stream.
 * 
 * @tparam T Type of the value to be read.
 * @param in Input file stream from which the value is read.
 * @param dst Pointer to the destination where the read value will be stored.
 */
template <typename T>
void read_binary_value(std::ifstream& in, T* dst) {
    in.read(as_char_ptr(dst), sizeof(T));
}

/**
 * @brief Reads an array of binary values from an input stream.
 * 
 * @tparam T Type of the values in the array.
 * @param in Input file stream from which the values are read.
 * @param dst Pointer to the destination array where the read values will be stored.
 * @param array_length Number of elements in the array.
 */
template <typename T>
void read_binary_array(std::ifstream& in, T* dst, size_t array_length) {
    size_t n_bytes = array_length * sizeof(T);
    in.read(as_char_ptr(dst), n_bytes);
}

/**
 * @brief Default constructor for the Geometry class.
 * 
 * Initializes the number of vertices to zero and sets default values for the reference area and speed.
 */
template <typename T>
Geometry<T>::Geometry() : num_vertices(0), A_ref(1.0), c_ref(1.0), bounding_box(6, 1) {}

/**
 * @brief Constructor that loads geometry data from a binary STL file.
 * 
 * This constructor reads an STL file, initializes the vertices, normals, and other properties,
 * and computes the bounding box for the geometry.
 * 
 * @param filename Name of the binary STL file to load.
 */
template <typename T>
Geometry<T>::Geometry(std::string filename) : A_ref(1.0), c_ref(1.0) {

    float p1_local[4], p2_local[4], p3_local[4], normal[4]; // Temporary arrays to hold vertex and normal data
    char dummy[3]; // Dummy variable to skip unused bytes in the STL format

    // Open the binary STL file
    auto in = open_binary_file(filename.c_str());
    
    // Read the STL header and the number of vertices (triangles)
    read_binary_array<char>(in, header, 80);
    header[80] = '\0'; // Null-terminate the header string
    read_binary_value<unsigned int>(in, &num_vertices);

    // Initialize matrices for vertices and normals based on the number of triangles
    p1[0] = Matrix<T>(1, num_vertices);
    p1[1] = Matrix<T>(1, num_vertices);
    p1[2] = Matrix<T>(1, num_vertices);
    p2[0] = Matrix<T>(1, num_vertices);
    p2[1] = Matrix<T>(1, num_vertices);
    p2[2] = Matrix<T>(1, num_vertices);
    p3[0] = Matrix<T>(1, num_vertices);
    p3[1] = Matrix<T>(1, num_vertices);
    p3[2] = Matrix<T>(1, num_vertices);

    norm_x = Matrix<T>(1, num_vertices);
    norm_y = Matrix<T>(1, num_vertices);
    norm_z = Matrix<T>(1, num_vertices);
    norm_mod = Matrix<T>(1, num_vertices);

    // Read the vertex and normal data for each triangle
    for (auto i = 0; i < num_vertices; i++) {
        read_binary_array<float>(in, normal, 3);
        read_binary_array<float>(in, p1_local, 3);
        read_binary_array<float>(in, p2_local, 3);
        read_binary_array<float>(in, p3_local, 3);
        in.read(dummy, 2); // Skip two bytes (STL format padding)
        
        // Store the vertex and normal data in the respective matrices
        p1[0](0, i) = static_cast<T>(p1_local[0]);
        p1[1](0, i) = static_cast<T>(p1_local[1]);
        p1[2](0, i) = static_cast<T>(p1_local[2]);
        p2[0](0, i) = static_cast<T>(p2_local[0]);
        p2[1](0, i) = static_cast<T>(p2_local[1]);
        p2[2](0, i) = static_cast<T>(p2_local[2]);
        p3[0](0, i) = static_cast<T>(p3_local[0]);
        p3[1](0, i) = static_cast<T>(p3_local[1]);
        p3[2](0, i) = static_cast<T>(p3_local[2]);

        norm_x(0, i) = static_cast<T>(normal[0]);
        norm_y(0, i) = static_cast<T>(normal[1]);
        norm_z(0, i) = static_cast<T>(normal[2]);

        // Calculate the magnitude of the normal vector and normalize it
        norm_mod(0, i) = std::sqrt(norm_x(0, i) * norm_x(0, i) + norm_y(0, i) * norm_y(0, i) + norm_z(0, i) * norm_z(0, i));
        norm_x(0, i) /= norm_mod(0, i);
        norm_y(0, i) /= norm_mod(0, i);
        norm_z(0, i) /= norm_mod(0, i);
    }

    // Compute the bounding box for the geometry
    compute_bounding_box();
}

/**
 * @brief Destructor for the Geometry class.
 * 
 * Resets the number of vertices to zero, effectively clearing the geometry data.
 */
template <typename T>
Geometry<T>::~Geometry() {
    num_vertices = 0;
}

/**
 * @brief Compute the bounding box of the geometry.
 * 
 * This method calculates the minimum and maximum coordinates of the vertices
 * in all three dimensions and stores them in the bounding box matrix.
 */
template <typename T>
void Geometry<T>::compute_bounding_box() {
    this->bounding_box = Matrix<T>(6, 1);
    T bb_x_1, bb_x_2, bb_y_1, bb_y_2, bb_z_1, bb_z_2;

    // Compute the min and max values for each dimension (x, y, z)
    bb_x_1 = std::min(p1[0].min(), std::min(p2[0].min(), p3[0].min()));
    bb_x_2 = std::max(p1[0].max(), std::max(p2[0].max(), p3[0].max()));
    bb_y_1 = std::min(p1[1].min(), std::min(p2[1].min(), p3[1].min()));
    bb_y_2 = std::max(p1[1].max(), std::max(p2[1].max(), p3[1].max()));
    bb_z_1 = std::min(p1[2].min(), std::min(p2[2].min(), p3[2].min()));
    bb_z_2 = std::max(p1[2].max(), std::max(p2[2].max(), p3[2].max()));

    // Store the bounding box values
    bounding_box(0, 0) = bb_x_1;
    bounding_box(1, 0) = bb_x_2;
    bounding_box(2, 0) = bb_y_1;
    bounding_box(3, 0) = bb_y_2;
    bounding_box(4, 0) = bb_z_1;
    bounding_box(5, 0) = bb_z_2;
}

/**
 * @brief Get the center point of the geometry's bounding box.
 * 
 * @return Matrix representing the center coordinates of the bounding box.
 */
template <typename T>
Matrix<T> Geometry<T>::get_center() {
    Matrix<T> center(3, 1);

    // Calculate the center point by averaging the min and max coordinates in each dimension
    center(0) = (bounding_box(0, 0) + bounding_box(1, 0)) / 2.0;
    center(1) = (bounding_box(3, 0) + bounding_box(2, 0)) / 2.0;
    center(2) = (bounding_box(5, 0) + bounding_box(4, 0)) / 2.0;

    return center;
}

/**
 * @brief Shift the geometry by an offset (array version).
 * 
 * This method shifts the positions of all vertices (p1, p2, p3) by a specified x, y, z offset.
 * It also updates the bounding box of the geometry.
 * 
 * @param offset Array containing the offset values for x, y, and z.
 */
template <typename T>
void Geometry<T>::shift(T* offset) {
    // Shift all vertices by the specified offset in x, y, and z directions.
    p1[0] += offset[0];
    p1[1] += offset[1];
    p1[2] += offset[2];

    p2[0] += offset[0];
    p2[1] += offset[1];
    p2[2] += offset[2];

    p3[0] += offset[0];
    p3[1] += offset[1];
    p3[2] += offset[2];

    // Update bounding box by adding offset to each corresponding dimension.
    bounding_box(0, 0) += offset[0];
    bounding_box(1, 0) += offset[0];
    bounding_box(2, 0) += offset[1];
    bounding_box(3, 0) += offset[1];
    bounding_box(4, 0) += offset[2];
    bounding_box(5, 0) += offset[2];

    // Add the current shift to the shift history
    if (record) {
        Matrix<T> new_shift(3, 1, {offset[0], offset[1], offset[2]});
        this->shift_history.push_back(new_shift);
    }
}

/**
 * @brief Shift the geometry by an offset (Matrix version).
 * 
 * This method shifts the positions of all vertices (p1, p2, p3) using a Matrix that contains x, y, z offsets.
 * It also updates the bounding box of the geometry.
 * 
 * @param offset Matrix containing the offset values for x, y, and z.
 */
template <typename T>
void Geometry<T>::shift(Matrix<T> offset) {
    // Shift all vertices by the offset Matrix values in x, y, and z directions.
    p1[0] += offset(0);
    p1[1] += offset(1);
    p1[2] += offset(2);

    p2[0] += offset(0);
    p2[1] += offset(1);
    p2[2] += offset(2);

    p3[0] += offset(0);
    p3[1] += offset(1);
    p3[2] += offset(2);

    // Update bounding box by adding offset to each corresponding dimension.
    bounding_box(0, 0) += offset(0);
    bounding_box(1, 0) += offset(0);
    bounding_box(2, 0) += offset(1);
    bounding_box(3, 0) += offset(1);
    bounding_box(4, 0) += offset(2);
    bounding_box(5, 0) += offset(2);

    // Add the current shift to the shift history
    if (record) {
        this->shift_history.push_back(offset);
    }
}

/**
 * @brief Scale the geometry by a factor.
 * 
 * This method scales all vertices (p1, p2, p3) and the bounding box of the geometry by the given factor.
 * 
 * @param factor The scaling factor to apply to the geometry.
 */
template <typename T>
void Geometry<T>::scale(T factor) {
    // Scale all vertices by the given factor.
    p1[0] *= factor;
    p1[1] *= factor;
    p1[2] *= factor;

    p2[0] *= factor;
    p2[1] *= factor;
    p2[2] *= factor;

    p3[0] *= factor;
    p3[1] *= factor;
    p3[2] *= factor;

    // Scale the bounding box by the given factor.
    bounding_box *= factor;

    // Add the current scale factor to the scale history
    if (record) {
        this->scale_history.push_back(factor);
    }
}

/**
 * @brief Rotates the geometry around the specified axis by the given angle.
 * 
 * This function rotates the vertices and normals of the geometry based on the given axis (0: x-axis, 1: y-axis, 2: z-axis)
 * and angle (in radians). After the rotation, the bounding box is recomputed.
 * 
 * @param axis The axis of rotation (0: x-axis, 1: y-axis, 2: z-axis).
 * @param angle The rotation angle (in radians).
 */
template <typename T>
void Geometry<T>::rotate(int axis, T angle) {
    if (axis == 0) {
        // Rotation around the x-axis
        Matrix<T> p1_y_new = (p1[1] * cos(angle) + p1[2] * sin(angle));
        Matrix<T> p1_z_new = (p1[2] * cos(angle) - p1[1] * sin(angle));
        Matrix<T> p2_y_new = (p2[1] * cos(angle) + p2[2] * sin(angle));
        Matrix<T> p2_z_new = (p2[2] * cos(angle) - p2[1] * sin(angle));
        Matrix<T> p3_y_new = (p3[1] * cos(angle) + p3[2] * sin(angle));
        Matrix<T> p3_z_new = (p3[2] * cos(angle) - p3[1] * sin(angle));

        Matrix<T> norm_y_new = norm_y * cos(angle) + norm_z * sin(angle);
        Matrix<T> norm_z_new = norm_z * cos(angle) - norm_y * sin(angle);

        // Update vertices and normals after the rotation.
        p1[1] = p1_y_new; p1[2] = p1_z_new;
        p2[1] = p2_y_new; p2[2] = p2_z_new;
        p3[1] = p3_y_new; p3[2] = p3_z_new;
        norm_y = norm_y_new; norm_z = norm_z_new;

        compute_bounding_box(); // Recompute the bounding box after rotation.
    } else if (axis == 1) {
        // Rotation around the y-axis
        Matrix<T> p1_x_new = p1[0] * cos(angle) - p1[2] * sin(angle);
        Matrix<T> p1_z_new = p1[0] * sin(angle) + p1[2] * cos(angle);
        Matrix<T> p2_x_new = p2[0] * cos(angle) - p2[2] * sin(angle);
        Matrix<T> p2_z_new = p2[0] * sin(angle) + p2[2] * cos(angle);
        Matrix<T> p3_x_new = p3[0] * cos(angle) - p3[2] * sin(angle);
        Matrix<T> p3_z_new = p3[0] * sin(angle) + p3[2] * cos(angle);

        Matrix<T> norm_x_new = (norm_x * cos(angle) - norm_z * sin(angle));
        Matrix<T> norm_z_new = (norm_x * sin(angle) + norm_z * cos(angle));

        // Update vertices and normals after the rotation.
        p1[0] = p1_x_new; p1[2] = p1_z_new;
        p2[0] = p2_x_new; p2[2] = p2_z_new;
        p3[0] = p3_x_new; p3[2] = p3_z_new;
        norm_x = norm_x_new; norm_z = norm_z_new;

        compute_bounding_box(); // Recompute the bounding box after rotation.
    } else if (axis == 2) {
        // Rotation around the z-axis
        Matrix<T> p1_x_new = p1[0] * cos(angle) + p1[1] * sin(angle);
        Matrix<T> p1_y_new = p1[1] * cos(angle) - p1[0] * sin(angle);
        Matrix<T> p2_x_new = p2[0] * cos(angle) + p2[1] * sin(angle);
        Matrix<T> p2_y_new = p2[1] * cos(angle) - p2[0] * sin(angle);
        Matrix<T> p3_x_new = p3[0] * cos(angle) + p3[1] * sin(angle);
        Matrix<T> p3_y_new = p3[1] * cos(angle) - p3[0] * sin(angle);

        Matrix<T> norm_x_new = norm_x * cos(angle) + norm_y * sin(angle);
        Matrix<T> norm_y_new = norm_y * cos(angle) - norm_x * sin(angle);

        // Update vertices and normals after the rotation.
        p1[0] = p1_x_new; p1[1] = p1_y_new;
        p2[0] = p2_x_new; p2[1] = p2_y_new;
        p3[0] = p3_x_new; p3[1] = p3_y_new;
        norm_x = norm_x_new; norm_y = norm_y_new;

        compute_bounding_box(); // Recompute the bounding box after rotation.
    } else {
        std::cout << "The axis must be either 0, 1, or 2!\n";
    }

    // Add the current rotation to the rotation history
    if (record) {
        this->rotation_axis_history.push_back(axis);
        this->rotation_history.push_back(angle);
    }
}

/**
 * @brief Resets the position, scale and orientation of the geometry.
 * 
 * This function cancels all of the shifts, scalings and rotations previously performed.
 * 
 */
template <typename T>
void Geometry<T>::reset() {

    // Stop recording geometry transformations
    this->record = false;

    // Shift the geometry back to its original location
    for(int i = this->shift_history.size() - 1; i  >= 0; i -- ) {
        this->shift(shift_history[i] * (-1.0));
    }

    // Scale the geometry back to its original size
    for(int i = this->scale_history.size() - 1; i  >= 0; i -- ) {
        this->scale(1.0 / scale_history[i]);
    }

    // Rotate the geometry back to its original orientation
    for(int i = this->rotation_history.size() - 1; i  >= 0; i -- ) {
        this->rotate(rotation_axis_history[i], rotation_history[i] * (-1.0));
    }

    // Clear out the history of translations, rotations and scalings
    this->shift_history.clear();
    this->scale_history.clear();
    this->rotation_axis_history.clear();
    this->rotation_history.clear();

    // Resume recording geometry transformations
    this->record = true;

    // Clear GSI and surface properties
    this->GSI_prop.clear();
    this->GSI_prop_names.clear();
    this->surface_prop.clear();
    this->surface_prop_names.clear();
}

/**
 * @brief Computes the area of a triangle at the specified index.
 * 
 * This function computes the area of the triangle at the specified index using the vertex coordinates.
 * 
 * @param index Index of the triangle.
 * @return The computed area of the triangle.
 */
template <typename T>
T Geometry<T>::compute_area(unsigned int index) {
    return 0.5 * fabs(p1[0](index) * (p2[1](index) * p3[2](index) - p3[1](index) * p2[2](index)) + 
                      p1[1](index) * (p3[0](index) * p2[2](index) - p2[0](index) * p3[2](index)) + 
                      p1[2](index) * (p2[0](index) * p3[1](index) - p3[0](index) * p2[1](index)));
}

/**
 * @brief Adds a scalar GSI (Gas-Surface Interaction) property to the geometry.
 * 
 * This function creates a GSI property with a scalar value and adds it to the list of GSI properties.
 * 
 * @param property_name The name of the GSI property.
 * @param value Scalar value of the GSI property.
 */
template <typename T>
void Geometry<T>::add_GSI_property(const char* property_name, T value) {
    Matrix<T> property(1, num_vertices);
    property += value;
    GSI_prop.push_back(property);
    GSI_prop_names.push_back(property_name);
}

/**
 * @brief Adds a matrix GSI (Gas-Surface Interaction) property to the geometry.
 * 
 * This function adds a matrix of values as a GSI property to the list of GSI properties.
 * 
 * @param property_name The name of the GSI property.
 * @param values Matrix of values for the GSI property.
 */
template <typename T>
void Geometry<T>::add_GSI_property(const char* property_name, const Matrix<T>& values) {
    GSI_prop.push_back(values);
    GSI_prop_names.push_back(property_name);
}

/**
 * @brief Adds a scalar surface property to the geometry.
 * 
 * This function creates a surface property with a scalar value and adds it to the list of surface properties.
 * 
 * @param property_name The name of the surface property.
 * @param value Scalar value of the surface property.
 */
template <typename T>
void Geometry<T>::add_surface_property(const char* property_name, T value) {
    Matrix<T> property(1, num_vertices);
    property += value;
    surface_prop.push_back(property);
    surface_prop_names.push_back(property_name);
}

/**
 * @brief Adds a matrix surface property to the geometry.
 * 
 * This function adds a matrix of values as a surface property to the list of surface properties.
 * 
 * @param property_name The name of the surface property.
 * @param values Matrix of values for the surface property.
 */
template <typename T>
void Geometry<T>::add_surface_property(const char* property_name, const Matrix<T>& values) {
    surface_prop.push_back(values);
    surface_prop_names.push_back(property_name);
}

/**
 * @brief Prints the details of the geometry, including vertices, normals, and properties.
 * 
 * This function outputs the details for each triangle in the geometry, such as the vertex coordinates (p1, p2, p3),
 * normal vectors, GSI properties, and surface properties. Each triangle is printed with its corresponding properties.
 */
template <typename T>
void Geometry<T>::print() {
    std::cout << "STL file with " << num_vertices << " triangles:\n\n";

    // Loop through each triangle and print its details
    for (auto i = 0; i < num_vertices; i++) {
        std::cout << "Triangle " << i << ":\n";
        std::cout << "p1: [" << p1[0](i) << " " << p1[1](i) << " " << p1[2](i) << "]\n";
        std::cout << "p2: [" << p2[0](i) << " " << p2[1](i) << " " << p2[2](i) << "]\n";
        std::cout << "p3: [" << p3[0](i) << " " << p3[1](i) << " " << p3[2](i) << "]\n";

        // Print the normal vector of the triangle
        std::cout << "normal: [" << norm_x(i) << " " << norm_y(i) << " " << norm_z(i) << "]\n";

        // Print the names of GSI properties
        std::cout << "GSI properties names: ";
        for (auto j = 0; j < GSI_prop_names.size(); j++) {
            std::cout << GSI_prop_names[j] << " ";
        }

        // Print the values of GSI properties for the current triangle
        std::cout << "\nGSI properties values: ";
        for (auto j = 0; j < GSI_prop.size(); j++) {
            std::cout << GSI_prop[j](i) << " ";
        }

        // Print the names of surface properties
        std::cout << "\nSurface properties names: ";
        for (auto j = 0; j < surface_prop_names.size(); j++) {
            std::cout << surface_prop_names[j] << " ";
        }

        // Print the values of surface properties for the current triangle
        std::cout << "\nSurface properties values: ";
        for (auto j = 0; j < surface_prop.size(); j++) {
            std::cout << surface_prop[j](i) << " ";
        }
        std::cout << "\n\n";
    }
}

/**
 * @brief Returns the bounding box of the geometry.
 * 
 * This function provides the bounding box matrix, which contains the minimum and maximum coordinates in x, y, and z directions.
 * 
 * @return Matrix<T> The bounding box of the geometry.
 */
template <typename T>
Matrix<T> Geometry<T>::get_bounding_box() {
    return bounding_box;
}
