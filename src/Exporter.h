/**
 * @brief Exports geometry and aerodynamic coefficients to a CSV file.
 * 
 * This function writes the geometric data of the satellite and its associated aerodynamic coefficients (pressure and shear) 
 * into a CSV file. The geometry includes vertex positions and normal vectors.
 *
 * @tparam T Type of the values (e.g., float, double).
 * @param filename The name of the output CSV file.
 * @param satellite The geometry object representing the satellite or object.
 * @param norm_coeff Array of pressure coefficients.
 * @param shear_x_coeff Array of X-component of the shear force coefficients.
 * @param shear_y_coeff Array of Y-component of the shear force coefficients.
 * @param shear_z_coeff Array of Z-component of the shear force coefficients.
 */
template <typename T>
void export_geometry(std::string filename, Geometry<T> satellite, T* norm_coeff, T* shear_x_coeff, T* shear_y_coeff, T* shear_z_coeff) {

    // Open file for writing
    std::ofstream write(filename);

    // Matrices to hold vertex positions and normal vectors
    Matrix<T> p1_x(1, satellite.get_num_vertices()), p1_y(1, satellite.get_num_vertices()), p1_z(1, satellite.get_num_vertices());
    Matrix<T> p2_x(1, satellite.get_num_vertices()), p2_y(1, satellite.get_num_vertices()), p2_z(1, satellite.get_num_vertices());
    Matrix<T> p3_x(1, satellite.get_num_vertices()), p3_y(1, satellite.get_num_vertices()), p3_z(1, satellite.get_num_vertices());
    Matrix<T> norm_x(1, satellite.get_num_vertices()), norm_y(1, satellite.get_num_vertices()), norm_z(1, satellite.get_num_vertices());

    // Retrieve vertex coordinates and normal vectors from the geometry object
    p1_x = satellite.get_p1_x(); p1_y = satellite.get_p1_y(); p1_z = satellite.get_p1_z();
    p2_x = satellite.get_p2_x(); p2_y = satellite.get_p2_y(); p2_z = satellite.get_p2_z();
    p3_x = satellite.get_p3_x(); p3_y = satellite.get_p3_y(); p3_z = satellite.get_p3_z();
    norm_x = satellite.get_norm_x(); norm_y = satellite.get_norm_y(); norm_z = satellite.get_norm_z();

    // Write the header line in the CSV file
    write << "Index; Point_1 X-value; Point_1 Y-value; Point_1 Z-value; "
          << "Point_2 X-value; Point_2 Y-value; Point_2 Z-value; "
          << "Point_3 X-value; Point_3 Y-value; Point_3 Z-value; "
          << "Normal X-value; Normal Y-value; Normal Z-value; "
          << "Pressure Coefficient; X-Shear Coefficient; Y-Shear Coefficient; Z-Shear Coefficient\n";

    // Write vertex positions, normals, and aerodynamic coefficients to the file
    for (auto i = 0; i < satellite.get_num_vertices(); i++) {
        write << i << "; " << p1_x(i) << "; " << p1_y(i) << "; " << p1_z(i) << "; "
              << p2_x(i) << "; " << p2_y(i) << "; " << p2_z(i) << "; "
              << p3_x(i) << "; " << p3_y(i) << "; " << p3_z(i) << "; "
              << norm_x(i) << "; " << norm_y(i) << "; " << norm_z(i) << "; "
              << norm_coeff[i] << "; " << shear_x_coeff[i] << "; " << shear_y_coeff[i] << "; " << shear_z_coeff[i] << "\n";
    }

    // Close the file after writing
    write.close();
}

/**
 * @brief Exports aerodynamic coefficients for a batch of parameters to a CSV file.
 * 
 * This function writes aerodynamic force and moment coefficients for a batch of parameters into a CSV file.
 *
 * @tparam T Type of the values (e.g., float, double).
 * @param foldername The name of the output folder.
 * @param param_1_name Name of the first parameter.
 * @param param_2_name Name of the second parameter.
 * @param param_3_name Name of the third parameter.
 * @param param_4_name Name of the fourth parameter.
 * @param param_5_name Name of the fifth parameter.
 * @param param_1 Vector of values for the first parameter.
 * @param param_2 Vector of values for the second parameter.
 * @param param_3 Vector of values for the third parameter.
 * @param param_4 Vector of values for the fourth parameter.
 * @param param_5 Vector of values for the fifth parameter.
 * @param force_coeff Vector of force coefficient matrices.
 * @param moment_coeff Vector of moment coefficient matrices.
 */
template <typename T>
void export_aerodynamic_coefficients_batch(std::string foldername, std::string param_1_name, std::string param_2_name, std::string param_3_name, std::string param_4_name, std::string param_5_name,
                                           std::vector<T> param_1, std::vector<T> param_2, std::vector<T> param_3, std::vector<T> param_4, std::vector<T> param_5, 
                                           std::vector<Matrix<T>> force_coeff, std::vector<Matrix<T>> moment_coeff) {
    
    // Open file for writing aerodynamic coefficients
    std::ofstream write(foldername + "/aerodynamic_coefficients.csv");
    write << "sep=;\n";  // Use semicolon as separator
    write << param_1_name << "; " << param_2_name << "; " << param_3_name << "; " << param_4_name << "; " << param_5_name 
          << "; Force Coefficient X-value; Force Coefficient Y-value; Force Coefficient Z-value; "
          << "Moment Coefficient X-value; Moment Coefficient Y-value; Moment Coefficient Z-value \n";

    // Write aerodynamic coefficients for each combination of parameters
    for (auto i = 0; i < param_1.size(); i++) {
        for (auto j = 0; j < param_2.size(); j++) {
            for (auto k = 0; k < param_3.size(); k++) {
                for (auto l = 0; l < param_4.size(); l++) {
                    for (auto m = 0; m < param_5.size(); m++) {
                        long index = i * param_2.size() * param_3.size() * param_4.size() * param_5.size() +
                                     j * param_3.size() * param_4.size() * param_5.size() +
                                     k * param_4.size() * param_5.size() +
                                     l * param_5.size() + m;

                        write << param_1[i] << " ; " << param_2[j] << " ; " << param_3[k] << " ; " << param_4[l] << " ; " << param_5[m] << " ; "
                              << force_coeff  << " ; " << force_coeff  << " ; " << force_coeff  << " ; "
                              << moment_coeff  << " ; " << moment_coeff  << " ; " << moment_coeff  << "\n";
                    }
                }
            }
        }
    }

    // Close the file after writing
    write.close();
}

/**
 * @brief Appends a single line of aerodynamic coefficients to an existing CSV file.
 * 
 * This function appends a line of aerodynamic force and moment coefficients for a given set of parameters to a CSV file.
 * 
 * @tparam T Type of the values (e.g., float, double).
 * @param foldername The name of the folder containing the output file.
 * @param param_1_name Name of the first parameter.
 * @param param_2_name Name of the second parameter.
 * @param param_3_name Name of the third parameter.
 * @param param_4_name Name of the fourth parameter.
 * @param param_5_name Name of the fifth parameter.
 * @param param_1 Value of the first parameter.
 * @param param_2 Value of the second parameter.
 * @param param_3 Value of the third parameter.
 * @param param_4 Value of the fourth parameter.
 * @param param_5 Value of the fifth parameter.
 * @param force_coeff Matrix of force coefficients.
 * @param moment_coeff Matrix of moment coefficients.
 * @param first_line If true, the header line is written.
 */
template <typename T>
void export_aerodynamic_coefficients_batch_line(std::string foldername, std::string param_1_name, std::string param_2_name, std::string param_3_name, std::string param_4_name, std::string param_5_name,
                                                T param_1, T param_2, T param_3, T param_4, T param_5, Matrix<T> force_coeff, Matrix<T> moment_coeff, bool first_line, int sim_index) {
    
    // Open file for appending aerodynamic coefficients
    std::stringstream ss;
    ss << sim_index;
    std::ofstream write(foldername + "/aerodynamic_coefficients_" + ss.str()  + ".csv", std::ios_base::app);

    // Write header line if it's the first line
    if (first_line) {
        write << "sep=;\n";
        write << param_1_name << "; " << param_2_name << "; " << param_3_name << "; " << param_4_name << "; " << param_5_name
              << "; Force Coefficient X-value; Force Coefficient Y-value; Force Coefficient Z-value; "
              << "Moment Coefficient X-value; Moment Coefficient Y-value; Moment Coefficient Z-value \n";
    }

    // Append aerodynamic coefficients for the given parameters
    write << param_1 << " ; " << param_2 << " ; " << param_3 << " ; " << param_4 << " ; " << param_5 << " ; "
          << force_coeff(0) << " ; " << force_coeff(1) << " ; " << force_coeff(2) << " ; "
          << moment_coeff(0) << " ; " << moment_coeff(1) << " ; " << moment_coeff(2) << "\n";

    // Close the file after writing
    write.close();
}

/**
 * @brief Exports a single set of aerodynamic coefficients to a CSV file.
 * 
 * This function writes a single set of aerodynamic force and moment coefficients into a CSV file.
 *
 * @tparam T Type of the values (e.g., float, double).
 * @param foldername The name of the output folder.
 * @param force_coeff Matrix of force coefficients.
 * @param moment_coeff Matrix of moment coefficients.
 */
template <typename T>
void export_aerodynamic_coefficients_single(std::string foldername, Matrix<T> force_coeff, Matrix<T> moment_coeff) {
    
    // Open file for writing aerodynamic coefficients
    std::ofstream write(foldername + "/aerodynamic_coefficients.csv");
    write << "sep=;\n";  // Use semicolon as separator
    write << "Force Coefficient X-value; Force Coefficient Y-value; Force Coefficient Z-value; "
          << "Moment Coefficient X-value; Moment Coefficient Y-value; Moment Coefficient Z-value\n";

    // Write force and moment coefficients to the file
    write << force_coeff(0) << " ; " << force_coeff(1) << " ; " << force_coeff(2) << " ; "
          << moment_coeff(0) << " ; " << moment_coeff(1) << " ; " << moment_coeff(2) << "\n";

    // Close the file after writing
    write.close();
}

/**
 * @brief Exports the results of a simulation and geometry data to a directory.
 * 
 * This function creates a directory for storing simulation results and exports aerodynamic coefficients and geometry data.
 *
 * @tparam T Type of the values (e.g., float, double).
 * @param foldername The name of the output folder.
 * @param sim The simulation object containing geometry and other relevant data.
 * @param force_coeff Matrix of force coefficients.
 * @param moment_coeff Matrix of moment coefficients.
 * @param norm_coeff Array of pressure coefficients.
 * @param shear_x_coeff Array of X-component shear force coefficients.
 * @param shear_y_coeff Array of Y-component shear force coefficients.
 * @param shear_z_coeff Array of Z-component shear force coefficients.
 */
template <typename T>
void export_single(std::string foldername, Simulation<T> sim, Matrix<T> force_coeff, Matrix<T> moment_coeff, T* norm_coeff, T* shear_x_coeff, T* shear_y_coeff, T* shear_z_coeff) {

    // Create directory for single simulation results
    create_directory_if_missing(foldername + "/Single_Results");

    // Export aerodynamic coefficients to a CSV file
    export_aerodynamic_coefficients_single(foldername + "/Single_Results", force_coeff, moment_coeff);

    // Export geometry data to a CSV file
    std::string geometry_name = "geometry.csv";
    export_geometry<T>(foldername + "/Single_Results/" + geometry_name, sim.get_geometry(), norm_coeff, shear_x_coeff, shear_y_coeff, shear_z_coeff);
}
