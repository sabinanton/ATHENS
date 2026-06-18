/**
 * @brief Skips a specified number of lines in a file.
 * 
 * This function reads and discards a specified number of lines from the input file stream.
 * 
 * @param file The input file stream from which lines are to be skipped.
 * @param skip_lines The number of lines to skip.
 */
void skip_lines(std::ifstream& file, int skip_lines) {
    std::string temp;
    for (auto i = 0; i < skip_lines; i++) {
        // Read and discard lines
        if (!std::getline(file, temp)) {
            std::cerr << "Failed to read or fewer lines in file." << std::endl;
        }
    }
}

/**
 * @brief Imports a batch analysis simulation configuration from the input folder.
 * 
 * This function reads input parameters and constructs a simulation for batch analysis. 
 * It configures the simulation's domain, satellite geometry, surface, and gas properties.
 * 
 * @tparam T Type of the simulation parameters (e.g., float, double).
 * @param foldername The folder containing input files.
 * @param inputs Vector of input values as strings.
 * @param satellite The satellite geometry object. 
 * @param index1 Index for alpha parameter.
 * @param index2 Index for s parameter.
 * @param index3 Index for angle of attack (AoA).
 * @param index4 Index for sideslip (SS).
 * @param index5 Index for sigma_t parameter.
 * @param num_procs Number of processors to use for parallel processing.
 * @return Simulation<T> The configured simulation object.
 */
template <typename T>
Simulation<T> import_batch_analysis(std::string foldername, std::vector<std::string> inputs, Geometry<T>& satellite, int index1, int index2, int index3, int index4, int index5, int num_procs) {
    // Convert input strings to corresponding simulation parameters
    T scale = static_cast<T>(std::atof(inputs[3].c_str()));
    T offset_x = static_cast<T>(std::atof(inputs[4].c_str()));
    T offset_y = static_cast<T>(std::atof(inputs[5].c_str()));
    T offset_z = static_cast<T>(std::atof(inputs[6].c_str()));
    T rotation_x = static_cast<T>(std::atof(inputs[7].c_str()));
    T rotation_y = static_cast<T>(std::atof(inputs[8].c_str()));
    T rotation_z = static_cast<T>(std::atof(inputs[9].c_str()));
    T A_ref = static_cast<T>(std::atof(inputs[10].c_str()));
    T c_ref = static_cast<T>(std::atof(inputs[11].c_str()));
    T surf_temp = static_cast<T>(std::atof(inputs[12].c_str()));

    T gas_molar_mass = static_cast<T>(std::atof(inputs[14].c_str()));
    T gas_temperature = static_cast<T>(std::atof(inputs[15].c_str()));
    T gas_speed = static_cast<T>(std::atof(inputs[16].c_str()));

    std::string kernel_name = inputs[18];
    void (*kernel)(T*, T*, T*, const Gas<T>&, Surface<T>&);

    // Select appropriate kernel function based on kernel name
    if (kernel_name == "DRIA") kernel = sample_DRIA;
    else if (kernel_name == "CLL") kernel = sample_CLL;
    else kernel = sample_rough_CLL;

    // Load Gas-Surface Interaction (GSI) parameters
    T GSI_parameters[3];
    GSI_parameters[0] = static_cast<T>(std::atof(inputs[19].c_str()));
    GSI_parameters[1] = static_cast<T>(std::atof(inputs[20].c_str()));
    T sigma_T = static_cast<T>(std::atof(inputs[21].c_str()));

    // Simulation parameters
    long num_rays = static_cast<long>(std::atof(inputs[23].c_str()));
    int ray_children = static_cast<int>(std::atof(inputs[24].c_str()));
    int max_bounces = static_cast<int>(std::atof(inputs[25].c_str()));

    // Define parameter ranges for simulation
    T param_1_lower = static_cast<T>(std::atof(inputs[30].c_str()));
    T param_1_upper = static_cast<T>(std::atof(inputs[31].c_str()));
    T param_1_step = static_cast<T>(std::atof(inputs[32].c_str()));
    T param_2_lower = static_cast<T>(std::atof(inputs[33].c_str()));
    T param_2_upper = static_cast<T>(std::atof(inputs[34].c_str()));
    T param_2_step = static_cast<T>(std::atof(inputs[35].c_str()));
    T param_3_lower = static_cast<T>(std::atof(inputs[36].c_str()));
    T param_3_upper = static_cast<T>(std::atof(inputs[37].c_str()));
    T param_3_step = static_cast<T>(std::atof(inputs[38].c_str()));
    T param_4_lower = static_cast<T>(std::atof(inputs[39].c_str()));
    T param_4_upper = static_cast<T>(std::atof(inputs[40].c_str()));
    T param_4_step = static_cast<T>(std::atof(inputs[41].c_str()));
    T param_5_lower = static_cast<T>(std::atof(inputs[42].c_str()));
    T param_5_upper = static_cast<T>(std::atof(inputs[43].c_str()));
    T param_5_step = static_cast<T>(std::atof(inputs[44].c_str()));

    // Pre-calculate the number of simulations based on parameter ranges
    int num_sims_1 = (int)((param_1_upper - param_1_lower) / param_1_step) + 1;
    int num_sims_2 = (int)((param_2_upper - param_2_lower) / param_2_step) + 1;
    int num_sims_3 = (int)((param_3_upper - param_3_lower) / param_3_step) + 1;
    int num_sims_4 = (int)((param_4_upper - param_4_lower) / param_4_step) + 1;
    int num_sims_5 = (int)((param_5_upper - param_5_lower) / param_5_step) + 1;

    // Initialize vectors to hold parameter values for simulations
    std::vector<T> AoA, SS, s, alpha_n, sigma_t;
    for (auto i = 0; i < num_sims_1; i++) AoA.push_back((param_1_lower + i * param_1_step) * M_PI / 180.0);
    for (auto i = 0; i < num_sims_2; i++) SS.push_back((param_2_lower + i * param_2_step) * M_PI / 180.0);
    for (auto i = 0; i < num_sims_3; i++) s.push_back((param_3_lower + i * param_3_step));
    for (auto i = 0; i < num_sims_4; i++) alpha_n.push_back((param_4_lower + i * param_4_step));
    for (auto i = 0; i < num_sims_5; i++) sigma_t.push_back((param_5_lower + i * param_5_step));

    // Load satellite geometry and surface properties
    // Geometry<T> satellite((foldername + "/surface_files/" + inputs[1]).c_str()); - deprecated
    Surface<T> surface((foldername + "/surface_files/" + inputs[2]).c_str());

    // Apply transformations to the satellite geometry
    T offset[4] = {offset_x, offset_y, offset_z};
    satellite.reset();
    satellite.scale(scale);
    satellite.rotate(0, rotation_x * M_PI / 180.0);
    satellite.rotate(1, rotation_y * M_PI / 180.0);
    satellite.rotate(2, rotation_z * M_PI / 180.0);
    satellite.shift(satellite.get_center() * (-1.0));
    satellite.shift(offset);
    satellite.set_A_ref(A_ref);
    surface.set_temperature(surf_temp);

    // Assign Gas-Surface Interaction (GSI) properties to the satellite
    if (kernel_name == "DRIA") {
        satellite.add_GSI_property("alpha", alpha_n[index1]);
    } else {
        satellite.add_GSI_property("alpha_n", alpha_n[index2]);
        satellite.add_GSI_property("sigma_t", sigma_t[index1]);
    }

    // Rotate the satellite based on AoA and SS
    satellite.rotate(1, AoA[index4]);
    satellite.rotate(2, SS[index5]);


    // Define gas properties and domain for the simulation
    Gas<T> atmosphere(1.0, gas_temperature, gas_molar_mass, s[index3] * std::sqrt(2.0 * R_gas * gas_temperature / gas_molar_mass));
    Matrix<double> bounding_box = satellite.get_bounding_box();
    double margin = 0.1;
    satellite.shift(satellite.get_center() * (-1.0));
    Domain<T> domain((bounding_box(1, 0) - bounding_box(0, 0)) + margin, 
                     (bounding_box(3, 0) - bounding_box(2, 0)) + margin,
                     (bounding_box(5, 0) - bounding_box(4, 0)) + margin, atmosphere, s[index3] * std::sqrt(2.0 * R_gas * gas_temperature / gas_molar_mass));

    // Create the simulation object
    Simulation<T> sim(domain, satellite, surface, kernel, (num_rays / num_procs) * num_procs, ray_children, max_bounces);

    return sim;
}

/**
 * @brief Imports a single simulation configuration from the input folder.
 * 
 * This function reads input parameters from a file and constructs a simulation configuration for a single run.
 * It configures the satellite geometry, surface, gas properties, and domain.
 * 
 * @tparam T Type of the simulation parameters (e.g., float, double).
 * @param foldername The folder containing input files.
 * @param num_procs Number of processors to use for parallel processing.
 * @return Simulation<T> The configured simulation object.
 */
template <typename T>
Simulation<T> import_single(std::string foldername, int num_procs) {

    // Open the input file containing simulation parameters
    std::string filename = foldername + "/Inputs.txt";
    std::ifstream file(filename);
    std::vector<std::string> inputs;
    std::string temp;
    std::string param_name, equal_sign, value;

    // Read the input file and store values
    for (int i = 0; i < 60; i++) {
        std::getline(file, temp);
        std::stringstream ss(temp);
        if (ss >> param_name >> equal_sign >> value) {
            inputs.push_back(value);
        }
    }
    // Close the file
    file.close();

    // Load satellite geometry and surface properties
    Geometry<T> satellite((foldername + "/surface_files/" + inputs[1]).c_str());
    Surface<T> surface((foldername + "/surface_files/" + inputs[2]).c_str());

    // Convert input strings to corresponding simulation parameters
    T scale = static_cast<T>(std::atof(inputs[3].c_str()));
    T offset_x = static_cast<T>(std::atof(inputs[4].c_str()));
    T offset_y = static_cast<T>(std::atof(inputs[5].c_str()));
    T offset_z = static_cast<T>(std::atof(inputs[6].c_str()));
    T rotation_x = static_cast<T>(std::atof(inputs[7].c_str()));
    T rotation_y = static_cast<T>(std::atof(inputs[8].c_str()));
    T rotation_z = static_cast<T>(std::atof(inputs[9].c_str()));
    T A_ref = static_cast<T>(std::atof(inputs[10].c_str()));
    T c_ref = static_cast<T>(std::atof(inputs[11].c_str()));
    T surf_temp = static_cast<T>(std::atof(inputs[12].c_str()));

    T gas_molar_mass = static_cast<T>(std::atof(inputs[14].c_str()));
    T gas_temperature = static_cast<T>(std::atof(inputs[15].c_str()));
    T gas_speed = static_cast<T>(std::atof(inputs[16].c_str()));

    // Determine the gas-surface interaction (GSI) kernel
    std::string kernel_name = inputs[18];
    void (*kernel)(T*, T*, T*, const Gas<T>&, Surface<T>&);

    // Choose the kernel function based on the kernel name
    if (kernel_name == "DRIA") kernel = sample_DRIA;
    else if (kernel_name == "CLL") kernel = sample_CLL;
    else kernel = sample_rough_CLL;

    // Load GSI parameters
    T GSI_parameters[3];
    GSI_parameters[0] = static_cast<T>(std::atof(inputs[19].c_str()));
    GSI_parameters[1] = static_cast<T>(std::atof(inputs[20].c_str()));
    T sigma_T = static_cast<T>(std::atof(inputs[21].c_str()));

    // Simulation parameters
    long num_rays = static_cast<long>(std::atof(inputs[23].c_str()));
    int ray_children = static_cast<int>(std::atof(inputs[24].c_str()));
    int max_bounces = static_cast<int>(std::atof(inputs[25].c_str()));
    T AoA = static_cast<T>(std::atof(inputs[26].c_str()));
    T SS = static_cast<T>(std::atof(inputs[27].c_str()));

    // Apply transformations to the satellite geometry
    T offset[4] = {offset_x, offset_y, offset_z};
    satellite.scale(scale);
    satellite.rotate(0, rotation_x * M_PI / 180.0);
    satellite.rotate(1, rotation_y * M_PI / 180.0);
    satellite.rotate(2, rotation_z * M_PI / 180.0);
    satellite.shift(satellite.get_center() * (-1.0));
    satellite.shift(offset);
    satellite.set_A_ref(A_ref);
    surface.set_temperature(surf_temp);

    // Apply rotations for angle of attack (AoA) and sideslip (SS)
    satellite.rotate(2, SS * M_PI / 180.0);
    satellite.rotate(1, AoA * M_PI / 180.0);

    // Assign GSI properties based on the kernel
    if (kernel_name == "DRIA") {
        satellite.add_GSI_property("alpha", GSI_parameters[0]);
    } else {
        satellite.add_GSI_property("alpha_n", GSI_parameters[0]);
        satellite.add_GSI_property("sigma_t", GSI_parameters[1]);
    }

    // Define gas properties and domain for the simulation
    Gas<T> atmosphere(1.0, gas_temperature, gas_molar_mass, gas_speed);
    Matrix<double> bounding_box = satellite.get_bounding_box();
    double margin = 0.1;
    satellite.shift(satellite.get_center() * (-1.0));
    Domain<T> domain((bounding_box(1, 0) - bounding_box(0, 0)) + margin,
                     (bounding_box(3, 0) - bounding_box(2, 0)) + margin,
                     (bounding_box(5, 0) - bounding_box(4, 0)) + margin, atmosphere, gas_speed);

    // Create and return the simulation object
    Simulation<T> sim(domain, satellite, surface, kernel, (num_rays / num_procs) * num_procs, ray_children, max_bounces);
    return sim;
}

/**
 * @brief Imports a batch simulation configuration from the input folder.
 * 
 * This function reads input parameters and constructs a batch simulation. 
 * It calls the import_batch_analysis function to handle the specific analysis setup.
 * 
 * @tparam T Type of the simulation parameters (e.g., float, double).
 * @param foldername The folder containing input files.
 * @param index1 Index for alpha parameter.
 * @param index2 Index for s parameter.
 * @param index3 Index for angle of attack (AoA).
 * @param index4 Index for sideslip (SS).
 * @param index5 Index for sigma_t parameter.
 * @param num_procs Number of processors to use for parallel processing.
 * @return Simulation<T> The configured batch simulation object.
 */
template <typename T>
Simulation<T> import_batch(std::string foldername, Geometry<T>& satellite, int index1, int index2, int index3, int index4, int index5, int num_procs) {

    // Open the input file containing simulation parameters
    std::string filename = foldername + "/Inputs.txt";
    std::ifstream file(filename);
    std::vector<std::string> inputs;
    std::string temp;
    std::string param_name, equal_sign, value;

    // Read the input file and store values
    for (int i = 0; i < 80; i++) {
        std::getline(file, temp);
        std::stringstream ss(temp);
        if (ss >> param_name >> equal_sign >> value) {
            inputs.push_back(value);
        }
    }

    // Close the file
    file.close();

    // Call the import_batch_analysis function to handle batch-specific setup
    Simulation<T> sim = import_batch_analysis<T>(foldername, inputs, satellite, index1, index2, index3, index4, index5, num_procs);

    std::string batch_type = inputs[29];
    std::string kernel_name = inputs[18];

    return sim;
}

/**
 * @brief Imports the satellite geometry.
 * 
 * This function reads geometry from an STL file.
 * 
 * @tparam T Type of the simulation parameters (e.g., float, double).
 * @param foldername The folder containing the geometry name.
 * @return The geometry object.
 */
template <typename T>
Geometry<T> import_geometry(std::string foldername) {

    // Open the input file containing simulation parameters
    std::string filename = foldername + "/Inputs.txt";
    std::ifstream file(filename);
    std::vector<std::string> inputs;
    std::string temp;
    std::string param_name, equal_sign, value;

    // Read the input file and store values
    for (int i = 0; i < 80; i++) {
        std::getline(file, temp);
        std::stringstream ss(temp);
        if (ss >> param_name >> equal_sign >> value) {
            inputs.push_back(value);
        }
    }

    // Close the file
    file.close();

    // Import the satellite geometry
    Geometry<T> satellite((foldername + "/surface_files/" + inputs[1]).c_str());

    return satellite;
}

/**
 * @brief Imports simulation parameters from the input folder for batch analysis.
 * 
 * This function reads input parameters and returns them as a tuple, including the parameter names and their ranges.
 * 
 * @tparam T Type of the simulation parameters (e.g., float, double).
 * @param foldername The folder containing input files.
 * @return A tuple containing parameter names and corresponding vectors of parameter values.
 */
template <typename T>
std::tuple<std::string, std::string, std::string, std::string, std::string, std::vector<T>, std::vector<T>, std::vector<T>, std::vector<T>, std::vector<T>> import_parameters(std::string foldername) {
    // Open the input file containing simulation parameters
    std::string filename = foldername + "/Inputs.txt";
    std::ifstream file(filename);
    std::vector<std::string> inputs;
    std::string temp;
    std::string param_name, equal_sign, value;

    // Read the input file and store values
    for (int i = 0; i < 80; i++) {
        std::getline(file, temp);
        std::stringstream ss(temp);
        if (ss >> param_name >> equal_sign >> value) {
            inputs.push_back(value);
        }
    }

    // Extract parameter ranges from the input values
    T param_1_lower = static_cast<T>(std::atof(inputs[30].c_str()));
    T param_1_upper = static_cast<T>(std::atof(inputs[31].c_str()));
    T param_1_step = static_cast<T>(std::atof(inputs[32].c_str()));
    T param_2_lower = static_cast<T>(std::atof(inputs[33].c_str()));
    T param_2_upper = static_cast<T>(std::atof(inputs[34].c_str()));
    T param_2_step = static_cast<T>(std::atof(inputs[35].c_str()));
    T param_3_lower = static_cast<T>(std::atof(inputs[36].c_str()));
    T param_3_upper = static_cast<T>(std::atof(inputs[37].c_str()));
    T param_3_step = static_cast<T>(std::atof(inputs[38].c_str()));
    T param_4_lower = static_cast<T>(std::atof(inputs[39].c_str()));
    T param_4_upper = static_cast<T>(std::atof(inputs[40].c_str()));
    T param_4_step = static_cast<T>(std::atof(inputs[41].c_str()));
    T param_5_lower = static_cast<T>(std::atof(inputs[42].c_str()));
    T param_5_upper = static_cast<T>(std::atof(inputs[43].c_str()));
    T param_5_step = static_cast<T>(std::atof(inputs[44].c_str()));

    // Pre-calculate the number of simulations based on parameter ranges
    int num_sims_1 = (int)((param_1_upper - param_1_lower) / param_1_step) + 1;
    int num_sims_2 = (int)((param_2_upper - param_2_lower) / param_2_step) + 1;
    int num_sims_3 = (int)((param_3_upper - param_3_lower) / param_3_step) + 1;
    int num_sims_4 = (int)((param_4_upper - param_4_lower) / param_4_step) + 1;
    int num_sims_5 = (int)((param_5_upper - param_5_lower) / param_5_step) + 1;

    // Initialize vectors to hold parameter values
    std::vector<T> AoA, SS, s, alpha_n, sigma_t;

    for (auto i = 0; i < num_sims_1; i++) AoA.push_back((param_1_lower + i * param_1_step) * M_PI / 180.0);
    for (auto i = 0; i < num_sims_2; i++) SS.push_back((param_2_lower + i * param_2_step) * M_PI / 180.0);
    for (auto i = 0; i < num_sims_3; i++) s.push_back((param_3_lower + i * param_3_step));
    for (auto i = 0; i < num_sims_4; i++) alpha_n.push_back((param_4_lower + i * param_4_step));
    for (auto i = 0; i < num_sims_5; i++) sigma_t.push_back((param_5_lower + i * param_5_step));

    std::string batch_type = inputs[29];
    std::string kernel_name = inputs[18];

    // Return the parameter names and corresponding vectors
    if (kernel_name == "DRIA") {
        std::string param_name_1 = "alpha";
        std::string param_name_2 = "NA";
        std::string param_name_3 = "s";
        std::string param_name_4 = "AoA";
        std::string param_name_5 = "SS";
        return std::make_tuple(param_name_1, param_name_2, param_name_3, param_name_4, param_name_5, alpha_n, sigma_t, s, AoA, SS);
    } else {
        std::string param_name_1 = "sigma_t";
        std::string param_name_2 = "alpha_n";
        std::string param_name_3 = "s";
        std::string param_name_4 = "AoA";
        std::string param_name_5 = "SS";
        return std::make_tuple(param_name_1, param_name_2, param_name_3, param_name_4, param_name_5, sigma_t, alpha_n, s, AoA, SS);
    }
}


template <typename T>
long find_index(std::string foldername, int sim_index) {

    long index = 0, num_lines = 0;

    std::stringstream ss;

    ss << sim_index;

     // Open the input file containing simulation parameters
    std::string filename = foldername + "/Batch_Results/aerodynamic_coefficients_" + ss.str() + ".csv";
    std::ifstream file(filename);
    std::string linestring;

    while (std::getline(file, linestring))
        ++num_lines;
    
    file.close();
    if(num_lines - 2 < 0) return 0;

    return num_lines - 2;
}
