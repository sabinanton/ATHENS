/**
 * @brief A template class that defines a 3D domain with an atmosphere and gas dynamics.
 * 
 * @tparam T The data type for the domain dimensions and gas speed.
 */
template <typename T>

class Domain {

    private:
        T domain_size_x, domain_size_y, domain_size_z; /**< Dimensions of the domain in x, y, and z directions. */
        Gas<T> atmosphere; /**< Atmosphere model for the domain. */
        T gas_speed_magnitude; /**< Magnitude of the gas speed. */

    public:

        /**
         * @brief Default constructor for the Domain class.
         * Initializes the domain size to zero.
         */

        Domain() {
            domain_size_x = 0.0;
            domain_size_y = 0.0;
            domain_size_z = 0.0;
        }
        
        /**
         * @brief Parameterized constructor for the Domain class.
         * 
         * @param size_x The size of the domain in the x-direction.
         * @param size_y The size of the domain in the y-direction.
         * @param size_z The size of the domain in the z-direction.
         * @param atmosphere The atmospheric gas in the domain.
         * @param gas_speed The magnitude of the gas speed.
         */
        Domain(T size_x, T size_y, T size_z, Gas<T> atmosphere, T gas_speed) 
        : domain_size_x(size_x), domain_size_y(size_y), domain_size_z(size_z), atmosphere(atmosphere), gas_speed_magnitude(gas_speed) {}
        
        /**
         * @brief Copy constructor for the Domain class.
         * 
         * @param other The Domain object to copy from.
         */
        Domain(const Domain<T>& other)
        : domain_size_x(other.domain_size_x), domain_size_y(other.domain_size_y), domain_size_z(other.domain_size_z), atmosphere(other.atmosphere), gas_speed_magnitude(other.gas_speed_magnitude) {}
        
        /**
         * @brief Move constructor for the Domain class.
         * 
         * @param other The Domain object to move from.
         */
        Domain(Domain<T>&& other)
        : domain_size_x(other.domain_size_x), domain_size_y(other.domain_size_y), domain_size_z(other.domain_size_z), atmosphere(other.atmosphere), gas_speed_magnitude(other.gas_speed_magnitude) {
            other.domain_size_x = 0.0;
            other.domain_size_y = 0.0;
            other.domain_size_z = 0.0;
            other.gas_speed_magnitude = 0.0;
        }

        /**
         * @brief Destructor for the Domain class.
         */
        ~Domain() {
            domain_size_x = 0.0;
            domain_size_y = 0.0;
            domain_size_z = 0.0;
        }

        /**
         * @brief Copy assignment operator for the Domain class.
         * 
         * @param other The Domain object to copy from.
         * @return Reference to the assigned Domain object.
         */
        Domain& operator=(const Domain& other);

        /**
         * @brief Move assignment operator for the Domain class.
         * 
         * @param other The Domain object to move from.
         * @return Reference to the assigned Domain object.
         */
        Domain& operator=(Domain&& other);

        /**
         * @brief Get the size of the domain in the x-direction.
         * 
         * @return The domain size in the x-direction.
         */
        T get_domain_size_x() {return domain_size_x;}

        /**
         * @brief Get the size of the domain in the y-direction.
         * 
         * @return The domain size in the y-direction.
         */
        T get_domain_size_y() {return domain_size_y;}

        /**
         * @brief Get the size of the domain in the z-direction.
         * 
         * @return The domain size in the z-direction.
         */
        T get_domain_size_z() {return domain_size_z;}

        /**
         * @brief Get the atmosphere of the domain.
         * 
         * @return The atmosphere model of the domain.
         */
        Gas<T> get_atmosphere() {return atmosphere;}

        /**
         * @brief Get the magnitude of the gas speed.
         * 
         * @return The magnitude of the gas speed.
         */
        T get_gas_speed_magnitude() {return gas_speed_magnitude;}

        /**
         * @brief Generate rays within the domain.
         * 
         * @param num_rays The number of rays to generate.
         * @param num_children The number of child rays.
         * @return A vector containing the generated rays.
         */
        std::vector<Ray<T>> generate_rays(unsigned long num_rays, int num_children);

        /**
         * @brief Generate rays in parallel within a specified range.
         * 
         * @param num_rays The total number of rays.
         * @param start_index The starting index for the rays to generate.
         * @param end_index The ending index for the rays to generate.
         * @param num_children The number of child rays.
         * @return A vector containing the generated rays within the specified range.
         */
        std::vector<Ray<T>> generate_rays_parallel(unsigned long num_rays, int start_index, int end_index, int num_children);
};

template <typename T>
Domain<T>& Domain<T>::operator=(const Domain<T>& other) {
    domain_size_x = other.domain_size_x;
    domain_size_y = other.domain_size_y;
    domain_size_z = other.domain_size_z;

    atmosphere = other.atmosphere;
    gas_speed_magnitude = other.gas_speed_magnitude;

    return *this;
}

template <typename T>
Domain<T>& Domain<T>::operator=(Domain<T>&& other) {
    domain_size_x = other.domain_size_x;
    domain_size_y = other.domain_size_y;
    domain_size_z = other.domain_size_z;

    atmosphere = other.atmosphere;
    gas_speed_magnitude = other.gas_speed_magnitude;

    other.domain_size_x = 0.0;
    other.domain_size_y = 0.0;
    other.domain_size_z = 0.0;
    other.gas_speed_magnitude = 0.0;

    return *this;
}

template <typename T>
std::vector<Ray<T>> Domain<T>::generate_rays(unsigned long num_rays, int num_children) {
    std::vector<Ray<T>> inflow_rays;

    unsigned long grid_length = static_cast<int>(sqrt(num_rays));

    for(auto i = 0; i < grid_length; i ++) {
        for(auto j = 0; j < grid_length; j ++) {

            T domain_sizes_1[6] = {domain_size_y, domain_size_x, domain_size_x, domain_size_x, domain_size_x, domain_size_y};
            T domain_sizes_2[6] = {domain_size_z, domain_size_z, domain_size_y, domain_size_z, domain_size_y, domain_size_z};
            T domain_sizes_3[6] = {domain_size_x, domain_size_y, domain_size_z, -domain_size_y, -domain_size_z, -domain_size_x};
            int indices_x[6] = {2, 0, 0, 0, 0, 2};
            int indices_y[6] = {0, 2, 1, 2, 1, 0};
            int indices_z[6] = {1, 1, 2, 1, 2, 1};
            int indices_faces[6] = {0, 1, 2, 1, 2, 0};

            for(auto k = 0; k < 6; k ++) {

                Matrix<T> ray_velocity(3, 1, {-gas_speed_magnitude, 0, 0});
                ray_velocity += atmosphere.get_thermal_velocity();
                T m_dot = 1.0 / grid_length / grid_length * (domain_sizes_1[k] * domain_sizes_2[k]) * fabs(ray_velocity(indices_faces[k], 0)) / ray_velocity.norm();
                T ray_pos[4] = {- domain_sizes_1[k] / 2.0 + (i + 0.5) * domain_sizes_1[k] / grid_length,
                            - domain_sizes_2[k] / 2.0 + (j + 0.5) * domain_sizes_2[k] / grid_length,
                            domain_sizes_3[k] / 2.0};

                Matrix<T> ray_origin(3, 1, {ray_pos[indices_x[k]], ray_pos[indices_y[k]], ray_pos[indices_z[k]]});
                if(!(k != 0 && k != 5 && ray_velocity(1, 0) * ray_origin(1, 0) >= 0.0 && ray_velocity(2, 0) * ray_origin(2, 0) >= 0.0)) {
                    Ray<T> ray(ray_origin, ray_velocity, m_dot, num_children, -1);
                    inflow_rays.push_back(ray);
                }
            }
        }
    }
    return inflow_rays;
}

template <typename T>
std::vector<Ray<T>> Domain<T>::generate_rays_parallel(unsigned long num_rays, int start_index, int end_index, int num_children) {
    std::vector<Ray<T>> inflow_rays;

    unsigned long grid_length = static_cast<int>(sqrt(num_rays));
    int start_i = start_index / grid_length, end_i = end_index / grid_length;

    for(auto i = start_i; i <= end_i; i ++) {
        for(auto j = 0; j < grid_length; j ++) {
            if (i * grid_length + j >= start_index && i * grid_length + j < end_index){

                T domain_sizes_1[6] = {domain_size_y, domain_size_x, domain_size_x, domain_size_x, domain_size_x, domain_size_y};
                T domain_sizes_2[6] = {domain_size_z, domain_size_z, domain_size_y, domain_size_z, domain_size_y, domain_size_z};
                T domain_sizes_3[6] = {domain_size_x, domain_size_y, domain_size_z, -domain_size_y, -domain_size_z, -domain_size_x};
                int indices_x[6] = {2, 0, 0, 0, 0, 2};
                int indices_y[6] = {0, 2, 1, 2, 1, 0};
                int indices_z[6] = {1, 1, 2, 1, 2, 1};
                int indices_faces[6] = {0, 1, 2, 1, 2, 0};

                for(auto k = 0; k < 6; k ++) {

                    Matrix<T> ray_velocity(3, 1, {-gas_speed_magnitude, 0, 0});
                    ray_velocity += atmosphere.get_thermal_velocity();
                    T m_dot = 1.0 / grid_length / grid_length * (domain_sizes_1[k] * domain_sizes_2[k]) * fabs(ray_velocity(indices_faces[k], 0)) / ray_velocity.norm();
                    T ray_pos[4] = {- domain_sizes_1[k] / 2.0 + (i + 0.5) * domain_sizes_1[k] / grid_length,
                                - domain_sizes_2[k] / 2.0 + (j + 0.5) * domain_sizes_2[k] / grid_length,
                                domain_sizes_3[k] / 2.0};

                    Matrix<T> ray_origin(3, 1, {ray_pos[indices_x[k]], ray_pos[indices_y[k]], ray_pos[indices_z[k]]});
                    if(!(k != 0 && k != 5 && ray_velocity(1, 0) * ray_origin(1, 0) >= 0.0 && ray_velocity(2, 0) * ray_origin(2, 0) >= 0.0)) {
                        Ray<T> ray(ray_origin, ray_velocity, m_dot, num_children, -1);
                        inflow_rays.push_back(ray);
                    }
                    
                }       
            }
        }
    }
    return inflow_rays;
}