
/**
 * @brief Represents a Ray object used in simulations for gas dynamics, collisions, and interactions.
 *
 * This class models a ray with properties such as origin, velocity, mass flux, and other relevant data
 * for propagation and collision detection within a geometric space.
 *
 * @tparam T Type of the values used for matrices, mass flux, and other parameters (e.g., float, double).
 */
template <typename T>
class Ray {

    private:
        Matrix<T> origin;             ///< Origin point of the ray (3D vector).
        Matrix<T> velocity;           ///< Velocity vector of the ray (3D vector).
        T mass_flux;                  ///< Mass flux associated with the ray.
        unsigned int children;        ///< Number of child rays generated from this ray (e.g., for splitting).
        int spawn_vertex_index;       ///< Index of the vertex from which the ray was spawned.
        bool failed = false;          ///< Flag indicating if the ray has failed (e.g., during propagation).

    public:
        /**
         * @brief Default constructor for Ray, initializes member variables.
         */
        Ray();

        /**
         * @brief Parameterized constructor to initialize the ray with specific values.
         *
         * @param origin The initial origin of the ray.
         * @param velocity The velocity vector of the ray.
         * @param m_dot The mass flux associated with the ray.
         * @param children The number of children rays spawned.
         * @param spawn_vertex_index Index of the vertex from which the ray is spawned.
         */
        Ray(Matrix<T>& origin, Matrix<T>& velocity, T m_dot, unsigned int children, int spawn_vertex_index);

        /**
         * @brief Copy constructor.
         *
         * @param other Another Ray object to copy from.
         */
        Ray(const Ray<T>& other);

        /**
         * @brief Move constructor.
         *
         * @param other Another Ray object to move from.
         */
        Ray(Ray<T>&& other);

        /**
         * @brief Destructor.
         */
        ~Ray();

        /**
         * @brief Copy assignment operator.
         *
         * @param other Another Ray object to copy from.
         * @return Ray& Reference to the updated object.
         */
        Ray& operator=(const Ray<T>& other);

        /**
         * @brief Move assignment operator.
         *
         * @param other Another Ray object to move from.
         * @return Ray& Reference to the updated object.
         */
        Ray& operator=(Ray<T>&& other);

        /**
         * @brief Getter for the velocity vector.
         *
         * @return Matrix<T>& Reference to the velocity vector.
         */
        Matrix<T>& get_velocity() { return this->velocity; }

        /**
         * @brief Setter for the velocity vector.
         *
         * @param new_velocity The new velocity vector to set.
         */
        void set_velocity(Matrix<T> new_velocity) { this->velocity = new_velocity; }

        /**
         * @brief Getter for the origin point.
         *
         * @return Matrix<T>& Reference to the origin point of the ray.
         */
        Matrix<T>& get_origin() { return this->origin; }

        /**
         * @brief Setter for the origin point.
         *
         * @param new_origin The new origin point to set.
         */
        void set_origin(Matrix<T> new_origin) { this->origin = new_origin; }

        /**
         * @brief Getter for the mass flux.
         *
         * @return T The current mass flux of the ray.
         */
        T get_mass_flux() { return mass_flux; }

        /**
         * @brief Setter for the mass flux.
         *
         * @param mass_flux The new mass flux value to set.
         */
        void set_mass_flux(T mass_flux) { this->mass_flux = mass_flux; }

        /**
         * @brief Checks if the ray collides with a vertex of a triangle surface.
         *
         * @param p1 First vertex of the triangle.
         * @param p2 Second vertex of the triangle.
         * @param p3 Third vertex of the triangle.
         * @param normal Normal vector of the triangle.
         * @param col_position Collision position, set upon detecting a collision.
         * @param distance Distance from the origin to the collision point.
         * @return true If the ray collides with the triangle vertex.
         * @return false If the ray does not collide.
         */
        bool check_vertex_collision(Matrix<T>& p1, Matrix<T>& p2, Matrix<T>& p3, Matrix<T>& normal, Matrix<T>& col_position, T& distance);

        /**
         * @brief Checks if the ray collides with a bounding box.
         *
         * @param bb Bounding box represented by a matrix.
         * @return true If the ray collides with the bounding box.
         * @return false If the ray does not collide.
         */
        bool check_bb_collision(Matrix<T>& bb);

        /**
         * @brief Propagates the ray through the geometry, potentially generating child rays upon collisions.
         *
         * @param geometry The geometry in which the ray is propagating.
         * @param surface The surface with which the ray may collide.
         * @param gas The gas object representing properties of the gas in the system.
         * @param kernel Pointer to the function modeling gas-surface interaction (e.g., CLL model).
         * @param collision_index Pointer to store the index of the collision point.
         * @return std::vector<Ray<T>> A vector of child rays resulting from the propagation.
         */
        std::vector<Ray<T>> propagate(Geometry<T>& geometry, Surface<T>& surface, const Gas<T>& gas, void (*kernel)(T*, T*, T*, const Gas<T>&, Surface<T>&), unsigned int* collision_index);
};


template <typename T>
Ray<T>::Ray() {
    origin = Matrix<T>(3, 1, {0.0, 0.0, 0.0});
    velocity =  Matrix<T>(3, 1, {0.0, 0.0, 0.0});
    mass_flux = 0.0;
    children = 0;
    spawn_vertex_index = -1;
}

template <typename T>
Ray<T>::Ray(Matrix<T>& origin, Matrix<T>& velocity, T m_dot, unsigned int children, int spawn_vertex_index) 
: origin(origin), velocity(velocity), mass_flux(m_dot), children(children), spawn_vertex_index(spawn_vertex_index) {}

template <typename T>
Ray<T>::Ray(const Ray<T>& other) 
: origin(other.origin), velocity(other.velocity), children(other.children), mass_flux(other.mass_flux), spawn_vertex_index(other.spawn_vertex_index) {}

template <typename T>
Ray<T>::Ray(Ray<T>&& other)
: origin(other.origin), velocity(other.velocity), children(other.children), mass_flux(other.mass_flux), spawn_vertex_index(other.spawn_vertex_index) {
    other.children = 0;
    other.mass_flux = 0.0;
    other.spawn_vertex_index = 0;
} 

template <typename T>
Ray<T>::~Ray() {
    children = 0;
    mass_flux = 0.0;
}

template <typename T>
Ray<T>& Ray<T>::operator=(const Ray<T>& other) {
    origin = other.origin;
    velocity = other.velocity;
    mass_flux = other.mass_flux;
    children = other.children;
    spawn_vertex_index = other.spawn_vertex_index;
    return *this;
}

template <typename T>
Ray<T>& Ray<T>::operator=(Ray<T>&& other) {
    origin = other.origin;
    velocity = other.velocity;
    mass_flux = other.mass_flux;
    children = other.children;
    spawn_vertex_index = other.spawn_vertex_index;

    other.children = 0;
    other.mass_flux = 0.0;
    other.children = 0;
    other.spawn_vertex_index = 0;
    return *this;
}

template <typename T>
bool Ray<T>::check_vertex_collision(Matrix<T>& p1, Matrix<T>& p2, Matrix<T>& p3, Matrix<T>& normal, Matrix<T>& col_position, T& distance) {
    
    const T EPSILON = 1e-15;
    Matrix<T> edge1(3, 1), edge2(3, 1), h(3, 1), s(3, 1), q(3, 1);
    Matrix<T> direction(3, 1);
    direction = velocity / velocity.norm(); 
    T a, f, u, v;
    edge1 = p2 - p1;
    edge2 = p3 - p1;
    h = direction.cross(edge2);
    a = edge1.tr().dot(h)(0, 0);
    if (a > -EPSILON && a < EPSILON){
        return false;    // This ray is parallel to this triangle.
    }

    if (normal.tr().dot(direction)(0, 0) > 0){
        return false;
    }

    f = 1.0 / a;
    s = (origin - p1);
    u = f * s.tr().dot(h)(0, 0);

    if (u < - EPSILON || u > 1.0 + EPSILON)
        return false;

    q = s.cross(edge1);
    v = f * direction.tr().dot(q)(0, 0);

    if (v < - EPSILON || u + v > 1.0 + EPSILON)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    T t = f * edge2.tr().dot(q)(0, 0);

    if (t > EPSILON) // ray intersection
    {
        distance = fabs(t);
        col_position = origin + direction * (t);
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}

template <typename T>
std::vector<Ray<T>> Ray<T>::propagate(Geometry<T>& geometry, Surface<T>& surface, const Gas<T>& gas, void (*kernel)(T*, T*, T*, const Gas<T>&, Surface<T>&), unsigned int* collision_index) {

    std::vector<Ray<T>> children_rays;
    Matrix<T> p1_col(3, 1), p2_col(3, 1), p3_col(3, 1), norm_col(3, 1), pos_col(3, 1);
    Matrix<T> p1_x = geometry.get_p1_x(), p1_y = geometry.get_p1_y(), p1_z = geometry.get_p1_z();
    Matrix<T> p2_x = geometry.get_p2_x(), p2_y = geometry.get_p2_y(), p2_z = geometry.get_p2_z();
    Matrix<T> p3_x = geometry.get_p3_x(), p3_y = geometry.get_p3_y(), p3_z = geometry.get_p3_z();
    Matrix<T> norm_x = geometry.get_norm_x(), norm_y = geometry.get_norm_y(), norm_z = geometry.get_norm_z();
    std::vector<Matrix<T>> GSI_properties(geometry.get_GSI_properties());
    std::vector<const char*> GSI_property_names(geometry.get_GSI_property_names());

    T* GSI_prop_col = new T[GSI_properties.size() + 1];

    T min_dist = 1e15;
    unsigned int vertex_index = 0;
    bool bounced = false;

    for(auto i = 0; i < geometry.get_num_vertices(); i ++) {
        Matrix<T> p1(3, 1, {p1_x(i), p1_y(i), p1_z(i)});
        Matrix<T> p2(3, 1, {p2_x(i), p2_y(i), p2_z(i)});
        Matrix<T> p3(3, 1, {p3_x(i), p3_y(i), p3_z(i)});
        Matrix<T> norm(3, 1, {norm_x(i), norm_y(i), norm_z(i)});
        Matrix<T> pos(3, 1);

        T dist = 0.0;
        bool collided = check_vertex_collision(p1, p2, p3, norm, pos, dist);
        if(i != spawn_vertex_index && collided && min_dist > dist) {
            p1_col = p1;
            p2_col = p2;
            p3_col = p3;
            norm_col = norm;
            pos_col = pos;
            vertex_index = i;
            min_dist = dist;
            bounced = true;
            
        }
    }

    if(bounced == false) return children_rays;

    Matrix<T> norm_direction = norm_col;
    Matrix<T> tan1_direction = (velocity - norm_direction * velocity.tr().dot(norm_direction)(0, 0)) / (velocity - norm_direction * velocity.tr().dot(norm_direction)(0, 0)).norm();
    Matrix<T> tan2_direction = norm_col.cross(tan1_direction);

    T vel_in_norm = (velocity.tr().dot(norm_direction))(0, 0);
    T vel_in_tan1 = (velocity.tr().dot(tan1_direction))(0, 0);
    T vel_in_tan2 = (velocity.tr().dot(tan2_direction))(0, 0);

    T vel_in[4] = {vel_in_tan1, vel_in_tan2, vel_in_norm};
    
    for(auto i = 0; i < GSI_properties.size(); i ++) {
        GSI_prop_col[i] = GSI_properties[i](0, vertex_index);
    }

    T v_z = 0.0, v_x = 0.0, angle = 0.0;
    for(auto i = 0; i < children; i ++) {

        T vel_refl[4] = {0.0, 0.0, 0.0};
        kernel(vel_in, vel_refl, GSI_prop_col, gas, surface);
        Matrix<T> child_velocity(3, 1);
        child_velocity = norm_direction * (-vel_refl[2]) * sgn<T>(vel_in[2]) + 
                        tan1_direction * (vel_refl[0]) + 
                        tan2_direction * (vel_refl[1]);

        Ray<T> child_ray(pos_col, child_velocity, mass_flux / children * velocity.norm() / child_velocity.norm(), 1, vertex_index);
        children_rays.push_back(child_ray); 
    }

    *collision_index = vertex_index;
    delete[] GSI_prop_col;

    return children_rays;
}