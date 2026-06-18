#include <string>
#include <stdexcept>
#include <cerrno>

#ifdef _WIN32
    #include <direct.h>
#else
    #include <sys/stat.h>
    #include <sys/types.h>
#endif

void create_directory_if_missing(const std::string& path)
{
#ifdef _WIN32
    int status = _mkdir(path.c_str());
#else
    int status = mkdir(path.c_str(), 0755);
#endif

    if (status != 0 && errno != EEXIST) {
        throw std::runtime_error("Could not create directory: " + path);
    }
}

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <exception>
#include <initializer_list>
#include <mpi.h>
#include <tuple>
#include <cmath>
#include <random>
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include <chrono>
#include <thread>
#include "Constants.h"
#include "Matrix.h"
#include "Hermite_Tools.h"
#include "Surface.h"
#include "Gas.h"
#include "PG_Kernel.h"
#include "Random_sampler.h"
#include "Geometry.h"
#include "Ray.h"
#include "Domain.h"
#include "Simulation.h"
#include "Importer.h"
#include "Exporter.h"



int main(int argc, char* argv[]) {

    bool batch = false;
    std::string project_name;
    int sim_index = 0;

    for(auto i = 0; i < argc; i ++) {
        std::string arg = argv[i];
        if(arg == "--batch") batch = true;
        if(arg.find("--project=") != std::string::npos) {
            project_name = arg.substr(10, arg.length());
        }
        if(arg.find("--id=") != std::string::npos) {
            sim_index = std::stoi(arg.substr(5, arg.length()));
        }
    }
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD , &size);

    if(batch) {

        std::tuple<std::string, std::string, std::string, std::string, std::string, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> batch_parameters = import_parameters<double>(project_name);
        std::string param_1_name = std::get<0>(batch_parameters);
        std::string param_2_name = std::get<1>(batch_parameters);
        std::string param_3_name = std::get<2>(batch_parameters);
        std::string param_4_name = std::get<3>(batch_parameters);
        std::string param_5_name = std::get<4>(batch_parameters);
        std::vector<double> param_1 = std::get<5>(batch_parameters);
        std::vector<double> param_2 = std::get<6>(batch_parameters);
        std::vector<double> param_3 = std::get<7>(batch_parameters);
        std::vector<double> param_4 = std::get<8>(batch_parameters);
        std::vector<double> param_5 = std::get<9>(batch_parameters);

        unsigned long s1 = param_1.size(), s2 = param_2.size(), s3 = param_3.size(), s4 = param_4.size(), s5 = param_5.size();


        bool first_line = true;

        if (rank == 0) create_directory_if_missing(project_name + "/Batch_Results");

        long start_index = find_index<double>(project_name, sim_index);

        if(start_index > 0) first_line = false;

        Geometry<double> satellite = import_geometry<double>(project_name);

        for(auto index = start_index; index < s1 * s2 * s3 * s4 * s5; index ++) {

            double start_time, end_time;
            start_time = MPI_Wtime();

            unsigned long i = index / (s2 * s3 * s4 * s5);
            unsigned long j = (index - i * (s2 * s3 * s4 * s5)) / (s3 * s4 * s5);
            unsigned long k = (index - i * (s2 * s3 * s4 * s5) - j * (s3 * s4 * s5)) / (s4 * s5);
            unsigned long l = (index - i * (s2 * s3 * s4 * s5) - j * (s3 * s4 * s5) - k * (s4 * s5)) / s5;
            unsigned long m = (index - i * (s2 * s3 * s4 * s5) - j * (s3 * s4 * s5) - k * (s4 * s5) - l * s5);

            Simulation<double> sim = import_batch<double>(project_name, satellite, i, j, k, l, m, size);

            unsigned long num_vertices = sim.get_geometry().get_num_vertices();

            long num_rays = sim.get_num_rays();
            long local_rays = num_rays / size;

            Matrix<double> force_coeff_vector_local, moment_coeff_vector_local;
            double* normal_coeff_local = new double[num_vertices + 1];
            double* shear_x_coeff_local = new double[num_vertices + 1];
            double* shear_y_coeff_local = new double[num_vertices + 1];
            double* shear_z_coeff_local = new double[num_vertices + 1];
            double* force_coeff = new double[4];
            double* moment_coeff = new double[4];
            double* norm_coeff = new double[num_vertices + 1];
            double* shear_x_coeff = new double[num_vertices + 1];
            double* shear_y_coeff = new double[num_vertices + 1];
            double* shear_z_coeff = new double[num_vertices + 1];

            std::tuple<Matrix<double>, Matrix<double>> aero_coefficients_local = sim.simulate_parallel(rank * local_rays, (rank + 1) * local_rays, normal_coeff_local, shear_x_coeff_local, shear_y_coeff_local, shear_z_coeff_local);
            
            end_time = MPI_Wtime();
            
            force_coeff_vector_local = std::get<0>(aero_coefficients_local);
            moment_coeff_vector_local = std::get<1>(aero_coefficients_local);
            
            MPI_Reduce( force_coeff_vector_local.getData() , force_coeff , 3 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
            MPI_Reduce( moment_coeff_vector_local.getData() , moment_coeff , 3 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
            MPI_Reduce( normal_coeff_local , norm_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
            MPI_Reduce( shear_x_coeff_local , shear_x_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
            MPI_Reduce( shear_y_coeff_local , shear_y_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
            MPI_Reduce( shear_z_coeff_local , shear_z_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

            if(rank == 0) {
                Matrix<double> force_coeff_vector(3, 1, {force_coeff[0], force_coeff[1], force_coeff[2]});
                Matrix<double> moment_coeff_vector(3, 1, {moment_coeff[0], moment_coeff[1], moment_coeff[2]});

                std::cout<<"Simulation for "<<param_1_name<<" = "<<param_1[i]<<", "
                <<param_2_name<<" = "<<param_2[j]<<", "
                <<param_3_name<<" = "<<param_3[k] <<", "
                <<param_4_name<<" = "<<param_4[l] * 180.0 / M_PI <<", "
                <<param_5_name<<" = "<<param_5[m] * 180.0 / M_PI <<" finished in "<<end_time - start_time<<" seconds!\n";
            
                std::ostringstream ss1, ss2, ss3, ss4, ss5;
                ss1 << param_1[i];
                ss2 << param_2[j];
                ss3 << param_3[k] * 180.0 / M_PI ;
                ss4 << param_4[l] * 180.0 / M_PI ;
                ss5 << param_5[m];
                std::string geometry_name = param_1_name + "_" + ss1.str() + "_" + param_2_name + "_" + ss2.str() + "_" + param_3_name + "_" + ss3.str() +  "_" + param_4_name + "_" + ss4.str() +  "_" + param_5_name + "_" + ss5.str() + ".csv";

                // export_geometry<double>(project_name + "/Batch_Results/" + geometry_name, sim.get_geometry(), norm_coeff, shear_x_coeff, shear_y_coeff, shear_z_coeff);

                export_aerodynamic_coefficients_batch_line<double>(project_name + "/Batch_Results/", param_1_name, param_2_name, param_3_name, param_4_name, param_5_name, 
                param_1[i], param_2[j], param_3[k], param_4[l], param_5[m], force_coeff_vector, moment_coeff_vector, first_line, sim_index);
                first_line = false;

            }
            delete[] normal_coeff_local;
            delete[] shear_x_coeff_local;
            delete[] shear_y_coeff_local;
            delete[] shear_z_coeff_local;
            delete[] norm_coeff;
            delete[] shear_x_coeff;
            delete[] shear_y_coeff;
            delete[] shear_z_coeff;
            delete[] force_coeff;
            delete[] moment_coeff;

            normal_coeff_local = nullptr;
            shear_x_coeff_local = nullptr;
            shear_y_coeff_local = nullptr;
            shear_z_coeff_local = nullptr;
            norm_coeff = nullptr;
            shear_x_coeff = nullptr;
            shear_y_coeff = nullptr;
            shear_z_coeff = nullptr;
            force_coeff = nullptr;
            moment_coeff = nullptr;
        }
    }
    else {
        Simulation<double> sim = import_single<double>(project_name, size);
        long num_rays = sim.get_num_rays();
        long local_rays = num_rays / size;

        unsigned long num_vertices = sim.get_geometry().get_num_vertices();

        Matrix<double> force_coeff_vector_local, moment_coeff_vector_local;
        double* normal_coeff_local = new double[num_vertices + 1];
        double* shear_x_coeff_local = new double[num_vertices + 1];
        double* shear_y_coeff_local = new double[num_vertices + 1];
        double* shear_z_coeff_local = new double[num_vertices + 1];
        double* force_coeff = new double[4];
        double* moment_coeff = new double[4];
        double* norm_coeff = new double[num_vertices + 1];
        double* shear_x_coeff = new double[num_vertices + 1];
        double* shear_y_coeff = new double[num_vertices + 1];
        double* shear_z_coeff = new double[num_vertices + 1];

        double start_time, end_time;
        start_time = MPI_Wtime();
        std::tuple<Matrix<double>, Matrix<double>> aero_coefficients_local = sim.simulate_parallel(rank * local_rays, (rank + 1) * local_rays, normal_coeff_local, shear_x_coeff_local, shear_y_coeff_local, shear_z_coeff_local);
        end_time = MPI_Wtime();
        
        
        force_coeff_vector_local = std::get<0>(aero_coefficients_local);
        moment_coeff_vector_local = std::get<1>(aero_coefficients_local);
        
        MPI_Reduce( force_coeff_vector_local.getData() , force_coeff , 3 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
        MPI_Reduce( moment_coeff_vector_local.getData() , moment_coeff , 3 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
        MPI_Reduce( normal_coeff_local , norm_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
        MPI_Reduce( shear_x_coeff_local , shear_x_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
        MPI_Reduce( shear_y_coeff_local , shear_y_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
        MPI_Reduce( shear_z_coeff_local , shear_z_coeff , num_vertices, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);

        if(rank == 0) {
            Matrix<double> force_coeff_vector(3, 1, {force_coeff[0], force_coeff[1], force_coeff[2]});
            Matrix<double> moment_coeff_vector(3, 1, {moment_coeff[0], moment_coeff[1], moment_coeff[2]});

            std::cout<<"Simulation finished in "<<end_time - start_time<<" seconds!\n";
            export_single(project_name, sim, force_coeff_vector, moment_coeff_vector, norm_coeff, shear_x_coeff, shear_y_coeff, shear_z_coeff);
        }
        delete[] normal_coeff_local;
        delete[] shear_x_coeff_local;
        delete[] shear_y_coeff_local;
        delete[] shear_z_coeff_local;
        delete[] norm_coeff;
        delete[] shear_x_coeff;
        delete[] shear_y_coeff;
        delete[] shear_z_coeff;
        delete[] force_coeff;
        delete[] moment_coeff;

        normal_coeff_local = nullptr;
        shear_x_coeff_local = nullptr;
        shear_y_coeff_local = nullptr;
        shear_z_coeff_local = nullptr;
        norm_coeff = nullptr;
        shear_x_coeff = nullptr;
        shear_y_coeff = nullptr;
        shear_z_coeff = nullptr;
        force_coeff = nullptr;
        moment_coeff = nullptr;
    }

    MPI_Finalize();

    return 0;
}