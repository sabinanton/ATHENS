import numpy as np
import sys
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mayavi import mlab

def plot_geometry_batch(folder_name, file_name, latitude, longitude, zoom):

    isExist = os.path.exists(folder_name + "/Batch_Plots")
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(folder_name + "/Batch_Plots")

    mesh_data = np.genfromtxt(folder_name + "/Batch_Results/" + file_name, skip_header=1, delimiter=";").T
    p1x, p1y, p1z = mesh_data[1], mesh_data[2], mesh_data[3]
    p2x, p2y, p2z = mesh_data[4], mesh_data[5], mesh_data[6]
    p3x, p3y, p3z = mesh_data[7], mesh_data[8], mesh_data[9]

    norm_array = mesh_data[13]
    shear_x_array = mesh_data[14]
    shear_y_array = mesh_data[15]
    shear_z_array = mesh_data[16]

    # Choose a colormap
    colormap = cm.viridis  # or another colormap you prefer

    norm1 = Normalize(vmin=np.mean(norm_array) - np.var(norm_array), vmax=np.mean(norm_array) + np.var(norm_array))
    norm2 = Normalize(vmin=shear_x_array.min(), vmax=shear_x_array.max())
    norm3 = Normalize(vmin=shear_y_array.min(), vmax=shear_y_array.max())
    norm4 = Normalize(vmin=shear_z_array.min(), vmax=shear_z_array.max())

    norms = [norm1, norm2, norm3, norm4]
    data = [norm_array, shear_x_array, shear_y_array, shear_z_array]
    names = [file_name[:-4] + "_norm.png", file_name[:-4] + "_shear_x.png", file_name[:-4] + "_shear_y.png", file_name[:-4] + "_shear_z.png"]
    labels = ["Pressure Coefficient", "X-Shear Coefficient", "Y-Shear Coefficient", "Z-Shear Coefficient"]

    # Create vertices from the provided arrays
    vertices = np.array([[[p1x[j], p1y[j], p1z[j]],
                        [p2x[j], p2y[j], p2z[j]],
                        [p3x[j], p3y[j], p3z[j]]] for j in range(len(p1x))])


    for i in range(4) :

        # Create a figure and a 3D axis
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

       
        # Create the polygons for coloring
        for j in range(len(vertices)):
            colored_polygon = Poly3DCollection([vertices[j]], 
                                            facecolors=colormap(norms[i](data[i][j])),
                                            edgecolor='black',
                                            linewidth=0.2) # Avoid drawing edges here
            ax.add_collection3d(colored_polygon)
            
        # Auto-scale the axis to the mesh size
        ax.autoscale_view()
        ax.view_init(elev=latitude, azim=longitude)

        frame_width = max(p1x.max(), p1y.max(), p1z.max(), p2x.max(), p2y.max(), p2z.max(), p3x.max(), p3y.max(), p3z.max()) - \
                    min(p1x.min(), p1y.min(), p1z.min(), p2x.min(), p2y.min(), p2z.min(), p3x.min(), p3y.min(), p3z.min())


        ax.set_xlim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary
        ax.set_ylim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary
        ax.set_zlim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary

        fig.tight_layout()
        
        m = cm.ScalarMappable(cmap=colormap)
        m.set_array(data[i])
        cax = fig.add_axes([0.1, 0.98, 0.8, 0.02])
        fig.colorbar(m, orientation='horizontal', cax=cax, label=labels[i])

        # Save the plot
        plt.savefig(folder_name + "/Batch_Plots/" + names[i], dpi=500)
        plt.clf()
        plt.close(fig)


def plot_geometry_single(folder_name, latitude, longitude, zoom):
    isExist = os.path.exists(folder_name + "/Single_Plots")
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(folder_name + "/Single_Plots")

    mesh_data = np.genfromtxt(folder_name + "/Single_Results/geometry.csv", skip_header=1, delimiter=";").T
    p1x, p1y, p1z = mesh_data[1], mesh_data[2], mesh_data[3]
    p2x, p2y, p2z = mesh_data[4], mesh_data[5], mesh_data[6]
    p3x, p3y, p3z = mesh_data[7], mesh_data[8], mesh_data[9]

    norm_array = mesh_data[13]
    shear_x_array = mesh_data[14]
    shear_y_array = mesh_data[15]
    shear_z_array = mesh_data[16]

    # Choose a colormap
    colormap = cm.viridis  # or another colormap you prefer

    norm1 = Normalize(vmin=0, vmax=0.05)
    norm2 = Normalize(vmin=0, vmax=0.05)
    norm3 = Normalize(vmin=0, vmax=0.05)
    norm4 = Normalize(vmin=0, vmax=0.05)

    norms = [norm1, norm2, norm3, norm4]
    data = np.array([norm_array, shear_x_array, shear_y_array, shear_z_array])
    names = ["geometry" + "_norm.png", "geometry" + "_shear_x.png", "geometry" + "_shear_y.png", "geometry" + "_shear_z.png"]
    labels = ["Pressure Coefficient", "X-Shear Coefficient", "Y-Shear Coefficient", "Z-Shear Coefficient"]

    # Consolidate vertex coordinates
    vertices = np.column_stack((p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z)).reshape(-1, 3)
    N = len(p1x)
    data_3d = np.column_stack((data[0], data[0], data[0])).flatten()

    # Create a connectivity array
    # For N triangles, the connectivity array has shape (N, 3)
    # Each row in the connectivity array represents a triangle,
    # with each value being an index into the vertices array
    # triangles = np.arange(3*N).reshape(-1, 3)
    # print(triangles.shape, data[0].shape)

    # # Plot the mesh
    # mlab.figure(bgcolor=(1, 1, 1))
    # mlab.triangular_mesh(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles, scalars=data_3d, colormap="viridis", vmin=0, vmax=2, line_width=0.5)
    # mlab.show()


    # Create vertices from the provided arrays
    vertices = np.array([[[p1x[j], p1y[j], p1z[j]],
                        [p2x[j], p2y[j], p2z[j]],
                        [p3x[j], p3y[j], p3z[j]]] for j in range(len(p1x))])


    for i in range(4) :

        # Create a figure and a 3D axis
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

       
        # Create the polygons for coloring
        for j in range(len(vertices)):
            colored_polygon = Poly3DCollection([vertices[j]], 
                                            facecolors=colormap(norms[i](data[i][j])),
                                            edgecolor='black',
                                            linewidth=0.1) # Avoid drawing edges here
            ax.add_collection3d(colored_polygon)
            
        # Auto-scale the axis to the mesh size
        ax.autoscale_view()
        ax.view_init(elev=latitude, azim=longitude)

        frame_width = max(p1x.max(), p1y.max(), p1z.max(), p2x.max(), p2y.max(), p2z.max(), p3x.max(), p3y.max(), p3z.max()) - \
                    min(p1x.min(), p1y.min(), p1z.min(), p2x.min(), p2y.min(), p2z.min(), p3x.min(), p3y.min(), p3z.min())


        ax.set_xlim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary
        ax.set_ylim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary
        ax.set_zlim([-frame_width / 2 / zoom, frame_width / 2 / zoom])  # assuming the points are all positive; adjust if necessary

        fig.tight_layout()

        m = cm.ScalarMappable(cmap=colormap)
        m.set_array([0, 0.001])
        cax = fig.add_axes([0.1, 0.98, 0.8, 0.02])
        fig.colorbar(m, orientation='horizontal', cax=cax, label=labels[i])

        # Save the plot
        plt.savefig(folder_name + "/Single_Plots/" + names[i], dpi=500)
        plt.clf()
        plt.close(fig)


def plot_geometries(folder_name, latitude, longitude, zoom):

    files = [f for f in listdir(folder_name + "/Batch_Results") if isfile(join(folder_name + "/Batch_Results", f))]
    for file in files:
        if file != "aerodynamic_coefficients.csv":
            plot_geometry_batch(folder_name, file, latitude, longitude, zoom)

if __name__ == "__main__":
    if sys.argv[2] == "True" or sys.argv[2] == "true":
        plot_geometries(sys.argv[1], float(sys.argv[3]), float(sys.argv[4]))
    else:
        plot_geometry_single(sys.argv[1], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

