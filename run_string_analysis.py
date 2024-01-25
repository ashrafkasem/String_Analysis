from core.utils import *
import argparse
import os
import glob
import pandas as pd
import networkx as nx

# Define the radius and the distance thresholds
radius = 7.5
distance = 25
gap = 35

global_distances = []
global_numbers = []
global_angles = []
global_angles_avg = []
global_name = ""
def parse_args():
    parser = argparse.ArgumentParser(description='Welcome to string analysis framework.')
    parser.add_argument("-data",'--data_dir', help='input directory', metavar='data_dir')
    parser.add_argument("-output",'--output_dir', help="output directory",  metavar='output_dir')
    parser.add_argument('-real', help="Real Data or modeling",  action='store_true')

    args = vars(parser.parse_args())
    return args

if __name__=='__main__':

    args = parse_args()
    
    r2_dir = args["data_dir"].replace("_coordinates","_angles")
    r2_plots = True
    if not os.path.exists(r2_dir): 
        r2_plots = False

    if args["real"]: 
        r2_plots = False
        if 'dark' in args["data_dir"].lower():
            global_name = "dark"
        elif 'light' in args["data_dir"].lower(): 
            global_name = "dark"
        else: 
            raise ValueError(f"Data directory does not contain 'dark' or 'light'.")
    else: 
        # Find the substring "mem_" and get the part after it
        start_index = args["data_dir"].find("mem_") + len("mem_")
        end_index = args["data_dir"].find("_coordinates")
        # Extract the desired substring
        global_name = args["data_dir"][start_index:end_index]

    list_of_files = sorted(glob.glob(f"{args['data_dir']}/*_PSIIs*.csv"))
    if r2_plots: 
        list_of_r2_files = sorted(glob.glob(f"{r2_dir}/*_PSIIs*.csv"))

    outdir = os.path.abspath(args['output_dir'])
    outdir = os.path.join(outdir,list_of_files[0].split("/")[-2])
    if not os.path.exists(f"{outdir}"):
        os.makedirs(f"{outdir}")

    df_dict = {}
    df_r2_dict = {}
    G_dict = {}
    lines_dict = {}
    new_lines_dict = {}
    selected_line_particles_dict = {}
    line_particles_dict = {}

    for csvfile in list_of_files: 
        globalkey = csvfile.split("/")[-1].replace(".csv","")
        # if not globalkey in ["mem_0_60c2s2m2_PSIIs_area67","mem_3_60c2s2m2_PSIIs_area67"]: continue
        df_dict[globalkey] = pd.read_csv(csvfile)

        # convert from Angestrom to nm for KW membrans
        if not args['real']:
            df_dict[globalkey] *= 0.1

        if r2_plots:
            r2_dir 
            csvr2 = csvfile.replace("_coordinates","_angles")
            if os.path.exists(csvr2):
                df_r2_dict[globalkey] = pd.read_csv(csvr2)
                df_r2_dict[globalkey]["Angle (rad.)"] = np.radians(df_r2_dict[globalkey]["Angle (deg.)"])
                df_r2_dict[globalkey]["Angle (rad.)"] = np.arctan2(np.sin(df_r2_dict[globalkey]["Angle (rad.)"]), np.cos(df_r2_dict[globalkey]["Angle (rad.)"]))
            else:
                print(f"{csvr2} file is not existing, will produce only the X-Y plots")
    
    # this part is used to make nx.graph for intial visualization of the membrane (not needed)
    for key, df in df_dict.items(): 
        # if not key in ["mem_0_60c2s2m2_PSIIs_area67","mem_3_60c2s2m2_PSIIs_area67"]: continue
        G_dict[key] = nx.Graph()
        # Add nodes to the graph
        for i in df.index.tolist():
            G_dict[key].add_node(i, pos=(df['X'][i], df['Y'][i]))
        plot_graph_w_connections(G_dict[key],df,output_dir= os.path.join(outdir,"1-graph"), name=key)

        # Loop over all particles
        lines_dict[key] = {}
        for i in df.index.tolist():
            # Select one anchor particle at a time
            anchor = df.loc[i]
            # Look at the next neighbors around the particle
            neighbors = df[((df.X - anchor.X)**2 + (df.Y - anchor.Y)**2) <= distance**2]
            # Fit straight lines between the anchor particle and each of the neighbors
            for j in neighbors.index:
                if i == j: continue
                # print(neighbors, j)
                neighbor = neighbors.loc[j]
                # Calculate the slope and the intercept of the line
                slope = (neighbor.Y - anchor.Y) / (neighbor.X - anchor.X)
                intercept = anchor.Y - slope * anchor.X
                # Store the line parameters in a tuple
                line = (slope, intercept)
                # Add the line to the list of lines
                lines_dict[key][(i,j)]= line

        # drop doublicates
        new_lines_dict[key] = drop_dict_dublicate(lines_dict[key])
        # # Look along each line and see if there are other particles that are at the line within offset of +/- radius
    selected_line_particles_dict[key] = {}
    line_particles_dict[key] = {}

    for key, new_lines in new_lines_dict.items():
        # if not key in ["mem_0_60c2s2m2_PSIIs_area67", "mem_3_60c2s2m2_PSIIs_area67"]: continue
        selected_line_particles = {}
        line_particles = {}
        df = df_dict[key]
        for line in new_lines.items():
            # for testing uncomment the next line
            # if not line[0] == (0, 7): continue
            # Get the line parameters
            slope, intercept = line[1]
            # Find the particles that are close to the line
            line_particles_ = df[[perpendicular_distance(slope, intercept, x, y) <= radius 
                                                    for x, y in zip(df['X'], df['Y'])]]
            # Sort other_points_indices based on the custom key functions
            # Define the lists of indices
            list1 = [line[0][0], line[0][1]]
            list2 = line_particles_.index.to_list()
            list2 = [x for x in list2 if not (x == line[0][0] or x==line[0][1] )]
            
            point1 = (df.loc[list1[0]].X, df.loc[list1[0]].Y)
            point2 = (df.loc[list1[1]].X, df.loc[list1[1]].Y)
            
            # Select the rows that correspond to the indices in list2, and create a new dataframe
            new_df = df.loc[list2]
            
            # Create two new columns that store the distance of each point from point1 and point2
            new_df["distance_from_point1"] = new_df.apply(lambda row: node_distance_(point1, (row.X, row.Y)), axis=1)
            new_df["distance_from_point2"] = new_df.apply(lambda row: node_distance_(point2, (row.X, row.Y)), axis=1)
            
            # Split the new dataframe into two parts, one for the points closer to point 1 than point 2, and one for the points closer to point 2 than point 1
            part1 = new_df[new_df["distance_from_point1"] < new_df["distance_from_point2"]]
            part2 = new_df[new_df["distance_from_point1"] > new_df["distance_from_point2"]]
            
            # Sort each part by the distance from point 1 in descending order, and the distance from point 2 in ascending order, respectively
            part1 = part1.sort_values(by="distance_from_point1", ascending=False)
            part2 = part2.sort_values(by="distance_from_point2", ascending=True)
                
            # Get the list of indices from the final dataframe
            final_list = part1.index.to_list() + list1 + part2.index.to_list()

                
            line_particles[line[0]] = final_list
            # Initialize an empty list to store the selected particles
            selected_line_particles[line[0]] = []
            # loop 1 
            selected_line_particles[line[0]].append(list1[0])
            # Loop over the line particles
            for k in part1.index.to_list():
                # Compare the particle with the last particle in the list
                last_particle = selected_line_particles[line[0]][0]
                # If the particle is gap nm or more far from the last particle, add it to the list
                if np.sqrt((df.loc[k].X - df.loc[last_particle].X)**2 + (df.loc[k].Y - df.loc[last_particle].Y)**2) <= gap:
                    selected_line_particles[line[0]].insert(0, k)

            # loop 2 
            selected_line_particles[line[0]].append(list1[1])
            # Loop over the line particles
            for k in part2.index.to_list():
                # Compare the particle with the last particle in the list
                last_particle = selected_line_particles[line[0]][-1]
                # If the particle is gap nm or more far from the last particle, add it to the list
                if np.sqrt((df.loc[k].X - df.loc[last_particle].X)**2 + (df.loc[k].Y - df.loc[last_particle].Y)**2) <= gap:
                    selected_line_particles[line[0]].append(k)
        selected_line_particles_dict[key] = clean_dict(selected_line_particles)
        line_particles_dict[key] = clean_dict(line_particles)

        # Plot all the selected particles in the df
        plt.scatter(df.X, df.Y, color="gray", label="All particles")
        for i in df.index:
            plt.text(df['X'][i], df['Y'][i] + 7.5, str(i), color='g', ha='center', va='bottom')

        for line in selected_line_particles_dict[key]: 
            slope, intercept = new_lines_dict[key][line]
            list_of_indexs = selected_line_particles_dict[key][line]
            line_particle_df = df.loc[list_of_indexs]
            
            # Plot the selected particles and the line
            plt.scatter(line_particle_df.X, line_particle_df.Y, color="red", label="Selected particles")
            if line_particle_df.shape[0] >= 4: 
                plt.plot(line_particle_df.X, line_particle_df.Y, linestyle='-', marker='o', color='blue', label='Line', markevery=[0, -1])
                for index, row in line_particle_df.iterrows():
                    circle = plt.Circle((row.X, row.Y), radius, color="green", fill=False)
                    plt.gca().add_patch(circle)
        os.makedirs(os.path.join(outdir,'2-graphs_selected_particles'), exist_ok=True)
        plt.savefig(f"{os.path.join(outdir,'2-graphs_selected_particles')}/{key}.png")
        plt.clf()

        # Plot all the particles in the df
        plt.scatter(df.X, df.Y, color="gray", label="All particles")
        for i in df.index:
            plt.text(df['X'][i], df['Y'][i] + 7.5, str(i), color='g', ha='center', va='bottom')

        for line in line_particles_dict[key]: 
            slope, intercept = new_lines_dict[key][line]
            list_of_indexs = line_particles_dict[key][line]
            line_particle_df = df.loc[list_of_indexs]
            
            # Plot the selected particles and the line
            plt.scatter(line_particle_df.X, line_particle_df.Y, color="red", label="Selected particles")
            if line_particle_df.shape[0] >= 4: 
                plt.plot(line_particle_df.X, line_particle_df.Y, linestyle='-', marker='o', color='blue', label='Line', markevery=[0, -1])
                for index, row in line_particle_df.iterrows():
                    circle = plt.Circle((row.X, row.Y), radius, color="green", fill=False)
                    plt.gca().add_patch(circle)
        
        os.makedirs(os.path.join(outdir,'3-graphs_all_particles'), exist_ok=True)
        plt.savefig(f"{os.path.join(outdir,'3-graphs_all_particles')}/{key}.png")
        plt.clf()

        distances = []

        for particle, indices in line_particles_dict[key].items():
            for i in range(len(indices) - 1):
                index1, index2 = indices[i], indices[i + 1]
                distance = node_distance(index1, index2, df, show=False)
                distances.append(distance)
                global_distances.append(distance)
                if distance < 15 : 
                    print(particle)

        # Step 3: Create a histogram
        hist_values, bin_edges, _ = plt.hist(distances, bins=40, density=True, histtype='stepfilled', edgecolor='k')

        # Define the range over which you want to calculate the integral
        lower_bound = 20
        upper_bound = 35

        # Find the indices corresponding to the specified range in the bin edges
        lower_index = np.searchsorted(bin_edges, lower_bound, side='right') - 1
        upper_index = np.searchsorted(bin_edges, upper_bound, side='right') - 1

        # Calculate the integral within the specified range
        integral = np.sum(hist_values[lower_index:upper_index] * np.diff(bin_edges)[lower_index:upper_index])

        plt.title('Distance Distribution between Neighboring Particles')
        plt.xlabel('Distance')
        plt.ylabel('Frequency')
        plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
        plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')

        # Add text annotation for the integral value under the x-axis with the integral symbol
        plt.annotate(fr'$\int_{{{lower_bound}}}^{{{upper_bound}}} f(x) \,dx = {integral:.2f}$',
                    xy=(0.35, 0.8),                # Coordinates for annotation (under x-axis)
                    xycoords='axes fraction',    # Use axes fraction for relative coordinates
                    ha='center', va='center',     # Positioning options
                    color='red')

        # Add arrows to point at the start and end points of the selected area
        plt.annotate('', xy=(upper_bound, max(hist_values)), xytext=(upper_bound, 0),
                    arrowprops=dict(arrowstyle='->', color='blue'))
        plt.annotate('', xy=(lower_bound, max(hist_values)), xytext=(lower_bound, 0),
                    arrowprops=dict(arrowstyle='->', color='blue'))
        # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
        
        os.makedirs(os.path.join(outdir,'4-distances_all_particles'), exist_ok=True)
        plt.savefig(f"{os.path.join(outdir,'4-distances_all_particles')}/{key}.png")
        plt.clf()

        numbers = []

        for particle, indices in line_particles_dict[key].items():
            numbers.append(len(indices))
            global_numbers.append(len(indices))
        # Step 3: Create a histogram
        hist_values, bin_edges, _ = plt.hist(numbers, bins=np.arange(min(numbers), max(numbers) + 2, step=1) , density=True, histtype='stepfilled', edgecolor='k')

        # Define the range over which you want to calculate the integral
        lower_bound = 4
        upper_bound = 9

        # Find the indices corresponding to the specified range in the bin edges
        lower_index = np.searchsorted(bin_edges, lower_bound, side='right') - 1
        upper_index = np.searchsorted(bin_edges, upper_bound, side='right') - 1

        # Calculate the integral within the specified range
        integral = np.sum(hist_values[lower_index:upper_index] * np.diff(bin_edges)[lower_index:upper_index])

        plt.title('Distance Distribution between Neighboring Particles')
        plt.xlabel('Nmuber of connected PSII')
        plt.ylabel('Frequency')
        plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 3')
        plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 4')

        # Add text annotation for the integral value under the x-axis with the integral symbol
        plt.annotate(fr'$\int_{{{lower_bound}}}^{{{upper_bound}}} f(x) \,dx = {integral:.2f}$',
                    xy=(0.65, 0.8),                # Coordinates for annotation (under x-axis)
                    xycoords='axes fraction',    # Use axes fraction for relative coordinates
                    ha='center', va='center',     # Positioning options
                    color='red')

        # Add arrows to point at the start and end points of the selected area
        plt.annotate('', xy=(upper_bound, max(hist_values)), xytext=(upper_bound, 0),
                    arrowprops=dict(arrowstyle='->', color='blue'))
        plt.annotate('', xy=(lower_bound, max(hist_values)), xytext=(lower_bound, 0),
                    arrowprops=dict(arrowstyle='->', color='blue'))
        # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
        os.makedirs(os.path.join(outdir,'5-numbers_all_particles'), exist_ok=True)
        plt.savefig(f"{os.path.join(outdir,'5-numbers_all_particles')}/{key}.png")
        plt.clf()
        angles = []
        if r2_plots and key in df_r2_dict:
            for line, indices in line_particles_dict[key].items():
                line_angles = []
                slope = new_lines_dict[key][line][0]
                line_rotation_to_xaxis = line_rotation_angle_radian(slope)
                for i in indices:
                    delta = df_r2_dict[key]["Angle (rad.)"][i] - line_rotation_to_xaxis
                    delta_norm = np.arctan2(np.sin(delta), np.cos(delta))
                    # index1, index2 = indices[i], indices[i + 1]
                    # angle = node_deltaTheta(index1, index2, df_r2_dict[key], show=False)
                    angles.append(delta_norm)
                    global_angles.append(delta_norm)
                    line_angles.append(delta_norm)
                if line_angles: 
                    global_angles_avg.append(np.mean(line_angles))
                else: 
                    global_angles_avg.append(-1000)
            
            # Step 3: Create a histogram
            hist_values, bin_edges, _ = plt.hist(angles, bins=20, density=True, histtype='stepfilled', edgecolor='k')

            # Define the range over which you want to calculate the integral
            lower_bound = 0
            upper_bound = np.pi/2.0

            plt.title('Angluar Distribution between Neighboring Particles')
            plt.xlabel(r'$\Delta \theta$ PSII axis [rad.]')
            plt.ylabel('Frequency')
            plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
            plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
            plt.axvline(x=upper_bound * 2, color='red', linestyle='--', label='Vertical Line at x = 20')
            plt.axvline(x=-upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
            plt.axvline(x=-upper_bound * 2, color='red', linestyle='--', label='Vertical Line at x = 20')

            # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
            os.makedirs(os.path.join(outdir,'6-angles_all_particles'), exist_ok=True)
            plt.savefig(f"{os.path.join(outdir,'6-angles_all_particles')}/{key}.png")
            plt.clf()
        
    """
    plot the global distance, numbers, angles
    """

    hist_values, bin_edges, _ = plt.hist(global_distances, bins=40, density=True, histtype='stepfilled', edgecolor='k')

    # Define the range over which you want to calculate the integral
    lower_bound = 20
    upper_bound = 35

    # Find the indices corresponding to the specified range in the bin edges
    lower_index = np.searchsorted(bin_edges, lower_bound, side='right') - 1
    upper_index = np.searchsorted(bin_edges, upper_bound, side='right') - 1

    # Calculate the integral within the specified range
    integral = np.sum(hist_values[lower_index:upper_index] * np.diff(bin_edges)[lower_index:upper_index])

    plt.title('Distance Distribution for Neighboring Particles')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
    plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')

    # Add text annotation for the integral value under the x-axis with the integral symbol
    plt.annotate(fr'$\int_{{{lower_bound}}}^{{{upper_bound}}} f(x) \,dx = {integral:.2f}$',
                xy=(0.35, 0.8),                # Coordinates for annotation (under x-axis)
                xycoords='axes fraction',    # Use axes fraction for relative coordinates
                ha='center', va='center',     # Positioning options
                color='red')

    # Add arrows to point at the start and end points of the selected area
    plt.annotate('', xy=(upper_bound, max(hist_values)), xytext=(upper_bound, 0),
                arrowprops=dict(arrowstyle='->', color='blue'))
    plt.annotate('', xy=(lower_bound, max(hist_values)), xytext=(lower_bound, 0),
                arrowprops=dict(arrowstyle='->', color='blue'))
    # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
    
    os.makedirs(os.path.join(outdir,'7-global_vars_all_particles'), exist_ok=True)
    plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/distance_{global_name}.png")
    plt.clf()

    # Step 3: Create a histogram
    hist_values, bin_edges, _ = plt.hist(global_numbers, bins=np.arange(min(global_numbers), max(global_numbers) + 2, step=1) , density=True, histtype='stepfilled', edgecolor='k')

    # Define the range over which you want to calculate the integral
    lower_bound = 4
    upper_bound = 9

    # Find the indices corresponding to the specified range in the bin edges
    lower_index = np.searchsorted(bin_edges, lower_bound, side='right') - 1
    upper_index = np.searchsorted(bin_edges, upper_bound, side='right') - 1

    # Calculate the integral within the specified range
    integral = np.sum(hist_values[lower_index:upper_index] * np.diff(bin_edges)[lower_index:upper_index])

    plt.title('Distance Distribution between Neighboring Particles')
    plt.xlabel('Nmuber of connected PSII')
    plt.ylabel('Frequency')
    plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 3')
    plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 4')

    # Add text annotation for the integral value under the x-axis with the integral symbol
    plt.annotate(fr'$\int_{{{lower_bound}}}^{{{upper_bound}}} f(x) \,dx = {integral:.2f}$',
                xy=(0.65, 0.8),                # Coordinates for annotation (under x-axis)
                xycoords='axes fraction',    # Use axes fraction for relative coordinates
                ha='center', va='center',     # Positioning options
                color='red')

    # Add arrows to point at the start and end points of the selected area
    plt.annotate('', xy=(upper_bound, max(hist_values)), xytext=(upper_bound, 0),
                arrowprops=dict(arrowstyle='->', color='blue'))
    plt.annotate('', xy=(lower_bound, max(hist_values)), xytext=(lower_bound, 0),
                arrowprops=dict(arrowstyle='->', color='blue'))
    # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
    plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/numbers_{global_name}.png")
    plt.clf()

    if r2_plots:
        # Step 3: Create a histogram
        hist_values, bin_edges, _ = plt.hist(global_angles, bins=20, density=True, histtype='stepfilled', edgecolor='k')

        # Define the range over which you want to calculate the integral
        lower_bound = 0
        upper_bound = np.pi/2.0

        plt.title('Angluar Distribution between PSII and strings w.r.t x-axis')
        plt.xlabel(r'$\Delta \theta$ (PSII - line ) [rad.]')
        plt.ylabel('Frequency')
        plt.axvline(x=lower_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
        plt.axvline(x=upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
        plt.axvline(x=upper_bound * 2, color='red', linestyle='--', label='Vertical Line at x = 20')
        plt.axvline(x=-upper_bound, color='red', linestyle='--', label='Vertical Line at x = 20')
        plt.axvline(x=-upper_bound * 2, color='red', linestyle='--', label='Vertical Line at x = 20')

        # print(np.sum(hist_values[0:40] * np.diff(bin_edges)[0:40]))
        plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/angles_{global_name}.png")
        plt.clf()

        # Create a 2D histogram
        plt.hist2d(global_angles_avg, global_numbers, bins=[20, np.arange(min(global_numbers), max(global_numbers) + 2)], cmap='viridis')

        # Add labels and a colorbar
        plt.xlabel(r'avarage $\Delta \theta$ (PSII - line ) [rad.]')
        plt.ylabel('Nmuber of connected PSII')
        plt.title('2D Histogram of Angles vs Distances')
        plt.colorbar(label='Frequency')
        plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/avg_angles_vs_numbers_{global_name}.png")
        plt.clf()

        # Create a 2D histogram
        plt.hist2d(global_angles, [value for value in global_numbers for _ in range(value)], bins=[20, np.arange(min(global_numbers), max(global_numbers) + 2)], cmap='viridis')
        # Add labels and a colorbar
        plt.xlabel(r'$\Delta \theta$ (PSII - line ) [rad.]')
        plt.ylabel('Nmuber of connected PSII')
        plt.title('2D Histogram of Angles vs Distances')
        plt.colorbar(label='Frequency')
        plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/angles_vs_numbers_{global_name}.png")
        plt.clf()

        # # Create a 2D histogram
        # plt.hist2d(global_angles, global_distances, bins=[20, 40], cmap='viridis')
        # # Add labels and a colorbar
        # plt.xlabel(r'$\Delta \theta$ (PSII - line ) [rad.]')
        # plt.ylabel('Distance')
        # plt.title('2D Histogram of Angles vs Distances')
        # plt.colorbar(label='Frequency')
        # plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/angles_vs_distance_{global_name}.png")
        # plt.clf()

    if args['real']: 
        total_number_of_PSII = []
        total_numner_of_PSII_in_strings = []
        for key, df in df_dict.items():
            total_number_of_PSII.append(df.shape[0])
            combined_list = [item for sublist in selected_line_particles_dict[key].values() if len(sublist) > 3 for item in sublist]
            total_numner_of_PSII_in_strings.append(len(list(set(combined_list))))
        
        plt.scatter(total_number_of_PSII, total_numner_of_PSII_in_strings, marker='o', color='blue', label='')

        # Add labels and a title
        plt.xlabel('Total number PSII')
        plt.ylabel('Number PSII in strings')
        plt.title('Sine Wave Plot')
        print(total_number_of_PSII)
        print(total_numner_of_PSII_in_strings)

        # Show the plot
        plt.savefig(f"{os.path.join(outdir,'7-global_vars_all_particles')}/total_PSII_strings_{global_name}.png")
        plt.clf()