
import matplotlib
matplotlib.use('Agg')  # 'Agg' is a non-interactive backend

import matplotlib.pyplot as plt
import numpy as np
import os

def plot_graph_w_connections(G,df,colors=['y','b','k-'], output_dir ="./", name=""):
    # Plotting the positions
    plt.scatter(df['X'], df['Y'], label='PSII Proteins')
    
    # Plotting circles around each point
    for i in range(len(df)):
        circle = plt.Circle((df['X'][i], df['Y'][i]), radius=7.5, color=colors[0], fill=False)
        plt.gca().add_patch(circle)
        plt.text(df['X'][i], df['Y'][i] + 7.5, str(i), color=colors[1], ha='center', va='bottom')
    
    # Plotting connections
    for edge in G.edges():
        node1, node2 = edge
        x_values = [df['X'][node1], df['X'][node2]]
        y_values = [df['Y'][node1], df['Y'][node2]]
        plt.plot(x_values, y_values, colors[2], alpha=0.5)
    
    # Set plot details
    plt.xlabel('X (nanometers)')
    plt.ylabel('Y (nanometers)')
    plt.title('Connected Points')
    # plt.legend()
    plt.grid(True)
    plt.axis('equal')
    if name: 
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(f"{output_dir}/{name}.png")
    # plt.show()
    plt.clf()

def node_distance(node1, node2, df, show=True):
    distance = np.sqrt((df['X'][node1] - df['X'][node2])**2 + (df['Y'][node1] - df['Y'][node2])**2)
    if show: 
        print(f"The distance between nodes {node1} and {node2} is: {distance}")
    return distance

def node_deltaTheta(node1, node2, df, show=True):
    delta = df["Angle (rad.)"][node1] - df["Angle (rad.)"][node2]
    delta_norm = np.arctan2(np.sin(delta), np.cos(delta))
    if show: 
        print(f"The angle between nodes {node1} and {node2} is: {delta_norm}")
    return delta_norm

def drop_dict_dublicate(old_dict):
    # Create a new dictionary without duplicate keys
    new_dict = {}
    for key in old_dict:
        reversed_key = tuple(reversed(key))
        if not ((reversed_key in new_dict) or (key in new_dict)):
            new_dict[key] = old_dict[key]
        else: continue
    return new_dict

def perpendicular_distance(m, b, x0, y0):
    distance = np.abs(m * x0 - y0 + b) / np.sqrt(m**2 + 1)
    return distance

def perpendicular_distance_df(df, m, b):
    df["distance"] = np.abs(m * df.X - df.Y + b) / np.sqrt(m**2 + 1)
    return df

# Define a function that calculates the distance between two points
def node_distance_(point1, point2):
    x1, y1 = point1
    x2, y2 = point2
    return ((x1 - x2)**2 + (y1 - y2)**2)**0.5

def clean_dict(d):
    # Loop over the items of the dictionary and store the keys and values in a list of tuples
    items = list(d.items())
    # Sort the list of tuples by the values
    items = sorted(items, key=lambda x: x[1])
    # Initialize an empty list to store the unique values, and an empty dictionary to store the final result
    unique = []
    result = {}
    # Loop over the sorted list of tuples
    for key, value in items:
        # Compare each value with the previous ones in the list
        if value in unique or value[::-1] in unique: continue
        found_subset = False
        for i, u in enumerate(unique):
            if is_subset(u, value):
                # Replace u with the current value
                unique[i] = value
                result[key] = value
                found_subset = True
                break
            elif is_subset(value, u):
                found_subset = True
                break

        # If the current value is not a subset of any existing value, add it to the unique list and result dictionary
        if not found_subset:
            unique.append(value)
            result[key] = value

    # Return the final dictionary
    return result

def is_subset(list1, list2):
    # Return True if list1 is a subset of list2, False otherwise
    return all(x in list2 for x in list1)