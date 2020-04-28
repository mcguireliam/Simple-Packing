from mpl_toolkits.mplot3d import Axes3D, axes3d
import matplotlib.pyplot as plt
import numpy as np
import json




def drawing_function(radius_list,loc_list):


    center = [loc_list[0][0],loc_list[1][0],loc_list[2][0]]
    radius = radius_list[0]
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x_list = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y_list = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z_list = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

    for k in range(1,len(radius_list)):
        center = [loc_list[0][k],loc_list[1][k],loc_list[2][k]]
        radius = radius_list[k]
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
        
        x_list = np.concatenate([x,x_list], axis = 1)
        y_list = np.concatenate([y,y_list], axis = 1)
        z_list = np.concatenate([z,z_list], axis = 1)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(x_list, y_list, z_list,  rstride=4, cstride=4, color='b')
    plt.show()





def get_json_and_plot(path):

    with open('packing.json') as f:
     packing_result = json.load(f)
     pdb_id = packing_result['optimal_result']['pdb_id']
     x_loc = packing_result['optimal_result']['x']
     y_loc = packing_result['optimal_result']['y']
     z_loc = packing_result['optimal_result']['z']


    radius_list = []
    center_list = [x_loc,y_loc,z_loc]
    for ii in range(len(pdb_id)):
        radius_list.extend([packing_result['boundary_shpere'][pdb_id[ii]]['radius'][0]])
    print(radius_list)
    print(center_list)

    drawing_function(radius_list,center_list)



######################### test part ####################
path='packing.json'
get_json_and_plot(path)
    
