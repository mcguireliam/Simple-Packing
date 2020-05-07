# coding:UTF-8
import os
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import sys
sys.path.append("..")
from packing_single_sphere import pdb2ball_single


def k_vaues(atom_number):
    '''
    :param atom_number: the atom number of each macromolecule
    :return: the value of k
    '''
    k_ori = int(math.ceil(atom_number/5000)) + 1
    if k_ori<3:
        k = 3
        # to ensure at lease 3 atom to calculate the rotate angel
    else:
        k = k_ori
    return k

def get_center_and_radius(atom_coord_array, kmeans_result, k, savepath, save_filename, saveORnot = 1, show_info = 1):
    '''

    :param atom_coord_array: the array of the coordinate of all atoms
    :param kmeans_result: the kmeans cluster result of this macrocolecules, using sklearn,cluster
    :param k: the number of cluster
    :param savepath: the path to save output figure and the simiplied pdb file
    :param save_filename: the name of the save file
    :param saveORnot: save the output file or not, 1: save, 0: do not save
    :return: the dictionary save the multiple sphere info of a macromolecules.
            the format of output is:


            sphere_info = {   0: {   'sphere_center': array([36.968, 47.965, -9.969], dtype=float32),
                                     'sphere_id': 0,
                                     'sphere_radius': 39.7004,
                                     'vmd_radius': 26.820240783691407},
                              1: {  'sphere_center': ***,
                                     'sphere_id': ***,
                                     'sphere_radius': ***,
                                     'vmd_radius': ****},
                              2: {}
                            }

    '''
    # get cluster center
    cluster_cent_array = kmeans_result.cluster_centers_

    sphere_info = {}
    i = 0 # cent_id
    while i < k:
        # 输出每个聚类中心处，球的半径（距离该点最远的距离）
        if show_info != 0:
            print('\ncluster number:', i)
            print('cluster_cent:', cluster_cent_array[i])
        cluster_array = get_cluster_data(atom_coord_array, kmeans_result, i)
        dist_array = pdb2ball_single.dist_Eur_array(cluster_array, cluster_cent_array[i])
        max_dist = dist_array['maxdist']
        if show_info != 0:
            print('radius (max_distance):', max_dist)

        sphere_info[i] = {}
        sphere_info[i]['sphere_id'] = i
        sphere_info[i]['sphere_center'] = cluster_array[i]
        sphere_info[i]['sphere_radius'] = max_dist
        sphere_info[i]['vmd_radius'] = (max_dist * 3) / 5 + 3

        # 存储粗粒化 pdb 结果
        if saveORnot != 0:
            save_center2pdb(cluster_cent_array,i, savepath, save_filename)

        i = i + 1
    if saveORnot != 0:
        with open(savepath + save_filename + '_kmeans.pdb', 'a') as f:
            f.write('END\n')
        if show_info != 0:
            print("write file done", savepath, save_filename)

    return sphere_info

def get_cluster_data(atom_coord_array, kmeans_result, cluster_id_number):
    '''

    :param atom_coord_array: the array of the coordinate of all atoms
    :param kmeans_result: the kmeans cluster result of this macrocolecules, using sklearn,cluster
    :param cluster_id_number: int, the id of custer in kmeans_result. It is the same the id of spheres in a macromolecule
    :return: array. the atom coodrinates in a specific cluster
    '''
    # 单独打印某个cluster
    cluster_list = []
    for atom_coo, cluster_id in zip(atom_coord_array, kmeans_result.labels_):
        if cluster_id == cluster_id_number:
            #print data_coo, cluster_id
            cluster_list.append(atom_coo)
    # print 'get cluster data successfully!'
    return np.array(cluster_list)

def save_center2pdb(cluster_array,cent_id, savepath, save_filename):
    '''

    :param cluster_array: array. the atom coodrinates in a specific cluster
    :param cent_id: int, the id of custer in kmeans_result. It is the same the id of spheres in a macromolecule
    :param savepath: the path to save the simplified pdb file. This file will only be used in NAMD simulation.
    :param save_filename: the filename of the saved simplified pdb file
    :return: True
    '''
    if cent_id == 0:
        write_mode = 'w'
    else:
        write_mode = 'a'
    record_type = 'ATOM  '
    atom_id = str(cent_id+1).rjust(5)
    #atom_name = '  N' + save_filename[3:] + str(cent_id+1).ljust(2)
    if save_filename[2:3].isdigit():
        atom_name =  save_filename[1:2] + save_filename[3:] + str(cent_id + 1)
        segment_name = '          ' + save_filename[1:2] + save_filename[3:]
    else:
        atom_name = save_filename[2:] + str(cent_id + 1)
        segment_name = '          ' + save_filename[2:]
    atom_name = atom_name.rjust(5)
    # print(atom_name)
    residue_name = ' ' + save_filename[0:3] + ' '
    residue_id = '    1    '
    xcoo = str(float('%.3f' % cluster_array[cent_id][0])).rjust(8)
    ycoo = str(float('%.3f' % cluster_array[cent_id][1])).rjust(8)
    zcoo = str(float('%.3f' % cluster_array[cent_id][2])).rjust(8)
    occupancy = '  1.00'
    temperature_factor = '  0.00'
    #segment_name = '          S' + save_filename[3:]
    linecont = record_type + atom_id + atom_name + residue_name + residue_id + xcoo + ycoo + zcoo + occupancy + temperature_factor + segment_name + '\n'
    linecont = linecont.upper()
    with open(savepath + save_filename + '_kmeans.pdb', write_mode) as f:
        f.write(linecont)
    # print "write line in:", savepath, save_filename
    return True

def scale_value(array_length):
    '''

    :param array_length: the atom number in a macromolecule
    :return: the value of scale
    '''
    if array_length < 5000:
        scale = 1
    elif array_length < 10000:
        scale = 0.5
    elif array_length < 20000:
        scale = 0.25
    else:
        scale = 0.1
    return scale

def draw_cluster(atom_coord_array, kmeans_result, k, scale, savepath, save_filename):
    '''
    :param atom_coord_array: the array of the coordinate of all atoms
    :param kmeans_result: the kmeans cluster result of this macrocolecules, using sklearn,cluster
    :param k: the number of cluster
    :param scale: int, the scale of a dot in figure.
    :param savepath: the path to save output figure
    :param save_filename: the name of the save file
    :return: True
    '''
    cluster_cent_array = kmeans_result.cluster_centers_

    x = cluster_cent_array[:, 0]
    y = cluster_cent_array[:, 1]
    z = cluster_cent_array[:, 2]
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, s=30, alpha=1, c='r', marker='^', label='center')

    colorlist = ['green', 'blue', 'brown', 'chocolate', 'y', 'darkviolet', 'dimgray', 'magenta', 'crimson',
                 'lightseagreen', 'black']

    i = 0  # cent_id
    while i < k:
        cluster_array = get_cluster_data(atom_coord_array, kmeans_result, i)
        x = cluster_array[:, 0]
        y = cluster_array[:, 1]
        z = cluster_array[:, 2]
        ax.scatter(x, y, z, s=scale, c=colorlist[i], alpha=0.3)
        i = i + 1

    # 添加坐标轴(顺序是Z, Y, X)
    ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red'})
    ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    plt.legend(loc='upper left')
    plt.savefig(savepath + save_filename + '.png', dpi=100)  # # save in specific resolution
    plt.show()

    return True


def pdb2ball_multiple(PDB_ori_path = '../IOfile/pdbtest/', savepath = '../IOfile/pdb_multi_sphere/', k_use = 1, k = 3, saveORnot = 1, show_info = 1):
    '''
     the main function

    :param PDB_ori_path: this is the path that save all the original pdb file
    :param savepath: the path to save output figure and the simiplied pdb file
    :param k_use: use calculated k or not? 1: calculated by atom number, 0: user define, all macromolecules will use the same k number
    :param k: defaule k value
    :return: the dictionary save the mulriple sphere info of a macromolecules.
            the format of output is:

            multi_sphere_info = {  '1f1b': {   0: {  'sphere_center': array([36.968, 47.965, -9.969], dtype=float32),
                                                     'sphere_id': 0,
                                                     'sphere_radius': 39.7004,
                                                     'vmd_radius': 26.820240783691407},

                                               1: {  'sphere_center': ***,
                                                     'sphere_id': ***,
                                                     'sphere_radius': ***,
                                                     'vmd_radius': ****},

                                               2: {  sphere_info }

                                               ...
                                         },

                                 'pdb_id': { 0:{}, 1:{}, 2:{}, ... n:{} },
                                 'pdb_id': { 0:{}, 1:{}, 2:{}, ... n:{} },
                                 ...
                               }

    '''
    multi_sphere_info = {}
    # read all pdb file in a path
    for file in os.listdir(PDB_ori_path):
        if file != '.DS_Store':
            print(file)
            save_filename = file[0:4]
            pdb_id = file[0:4]

            # get all atom coordinate
            atom_coord_array = pdb2ball_single.get_coord_array(PDB_ori_path, file)
            atom_number = len(atom_coord_array)

            # calculate the value of k
            if k_use == 1:
                k = k_vaues(atom_number)
            else:
                pass

            # obtain k-means result
            kmeans_result = KMeans(n_clusters=k).fit(atom_coord_array)
            print('k-means cluster DONE!')

            cluster_cent_array = kmeans_result.cluster_centers_

            #get center and radius, then save as pdb file
            sphere_info = get_center_and_radius(atom_coord_array, kmeans_result, k, savepath, save_filename, saveORnot = saveORnot, show_info= show_info)
            multi_sphere_info[pdb_id] = sphere_info

            # draw figure
            scale = scale_value(len(atom_coord_array))
            draw_cluster(atom_coord_array, kmeans_result, k, scale, savepath, save_filename)


            print('file', file, 'done!')
        else:
            pass

    print('All File Done!')

    return multi_sphere_info


if __name__ == '__main__':
    pdb2ball_multiple()

