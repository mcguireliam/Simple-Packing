# coding:UTF-8
import os
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import pprint


# get atom coord as an array
def get_coord_array (path, file_name):
    '''

    Function: get coord array of all atoms in a pdb file

    :param path: the path of pdb file of all proteins
    :param file_name: the file name of ****.pdb

    :return: atom coord array

                [[x0,y0,z0],
                 [x1,y1,z1],
                    ... ,
                 [xn,yn,zn]]
    '''
    parser = PDBParser(PERMISSIVE=1)
    structure_id = file_name.split('.')[0]
    path_file_name = path + file_name
    structure = parser.get_structure(structure_id, path_file_name)

    atom_coord_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_coord = atom.get_coord()
                    atom_coord_list.append(atom_coord)
    atom_coord_array = np.array(atom_coord_list)
    # print file_name,'atom_coo_array\n',atom_coord_array,'\n'
    # print('get_coord_array DONE!\t', path, ": ", file_name)
    return atom_coord_array

def dist_Eur(vecA,vecB):
    return np.sqrt(sum(np.power((vecA - vecB),2)))

def dist_Eur_array(array,origin):
    power_result = np.power((array - origin), 2)
    sum_power = power_result.sum(axis=1) # axis = 1 row; axis = 0 column
    dist_array = np.sqrt(sum_power)
    maxdist = dist_array.max()
    dist_Eur_dic = {}
    dist_Eur_dic['dist_array'] = dist_array
    dist_Eur_dic['maxdist'] = maxdist
    return dist_Eur_dic



PDB_ori_path = './pdbfile/'
pdb_dict = {}

for file in os.listdir(PDB_ori_path):
    if file != '.DS_Store':

        # get pdb id of each protein
        # print(file)
        pdb_id = file[0:4]

        # get atom number of each protein, may be useful when calculating mass
        atom_coord_array = get_coord_array(PDB_ori_path, file)
        atom_number = len(atom_coord_array[0])

        # get center coordinate of each protein, this is the center of the boundary ball
        xyz_coord_array = atom_coord_array.T
        center = np.mean(xyz_coord_array, axis=1)
        '''
        atom_coord_array              xyz_coord_array
           
        [[x0,y0,z0],               [[x0,x1,x2, ... , xn],
         [x1,y1,z1],       TO       [y0,y1,y2, ... , yn],
             ... ,                  [z0,z1,z2, ... , zn]]
         [xn,yn,zn]]
    
        '''

        # get radius of the boundary ball
        radius = dist_Eur_array(atom_coord_array, center)['maxdist']

        tmp_dict = {}
        tmp_dict['pdb_id'] = pdb_id
        tmp_dict['atom_number'] = atom_number
        tmp_dict['center'] = center
        tmp_dict['radius'] = radius
        '''
        tmp_dict format:
        
        {   'pdbid': ****, 
            'atom_number': ****, 
            'center': ****, 
            'radius': ****
        }
        '''
        # print('pdb 2 single ball: DONE!\t', tmp_dict)
        # print('==========================================\n\n\n')

        # save in a pdb_dict
        pdb_dict[file] = tmp_dict
        '''
        pdb_dict format:

        {
            'protein1':  {'pdbid': ****, 'atom_number': ****, 'center': ****, 'radius': ****}
            'protein2':  {'pdbid': ****, 'atom_number': ****, 'center': ****, 'radius': ****}
            ...
        }
        '''

    else:
        pass


# print the result
dic_print = pprint.PrettyPrinter(indent=4)
dic_print.pprint(pdb_dict)
print('All File Done!')