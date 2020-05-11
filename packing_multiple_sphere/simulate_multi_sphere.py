import pdb2ball_multiple as P2BM
import sys
sys.path.append("..")
import pprint

op_p2mb = { 'target_protein': '1bxn',
            'PDB_ori_path': '../IOfile/pdbtest/',
            'savepath' :'../IOfile/pdb_multi_sphere/',
            'k_use' : 1,
            'k' : 3,
            'saveORnot' : 1,
            'show_info': 1}




def packing_with_target_mtsp( op_p2mb ):

    # convert pdb file into single ball and get the center and radius of this ball.
    multi_sphere_info = P2BM.pdb2ball_multiple(op_p2mb)
    dic_print = pprint.PrettyPrinter(indent=4)
    dic_print.pprint(multi_sphere_info)

    # set target protein
    print('target protein is', op_p2mb['target_protein'],'\n\n')
    protein_name = []
    protein_name.append(op_p2mb['target_protein'])
    sphere_list = []
    sphere_list.append(multi_sphere_info[protein_name[0]])

    print(sphere_list)

    # # select random proteins
    # random_protein = RS.get_random_protein(boundary_shpere,protein_number = random_protein_number)
    #
    # # get important info
    # info = RS.get_radius_and_id(random_protein, radii_list = radii_list, protein_name = protein_name, show_log = show_log)
    # radius_list = info['radius_list']
    # protein = info['protein_key']

    # set box

    # initialization

    # packing


if __name__ == '__main__':
    packing_with_target_mtsp(op_p2mb)



