import packing_single_sphere.simulate as PSSS
import map_tomo.simu_subtomo as SS

num = 1

op = {
    'map':{'situs_pdb2vol_program':'/shared/opt/local/img/em/et/util/situs/Situs_2.7.2/bin/pdb2vol',
           'spacing_s': [10.0], 'resolution_s':[10.0],
           'pdb_dir':'IOfile/pdbfile/',
           'out_file':'/IOfile/map_single/situs_maps.pickle',
           'map_single_path': './IOfile/map_single'}, # read density map from mrc
    'tomo':{'model':{'missing_wedge_angle':30, 'SNR':500000000},
            'ctf':{'pix_size':1.0, 'Dz':-5.0, 'voltage':300, 'Cs':2.0, 'sigma':0.4}},
    'target_size':30
    }


packing_op = {'target': '1bxn',
              'random_protein_number': 4,
              'PDB_ori_path': op['map']['pdb_dir'],
              'iteration':5001,
              'step':1,
              'show_img': 1,
              'show_log': 1
              }

output = {
    'initmap':{
        'mrc':'IOfile/initmap/mrc/initmap{}.mrc'.format(num),
        'png':'IOfile/initmap/png/initmap{}.png'.format(num),
        'trim':'Ofile/initmap/trim/initmap{}T.mrc'.format(num)},
    'packmap':{
        'mrc':'IOfile/packmap/mrc/packmap{}.mrc'.format(num),
        'png':'IOfile/packmap/png/packmap{}.png'.format(num),
        'trim':'IOfile/packmap/trim/packmap{}T.mrc'.format(num),
        'target':{
            'mrc':'IOfile/packmap/target/mrc/packtarget{}.mrc'.format(num),
            'png':'IOfile/packmap/target/png/packtarget{}.png'.format(num)}},
    'tomo':{
        'mrc':'IOfile/tomo/mrc/tomo{}.mrc'.format(num),
        'png':'IOfile/tomo/png/tomo{}.png'.format(num),
        'trim':'IOfile/tomo/trim/tomo{}T.mrc'.format(num),
        'target':{
            'mrc':'IOfile/tomo/target/mrc/tomotarget{}.mrc'.format(num),
            'png':'IOfile/tomo/target/png/tomotarget{}.png'.format(num)}},
    'json':{
        'pack':'../IOfile/json/packing{}.json'.format(num),
        'target':'../IOfile/json/target{}.json'.format(num)}}

target_simu_tomo = SS.simu_subtomo(op, packing_op, output, save_tomo = 0, save_target = 1, save_tomo_slice = 0)

'''
target_simu_tomo = { 'tomo' : subtomogram of target macromolecule, .mrc file
                     'density_map': density map, .mrc file
                     'info': {  'loc' : coordinate of target macromolecule
                                'rotate' : the rotation angle of this macromolecule (ZYZ, Euler angle)
                              }
                    }
                    
the json file of packing result and target will be saved in '../IOfile/json'

save_tomo : 1 save, 0 not save, if the large subtomogram.mrc of the whole packing scene(including 5-10 macromolecules) will be saved
save_target:  1 save, 0 not save, if the subtomogram.mrc of the target will be saved
save_tomo_slice:  1 save, 0 not save, if the sliced image will be saved.
'''