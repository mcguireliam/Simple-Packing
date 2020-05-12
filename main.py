import packing_single_sphere.simulate as PSSS
import map_tomo.simu_subtomo as SS

c

target_simu_tomo = SS.simu_subtomo(op, packing_op, output, save_tomo = 0, save_target = 1, save_tomo_slice = 0)

'''
target_simu_tomo = { 'tomo' : subtomogram of target macromolecule, .mrc file
                     'density_map': density map, .mrc file
                     'info': {  'loc' : coordinate of target macromolecule
                                'rotate' : the rotation angle of this macromolecule (ZYZ, Euler angle)
                                'name' : the name of a macromolecule
                              }
                    }
                    
the json file of packing result and target will be saved in '../IOfile/json'

save_tomo : 1 save, 0 not save, if the large subtomogram.mrc of the whole packing scene(including 5-10 macromolecules) will be saved
save_target:  1 save, 0 not save, if the subtomogram.mrc of the target will be saved
save_tomo_slice:  1 save, 0 not save, if the sliced image will be saved.
'''