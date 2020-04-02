import sys
import copy
import numpy as N
from scipy import stats


import IMP as I
import IMP.core as IC
import IMP.container as ICO
import IMP.algebra as IAL
import IMP.atom as IAT


def do_packing(conf, op, pymol_file_name, verbose=True):
    m = I.Model()

    box_size = [op['box']['x'], op['box']['y'], op['box']['z']]
    corner0 = IAL.Vector3D(0, 0, 0)
    corner1 = IAL.Vector3D(box_size[0], box_size[1], box_size[2])
    box = IAL.BoundingBox3D(corner0, corner1)

    expand_n = 1
    corner0_expand = IAL.Vector3D(0, 0, 0)
    corner1_expand = IAL.Vector3D(box_size[0] * expand_n, box_size[1] * expand_n, box_size[2] * expand_n)
    box_expand = IAL.BoundingBox3D(corner0_expand, corner1_expand)

    ps = []

    for i, b in enumerate(conf):
        p = I.Particle(m)

        v = IAL.get_random_vector_in(box_expand)  # lsn: radom posotion

        s = IAL.Sphere3D(v, b['r'])  # lsn: position, radius

        pt = IC.XYZR.setup_particle(p, s)  # lsn: particle, parameters

        pt = IAT.Mass.setup_particle(p, b['mass'])

        p.add_attribute(I.FloatKey('vx'), 0.0)  # initial velocity attribute is needed for molecular dynamics simulation
        p.add_attribute(I.FloatKey('vy'), 0.0)  # initial velocity attribute is needed for molecular dynamics simulation
        p.add_attribute(I.FloatKey('vz'), 0.0)

        ps.append(p)

    '''for i, p in enumerate(ps):
        xv = m.get_attribute(I.FloatKey('vx'), p)
        yv = m.get_attribute(I.FloatKey('vy'), p)
        zv = m.get_attribute(I.FloatKey('vz'), p)

        if N.isnan(xv):
            m.set_attribute(I.FloatKey('vx'), p, 0.0)

        if N.isnan(yv):
            m.set_attribute(I.FloatKey('vy'), p, 0.0)

        if N.isnan(zv):
            m.set_attribute(I.FloatKey('vz'), p, 0.0)

        ps[i] = p'''
        
        
    #TODO: Group spheres together for multiBall packing

    lsc_ps = ICO.ListSingletonContainer(m, ps)

    sscell = IC.BoundingBox3DSingletonScore(IC.HarmonicUpperBound(0, 1), box)


    #TODO: New scoring method
    """sf = IC.RestraintsScoringFunction([ICO.SingletonsRestraint(sscell, lsc_ps), IC.ExcludedVolumeRestraint(lsc_ps)])"""

    temprature = op['temprature']

    o = IAT.MolecularDynamics(m)
    o.set_particles(ps)
    o.set_scoring_function(sf)
    o.assign_velocities(temprature)

    s = N.inf

    recent_scores = []

    while (s > op['min score']) and (temprature > 0):
        o.assign_velocities(temprature)
        md = IAT.VelocityScalingOptimizerState(m, ps, temprature)
        o.add_optimizer_state(md)
        s = o.optimize(op['step'])

        cx = N.array(get_center_coordinates(ps))

        is_in = N.array([True] * len(ps))

        for dim_i in range(3):
            is_in[cx[:, dim_i].flatten() < 0] = False
            is_in[cx[:, dim_i].flatten() >= box_size[dim_i]] = False

        print ('\r', 'score:', s, '    ', 'temprature:', temprature, '       ', 'center inside box particles: ', is_in.sum(), '              ',)

        if True:
            recent_scores.append(s)
            while len(recent_scores) > op['recent_scores number']:
                recent_scores.pop(0)

            if len(recent_scores) == op['recent_scores number']:
                x = N.array(range(len(recent_scores)))
                y = N.array(recent_scores)
                recent_scores_slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

                if N.abs(recent_scores_slope) < op['recent_scores slope min']:
                    temprature *= 1-op['temprature decrease']
                    recent_scores = []
        sys.stdout.flush()

    cx = get_center_coordinates(ps)

    conf = copy.deepcopy(conf)
    for i, b in enumerate(conf):
        b['x'] = cx[i]

    return {'conf': conf, 'score': s, 'temprature': temprature, 'inside_box_num': is_in.sum()}


def get_center_coordinates(ps):
    x = []
    for p in ps:
        p = IC.XYZ(p)
        x.append([p.get_x(), p.get_y(), p.get_z()])

    return x
