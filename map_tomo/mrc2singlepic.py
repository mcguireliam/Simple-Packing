import pickle
import numpy as N
import os
import mrcfile

# convert a 3D cube to a 2D image of slices
def cub_img(v, view_dir=2):
    if view_dir == 0:
        vt = N.transpose(v, [1, 2, 0])
    elif view_dir == 1:
        vt = N.transpose(v, [2, 0, 1])
    elif view_dir == 2:
        vt = v

    row_num = vt.shape[0] + 1
    col_num = vt.shape[1] + 1
    slide_num = vt.shape[2]
    disp_len = int(N.ceil(N.sqrt(slide_num)))

    slide_count = 0
    im = N.zeros((row_num * disp_len, col_num * disp_len)) + float('nan')
    for i in range(disp_len):
        for j in range(disp_len):
            im[(i * row_num): ((i + 1) * row_num - 1), (j * col_num): ((j + 1) * col_num - 1)] = vt[:, :, slide_count]
            slide_count += 1

            if (slide_count >= slide_num):
                break

        if (slide_count >= slide_num):
            break

    im_v = im[N.isfinite(im)]

    if im_v.max() > im_v.min():
        im = (im - im_v.min()) / (im_v.max() - im_v.min())

    return {'im': im, 'vt': vt}

# format a 2D array for png saving
def format_png_array(m, normalize=True):
    m = N.array(m, dtype=N.float)

    mv = m[N.isfinite(m)]
    if normalize:
        # normalize intensity to 0 to 1
        if mv.max() - mv.min() > 0:
            m = (m - mv.min()) / (mv.max() - mv.min())
        else:
            m = N.zeros(m.shape)
    else:
        assert mv.min() >= 0
        assert mv.max() <= 1

    m = N.ceil(m * 65534)
    m = N.array(m, dtype=N.uint16)
    # print("max, min:", m.max(), m.min())
    # print("=========")
    # print(m.shape)
    # print("=========")
    return m


def save_png(m, name, normalize=True, verbose=False):
    if verbose:
        print('save_png()')
        print('unique values', sorted(set(m.flatten())))

    m = format_png_array(m, normalize=normalize)

    import png  # in pypng package
    png.from_array(m, mode='L;16').save(name)

def mrc2singlepic(mrcfile, pngdir, pngname='', view_dir=1):
    import iomap as IM
    data = IM.readMrcMap(mrcfile)
    if view_dir == 0:
        data = N.transpose(data, [1, 2, 0])
    elif view_dir == 1:
        data = N.transpose(data, [2, 0, 1])
    elif view_dir == 2:
        data = data

    shape = data.shape
    for j in range(shape[0]):
        d = data[j]
        name = pngdir + pngname + '_' + str(j) + '.png'
        if not os.path.exists(pngdir):
            os.makedirs(pngdir)
        save_png(d, name)
        print('save' + name)

output = {
    'initmap':{'mrc':'IOfile/initmap/mrc/initmap1.mrc','png':'IOfile/initmap/png/initmap1.png','trim':'IOfile/initmap/trim/initmap1T.mrc'},
    'packmap':{'mrc':'IOfile/packmap/mrc/packmap1.mrc','png':'IOfile/packmap/png/packmap1.png','trim':'IOfile/packmap/trim/packmap1T.mrc'},
    'tomo':{'mrc':'IOfile/tomo/mrc/tomo1.mrc','png':'IOfile/tomo/png/tomo1.png','trim':'IOfile/tomo/trim/tomo1T.mrc'},
    'json':'IOfile/json/packing1.json'}

if __name__ == '__main__':
    # in_dirs = {'initmap':'output/initmap/mrc', 'packmap':'output/packmap/mrc', 'tomo':'output/tomo/mrc'}
    # out_dirs = {'initmap':'output/initmap/png', 'packmap':'output/packmap/png', 'tomo':'output/tomo/png'}
    # goal = 'packmap'
    # rootdir = in_dirs[goal]
    # list = os.listdir(rootdir) 
    # for i in range(0,len(list)):
    #     path = os.path.join(rootdir,list[i])
    #     if os.path.isfile(path):
    #         (filename, extension) = os.path.splitext(list[i])
    #         with mrcfile.open(path) as mrc:
    #             data = mrc.data
    #             view_dir = 1
    #             if view_dir == 0:
    #                 data = N.transpose(data, [1, 2, 0])
    #             elif view_dir == 1:
    #                 data = N.transpose(data, [2, 0, 1])
    #             elif view_dir == 2:
    #                 data = data

    #             shape = data.shape
    #             for j in range(shape[0]):
    #                 d = data[j]
    #                 pngdir = out_dirs[goal] + '/' + filename
    #                 pngname = pngdir + '/' + filename + '_' + str(j) + '.png'
    #                 if not os.path.exists(pngdir):
    #                     os.makedirs(pngdir)
    #                 save_png(d, pngname)
    #                 print ('save' + pngname)

    mrc2singlepic('IOfile/tomo/mrc/tomo_SNR04.mrc', 'IOfile/tomo/png/tomo2/SNR04/', 'SNR04')
