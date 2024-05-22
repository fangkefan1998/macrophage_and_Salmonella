# BF的意思是这个程序是专门用来进行明场的阈值分割并且产生Mask的
# Import dependencies准备导入必要的库
import numpy as np
from cellpose_omni import models, core

# for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

# 使用jit编译器加快运行速度
# from numba import jit
# @jit
def omni_seg(file_list):
    mpl.rcParams['figure.dpi'] = 1000  # 这个参数可以更改
    plt.style.use('dark_background')

    # 首先建立起来关于目标图像的路径列表，需要将目标图像路径按照列表的形式呈现###########################################################
    files = file_list

    # 开始准备一下读取图像信息
    from cellpose_omni import io, transforms
    from omnipose.utils import normalize99

    imgs = [io.imread(f) for f in files]

    # print some info about the images.
    for i in imgs:
        print('Original image shape:', i.shape)
        print('data type:', i.dtype)
        print('data range:', i.min(), i.max())
    nimg = len(imgs)
    print('number of images:', nimg)

    # fig = plt.figure(figsize=[40] * 2)  # initialize figure
    print('new shape:')
    for k in range(len(imgs)):
        img = transforms.move_min_dim(imgs[k])  # move the channel dimension last
        if len(img.shape) > 2:
            # imgs[k] = img[:,:,1] # could pick out a specific channel
            imgs[k] = np.mean(img, axis=-1)  # or just turn into grayscale

        imgs[k] = normalize99(imgs[k])
        print(imgs[k].shape)


    # 开始准备正式进行分割运算 #################################################################################################
    from cellpose_omni import models
    from cellpose_omni.models import MODEL_NAMES

    MODEL_NAMES
    # 这里可以调整是识别BF还是荧光场的不同的模式 Bact_fluor_omni是荧光场 bact_phase_omni是BF
    model_name = 'bact_phase_omni'
    model = models.CellposeModel(model_type=model_name) # 此处应将GPU去掉，除非电脑上有nVidia的显卡并安装CUDA，否则加上GPU参数会报错


    chans = [0,0] #this means segment based on first channel, no second channel

    n = [0] # make a list of integers to select which images you want to segment
    n = range(nimg) # or just segment them all

    # define parameters
    mask_threshold = -1
    verbose = 0 # turn on if you want to see more output
    #use_gpu = use_GPU #defined above
    transparency = True # transparency in flow output
    rescale=None # give this a number if you need to upscale or downscale your images
    omni = True # we can turn off Omnipose mask reconstruction, not advised
    flow_threshold = 0 # default is .4, but only needed if there are spurious masks to clean up; slows down output
    resample = True #whether or not to run dynamics on rescaled grid or original grid
    cluster=True # use DBSCAN clustering

    masks, flows, styles = model.eval([imgs[i] for i in n],channels=chans,rescale=rescale,mask_threshold=mask_threshold,transparency=transparency,flow_threshold=flow_threshold,omni=omni,cluster=cluster, resample=resample,verbose=verbose)



    io.save_masks(imgs, masks, flows, files,
                  png=False,
                  tif=True, # whether to use PNG or TIF format
                  suffix='', # suffix to add to files if needed
                  save_flows=False, # saves both RGB depiction as *_flows.png and the raw components as *_dP.tif
                  save_outlines=False, # save outline images
                  dir_above=0, # save output in the image directory or in the directory above (at the level of the image directory会保存在样本所在的文件夹)
                  in_folders=True, # save output in folders (recommended)
                  save_txt=False, # txt file for outlines in imageJ
                  save_ncolor=False) # save ncolor version of masks for visualization and editing
    return masks