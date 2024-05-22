import cv2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import normalization
import ATPsensor
import Judge
import pandas as pd
# 封装为函数，方便批量调用处理，输入：640、405、488三通道；输出：405/488的比值
# 使用jit编译器加快运行速度
# from numba import jit
# @jit(nopython=True)
def Batch_ATP(omni_file="D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\640\\640-back.tif",c405="D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\405\\405.tif",c488="D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\488\\488.tif"):
    # 在进行分析Hslu通道之前先用imageJ进行背景去除，参数设置为半径=10！！！这样可以增强小黑点分析的效果！！
    # 指定待分割的图像的全路径，这里一定要加上[]因为那个函数输入的必须是一个列表！！！直接输入字符串会出错！！！
    file = [omni_file]

    # 得到整个640通道的所有细菌的mask，需要调用已经写好的ATPsensor的函数
    mask_640 = ATPsensor.omni_seg(file_list=file)
    mask_640 = mask_640[0]

    # 分别读入405和488的图像
    chan405 = matplotlib.image.imread(c405)
    chan488 = matplotlib.image.imread(c488)
    chan640 = matplotlib.image.imread(omni_file)

    # 计算405和488的背景荧光强度
    background_mask = mask_640.copy()
    background_mask[background_mask == 0] = 255
    background_mask[background_mask < 255] = 0
    background_mask[background_mask == 255] = 1

    background_405 = np.mean(chan405*background_mask)
    background_488 = np.mean(chan488*background_mask)


    # 遍历mask中的所有细菌，一个一个的判定，
    # 首先找到一共有多少个细菌
    number = mask_640.max()

    # 建立起储存foci或者even的atp ratio的列表
    foci_list = []
    even_list = []

    # 这里建立数组一定是二维的[[]]!!!要不然没法拼接
    foci_channel = np.array([[405, 488]])
    even_channel = np.array([[405, 488]])

    # 建立腐蚀核
    kernel = np.ones((3, 3), np.uint8)

    for n in range(1, number+1):
        # 提示运行进度
        print(n, "/", number)
        # 取出其中一个细菌单独mask，称之为bacteria mask
        bacteria_mask = mask_640*(mask_640 == n)
        # 将其赋值为1，便于接下来的运算
        bacteria_mask[bacteria_mask > 0] = 1
        # 取出一个细菌的640通道，其它部分的细菌都是0！
        bacteria640 = chan640*bacteria_mask
        # 判断是否具有小黑点 并且输出result作为判断的标准
        result = Judge.Judge_Foci(array=bacteria640, fold=2, adjacent_num=4)  # 可以对小黑点的筛选标准进行修改，通过改这两个参数
        # 在得到单个mask之后，将其copy一下然后进行腐蚀（进行小黑点判定的mask不需要进行腐蚀！！本来就挺好的）
        bacteria_mask2 = bacteria_mask.copy()
        erosion = cv2.erode(bacteria_mask2, kernel)  # 进行一次腐蚀

        # 将腐蚀之后的mask与各自的荧光通道相乘计算荧光强度，通过更改erosion和bacteriamask控制是否进行腐蚀！！！！
        bacteria405 = chan405*erosion
        bacteria488 = chan488*erosion

        # 计算出细菌的两个通道的各自的平均荧光强度，用于后续判断与计算
        intensity_405 = np.mean(bacteria405[bacteria405 > 1])
        intensity_488 = np.mean(bacteria488[bacteria488 > 1])
        # ATPratio = (np.mean(bacteria405[bacteria405 > 1]) - background_405)/(np.mean(bacteria488[bacteria488 > 1]) - background_488)

        # 加上这个用于去除405或者488太低的细菌，可以使用这个参数控制荧光的信噪比，然后再限制一下mask的面积把那些特别小的地方去掉
        if (intensity_405 <= 2*background_405) or (intensity_488 <= 2*background_488):
            continue  # 代表跳出本次循环不在把这种细菌计入！！！
        elif bacteria_mask.sum() <= 220:  # 面积小于多少的细菌认为是太小了不计入
            continue
        else:
            pass

        ATPratio = (intensity_405 - background_405) / (intensity_488 - background_488)

        # 等到foci的条件成熟了可以不进行图片的输出，直接算会快一些
        if result == 1:
            # path = "D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\" + str(n) + "foci.png"
            # single_normal = normalization.normalize(bacteria640)
            # im = Image.fromarray(single_normal * 255.0)
            # im.convert('L').save(path)
            foci_list.append(ATPratio)
            channel_intensity = np.array([[intensity_405, intensity_488]])
            foci_channel = np.append(foci_channel, channel_intensity, axis=0)
        elif result == 0:
            # path = "D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\" + str(n) + "even.png"
            # single_normal = normalization.normalize(bacteria640)
            # im = Image.fromarray(single_normal * 255.0)
            # im.convert('L').save(path)
            even_list.append(ATPratio)
            channel_intensity = np.array([[intensity_405, intensity_488]])
            even_channel = np.append(even_channel, channel_intensity, axis=0)
        else:
            pass

    # 将最后的结果储存在一个字典里面并将其输出用作进一步的分析
    final = dict()
    final['foci_list'] = foci_list
    final['even_list'] = even_list
    final['foci_channel'] = foci_channel
    final['even_channel'] = even_channel

    return final


# foci_output = pd.DataFrame(data=foci_list)
# foci_output.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\Foci_atp_ratio.csv", encoding="gbk")
# foci_output1 = pd.DataFrame(data=foci_channel)
# foci_output1.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\Foci_atp_channel.csv", encoding="gbk")

# even_output = pd.DataFrame(data=even_list)
# even_output.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\Even_atp_ratio.csv", encoding="gbk")
# even_output1 = pd.DataFrame(data=even_channel)
# even_output1.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\Even_atp_channel.csv", encoding="gbk")











