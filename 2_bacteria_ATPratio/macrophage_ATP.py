import cv2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import normalization
import ATPsensor
import Judge
import pandas as pd
# 在进行分析Hslu通道之前先用imageJ进行背景去除，参数设置为半径=10！！！这样可以增强小黑点分析的效果！！
# 指定待分割的图像的全路径
file = ["D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\640\\640-back.tif"]

# 得到整个640通道的所有细菌的mask
mask_640 = ATPsensor.omni_seg(file_list=file)
mask_640 = mask_640[0]

# 分别读入405和488的图像
chan405 = matplotlib.image.imread("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\405\\405.tif")
chan488 = matplotlib.image.imread("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\488\\488.tif")
chan640 = matplotlib.image.imread("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\640\\640-back.tif")

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
kernel = np.ones((1, 1), np.uint8)

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
    result = Judge.Judge_Foci(array=bacteria640, fold=2, adjacent_num=4)
    # 在得到单个mask之后，将其copy一下然后进行腐蚀（进行小黑点判定的mask不需要进行腐蚀！！本来就挺好的）
    bacteria_mask2 = bacteria_mask.copy()
    erosion = cv2.erode(bacteria_mask2, kernel)  # 进行一次腐蚀

    # 将腐蚀之后的mask与各自的荧光通道相乘计算荧光强度
    bacteria405 = chan405*erosion
    bacteria488 = chan488*erosion

    intensity_405 = np.mean(bacteria405[bacteria405 > 1])
    intensity_488 = np.mean(bacteria488[bacteria488 > 1])
    ATPratio = (np.mean(bacteria405[bacteria405 > 1]) - background_405)/(np.mean(bacteria488[bacteria488 > 1]) - background_488)

    if result == 1:
        path = "D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\" + str(n) + "foci.png"
        single_normal = normalization.normalize(bacteria640)
        im = Image.fromarray(single_normal * 255.0)
        im.convert('L').save(path)
        foci_list.append(ATPratio)
        channel_intensity = np.array([[intensity_405, intensity_488]])
        foci_channel = np.append(foci_channel, channel_intensity, axis=0)
    elif result == 0:
        path = "D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\" + str(n) + "even.png"
        single_normal = normalization.normalize(bacteria640)
        im = Image.fromarray(single_normal * 255.0)
        im.convert('L').save(path)
        even_list.append(ATPratio)
        channel_intensity = np.array([[intensity_405, intensity_488]])
        even_channel = np.append(even_channel, channel_intensity, axis=0)
    else:
        pass

foci_output = pd.DataFrame(data=foci_list)
foci_output.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\Foci_atp_ratio.csv", encoding="gbk")
foci_output1 = pd.DataFrame(data=foci_channel)
foci_output1.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\foci\\Foci_atp_channel.csv", encoding="gbk")

even_output = pd.DataFrame(data=even_list)
even_output.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\Even_atp_ratio.csv", encoding="gbk")
even_output1 = pd.DataFrame(data=even_channel)
even_output1.to_csv("D:\\Documents\\PythonProject\\Segmentation\\Data\\ATPsensor\\even\\Even_atp_channel.csv", encoding="gbk")











