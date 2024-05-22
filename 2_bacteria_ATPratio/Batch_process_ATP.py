import cv2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import normalization
import ATPsensor
import Judge
import pandas as pd
import Function_macrophage_ATP

# 先写出每个图像的文件名称，一定按照顺序！！！
# 注意这个640通道代表小黑点通道，对于514和561全部适用！！！
list_640 = ["561-back.tif"]
list_405 = ["405.tif"]
list_488 = ["488.tif"]

file_640 = "D:\\DataAnalysis\\Python\\test\\561\\"
file_405 = "D:\\DataAnalysis\\Python\\test\\405\\"
file_488 = "D:\\DataAnalysis\\Python\\test\\488\\"

# 先建立一个大的list用来存放计算得到的所有细菌的ATPratio
foci_ATP = []
even_ATP = []
foci_channel = []
even_channel = []

for n in range(0, list_640.__len__()):
    # 先写出需要分析的某一张图片的全路径
    path_640 = file_640 + list_640[n]
    path_405 = file_405 + list_405[n]
    path_488 = file_488 + list_488[n]

    # 记住final-result输出的是一个字典！按照字典的格式进行调用！
    final_result = Function_macrophage_ATP.Batch_ATP(omni_file=path_640, c405=path_405, c488=path_488)
    foci_ATP = foci_ATP + final_result['foci_list']
    even_ATP = even_ATP + final_result['even_list']
    foci_channel.append(final_result['foci_channel'])
    even_channel.append(final_result['even_channel'])

# 接下来保存那些需要保存的数据与文件
foci_output = pd.DataFrame(data=foci_ATP)
foci_output.to_csv("D:\\DataAnalysis\\Python\\test\\Foci.csv", encoding="gbk")

even_output = pd.DataFrame(data=even_ATP)
even_output.to_csv("D:\\DataAnalysis\\Python\\test\\Even.csv", encoding="gbk")