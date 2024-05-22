import cv2
from skimage import measure
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import Judge as jd
from PIL import Image
import normalization


img = cv2.imread("D:\\Foci_Identify\\example\\macrophage_mEcitrine\\macro_mEci1.tif")

# 灰度化
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
# 高斯模糊处理:去噪(效果最好)
blur = cv2.GaussianBlur(gray, (9, 9), 0)
# Sobel计算XY方向梯度
gradX = cv2.Sobel(gray, ddepth=cv2.CV_32F, dx=1, dy=0)
gradY = cv2.Sobel(gray, ddepth=cv2.CV_32F, dx=0, dy=1)
# 计算梯度差
gradient = cv2.subtract(gradX, gradY)
# 绝对值
gradient = cv2.convertScaleAbs(gradient)
# 高斯模糊处理:去噪(效果最好)
blured = cv2.GaussianBlur(gradient, (9, 9), 0)
# 二值化
_, dst = cv2.threshold(blured, 10, 255, cv2.THRESH_BINARY)
# 滑动窗口
kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
# 形态学处理:形态闭处理(腐蚀)
closed = cv2.morphologyEx(dst, cv2.MORPH_CLOSE, kernel)
# 腐蚀与膨胀迭代,迭代的次数可人为调整
closed = cv2.erode(closed, None, iterations=3)
closed = cv2.dilate(closed, None, iterations=1)

# 统计联通区域，lable使用数字标记每个细菌
label, num = measure.label(closed, connectivity=2, background=0, return_num=True) #连通域面积筛选控制在50-250

# 读取原始图像作为矩阵，读取原始的16bit图像
imraw = matplotlib.image.imread('D:\\Foci_Identify\\example\\macrophage_mEcitrine\\macro_mEci1.tif')

# 进行面积筛选
# 生成一个0矩阵作为接下来的mask使用
mask = np.zeros([label.shape[0], label.shape[1]])

# 提示运行进度
#print("Before Mask")

# 新建一个储存每个细菌mask的数组storage，并将每个细菌的mask存储在里面
storage = []
for m in range(1,num+1):  # 如果要计算全部细菌将数字替换为num+1！！！
    if np.sum(label == m) >= 50 and np.sum(label == m) <=175:

        for l in range(0,label.shape[0]):
            for w in range(0,label.shape[1]):
                if label[l,w] == m:
                    mask[l,w] = int(1)
                else:
                    mask[l,w] = int(0)
        storage.append(mask)
    mask = np.zeros([label.shape[0], label.shape[1]])


# 开始对每个符合筛选要求的mask进行小黑点判定，并且用res_list将其储存
res_list = []
for k in range(0, storage.__len__()):
    single = imraw * storage[k]
    # 调用判断小黑点函数,s为结果判断，1代表小黑点，0代表正常，调用jd函数的后面两个参数均可调整
    s = jd.Judge_Foci(array=single, fold=2, adjacent_num=3)
    res_list.append(s)

    if s == 1:
        path = "D:\\Foci_Identify\\example\\macrophage_mEcitrine\\foci\\" + str(k) + "foci.png"
        single_normal = normalization.normalize(single)
        im = Image.fromarray(single_normal * 255.0)
        im.convert('L').save(path)
    elif s == 0:
        path = "D:\\Foci_Identify\\example\\macrophage_mEcitrine\\even\\" + str(k) + "even.png"
        single_normal = normalization.normalize(single)
        im = Image.fromarray(single_normal * 255.0)
        im.convert('L').save(path)
    else:
        pass

# 计算foci细胞的比例
if ([x for x in res_list if x == 1].__len__() + [x for x in res_list if x == 0].__len__()) == 0:
    print("There is no available cells in this image!!!")
elif ([x for x in res_list if x == 1].__len__() + [x for x in res_list if x == 0].__len__()) > 0:
    ratio = sum([x for x in res_list if x == 1])/([x for x in res_list if x == 1].__len__() + [x for x in res_list if x == 0].__len__())
    #ratio = sum(res_list)/res_list.__len__()
    print("The Ratio of foci:", ratio)
else:
    pass










"""
plot1 = imraw * storage[1]
print(jd.Judge_Foci(plot1))
"""

"""
#进入判断小黑点程序
res = 'no'
average = np.sum(plot1)/np.sum(plot1 > 0)

if np.max(plot1) >= average:
    binary = np.int64(plot1 >= 2*average) #形成小黑点的二值图像，这里的几倍于平均值可以调整，控制小黑点筛选的严格程度
    tag, quantity = measure.label(binary, connectivity=2, background=0, return_num=True) #对小黑点的二至图像进行连同区域判断

    for n in range(0, quantity+1):
        if np.sum(tag == n) >= 4 and np.sum(tag == n) <=50: #这里的第一个数字可以调整，第二个数字是为了把大背景去掉用的不必调整
            res = 'yes'
            break
        else:
            pass
print(res)
"""


"""
#找一个图像进行显示
plot1 = imraw * storage[0]
#保存输出这个图像
matplotlib.image.imsave('D:\\Foci_Identify\\out.png', plot1)
"""

"""#显示图像
plt.imshow(plot1, interpolation='nearest', cmap='bone', origin='lower')
plt.colorbar(shrink=.92)  # 设置一个颜色指引条

plt.xticks(())
plt.yticks(())
plt.show()
"""

#显示图像时需要加入保持图像代码
"""
cv2.imshow("hslu-meEcitrine", dst)
cv2.waitKey(0)
cv2.destroyAllWindows()
"""


