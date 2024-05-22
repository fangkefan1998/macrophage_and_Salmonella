#定义判断foci函数，这里的array可以看作plot1
import numpy as np
from skimage import measure

# 使用jit编译器加快运行速度
# from numba import jit
# @jit
def Judge_Foci(array, fold = 2, adjacent_num = 4):
    res = 0 #No!!!

    # 计算细菌部分的平均荧光强度
    average = np.sum(array)/np.sum(array > 0)

    # 如果整个矩阵的最大值超过fold倍的平均荧光强度，才有资格进行小黑点判断，否则一律没有
    if np.max(array) >= fold * average:
        binary = np.int64(array >= fold * average)  # 形成小黑点的二值图像，这里和上面的几倍于平均值可以调整，控制小黑点筛选的严格程度
        tag, quantity = measure.label(binary, connectivity=2, background=0, return_num=True)  # 对小黑点的二至图像进行连同区域判断

        for n in range(0, quantity + 1):
            if np.sum(tag == n) >= adjacent_num and np.sum(tag == n) <= 15:  # 这里的第一个数字可以调整，第二个数字是为了把大背景去掉用的不必调整
                res = 1  # Yes!!!
                break
            elif np.sum(tag == n) < adjacent_num:
                res = 0  # 表示不确定是否有小黑点：unknow!!!
            else:
                res = 0.5  # 表示不确定是否有小黑点：unknow!!!
    else:
        pass
    return res


