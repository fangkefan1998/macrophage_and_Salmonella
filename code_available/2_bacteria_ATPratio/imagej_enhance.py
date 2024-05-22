import imagej

# initialize ImageJ2
ij = imagej.init(mode='headless')
print(f"ImageJ2 version: {ij.getVersion()}")

dataset = ij.IJ.openImage('E:\\Documents\\Python\\Data\\C11_bodipy\\BF\\MG1655-control-BF.tif')


ij.IJ.run(dataset, "Subtract Background...", "rolling=5 light")
ij.IJ.run(dataset, "Enhance Contrast", "saturated=10")

ij.IJ.saveAs(dataset, "Jpeg", "E:/Documents/Python/Data/C11_bodipy/BF/3-6.jpg")


