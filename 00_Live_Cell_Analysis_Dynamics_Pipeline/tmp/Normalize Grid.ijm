rename("original");
run("Grouped Z Project...", "projection=[Min Intensity] group=169");
imageCalculator("Subtract create stack", "original","MIN_original");
selectWindow("Result of original");
run("Duplicate...", "title=filter duplicate");
run("Maximum...", "radius=51 stack");
run("Gaussian Blur...", "sigma=300 stack");
run("Calculator Plus", "i1=original i2=filter operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=500 k2=0 create");
run("Image Sequence... ", "dir=/Users/u_deliz/Desktop/test/ format=TIFF name=[] start=1");
