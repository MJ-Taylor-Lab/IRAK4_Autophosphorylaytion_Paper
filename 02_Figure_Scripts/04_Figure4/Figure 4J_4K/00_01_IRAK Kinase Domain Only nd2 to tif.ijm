// Define the directory and image file
IRAK1_KinaseDomainOnly_Directory = "/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK2KinaseDomainOnly/Inhibitor20um/20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM/20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_001";
IRAK1_KinaseDomainOnly_IMAGE = "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_001.nd2";

// Construct the full path to the image file
IRAK1_KinaseDomainOnly_Path = IRAK1_KinaseDomainOnly_Directory + "/" + IRAK1_KinaseDomainOnly_IMAGE;

// Open the image
open(IRAK1_KinaseDomainOnly_Path);

// Adjust the window title for further operations
windowTitle = IRAK1_KinaseDomainOnly_IMAGE + " - C=0";
selectWindow(windowTitle);

// Convert to grayscale and invert LUT
run("Grays");
run("Invert LUT");

// Adjust brightness and contrast
setMinAndMax(100, 250); // Adjust these values as needed

// Define the path for saving the processed image
savePath = IRAK1_KinaseDomainOnly_Directory + "/IRAK1 Kinase Domain only.tif";

// Save the processed image
saveAs("Tiff", savePath);

// Close the current image window
close();