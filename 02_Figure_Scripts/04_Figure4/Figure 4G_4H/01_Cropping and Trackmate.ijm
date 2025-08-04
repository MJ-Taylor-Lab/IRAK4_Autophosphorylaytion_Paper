// Define the directory and image file
Directory = "/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK2KinaseDomainOnly/Inhibitor20um/20250530 grid03_20nM_cl500_IRAK2KinaseDomain_inhib20uM/20250530 grid03_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007";
IRAK1_KinaseDomainOnly_IMAGE = "20250530 grid03_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007.nd2";

// Construct the full path to the image file
IRAK1_KinaseDomainOnly_Path = Directory + "/" + IRAK1_KinaseDomainOnly_IMAGE;

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
savePath = Directory + "/IRAK2 Kinase Domain only.tif";				//change in line 34 as well

// Save the processed image
saveAs("Tiff", savePath);

// Close the current image window
close();

// ImageJ macro script to process ROIs and TIF files with manual input

roiFile = Directory + "/RoiSet.zip";
tifFile = Directory + "/IRAK2 Kinase Domain only.tif"; 
outputDir = Directory;

// Open the ROI file
run("ROI Manager...");
roiManager("Open", roiFile);

// Open the TIF file
open(tifFile);

// Get the list of ROIs
roiCount = roiManager("count");

// Loop through each ROI to create and save duplicates
for (i = 0; i < roiCount; i++) {
    // Select the ROI
    roiManager("Select", i);
    
    // Create a name for the ROI based on its index
    roiName = "Cell" + (i + 1);  // Adding 1 to start numbering from 1 instead of 0
    
    // Create a directory for the ROI
    roiDir = outputDir + File.separator + roiName;
    File.makeDirectory(roiDir);
    
    // Duplicate all frames of the ROI from the TIF file
    run("Duplicate...", "title=" + roiName + " duplicate range=1-");
    
    // Save the duplicate in the specific folder
    saveAs("Tiff", roiDir + File.separator + roiName + ".tif");
    
        // Run TrackMate
    run('TrackMate', "use_gui=false "
                    + "save_to=["+ roiDir + File.separator + roiName + ".xml] "
                    + "display_results=false "
                    + "radius=.5 "    //Used .3 for IRAK1 Kinase Domain Only
                    + "threshold=10 "    //Used 10 for IRAK1 Kinase Domain Only
                    + "subpixel=true "
                    + "median=true "
                    + "channel=0 "
                    + "max_distance=.5 "
                    + "max_gap_distance=.5 "
                    + "max_frame_gap=1 ");
    
    wait(1000);
    
    // Check if TrackMate results exist
    if (!File.exists(roiDir + File.separator + roiName + ".xml")) {
    	wait(500);
        // If TrackMate results do not exist, skip processing and move to the next ROI
        close(); // Ensure the current image is closed before continuing
        continue;
    }
    
    // Close the duplicate
    close();
}

// Close all open images
run("Close All");