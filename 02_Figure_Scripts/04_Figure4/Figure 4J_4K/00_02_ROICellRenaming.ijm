// Get the current ROI Manager
roiManager("Select All");
roiCount = roiManager("Count");

// Loop through the ROIs in the ROI Manager and rename them
for (i = 0; i < roiCount; i++) {
    roiManager("Select", i);
    newName = "Cell" + (i + 1);
    roiManager("Rename", newName);
}

print(roiCount);
print("Cells Done");
roiManager("Select All");