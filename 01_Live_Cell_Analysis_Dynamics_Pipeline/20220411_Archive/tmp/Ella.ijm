//Converts nd2 files to 16-bit TIFF. Also opens up as a dialog.

// --- INPUT ---

//Folder containing images
#@ File (label = "Input/Output Directory", style = "directory", value = "~/ImageAnalysisWorkflow/01_TIFF-Subtract") dir

// --- ND2 TO TIFF CONVERSION ---
dir = dir + "/";

ext ="tif";
inList = getFileList(dir);
list = getFromFileList(ext, inList);

setBatchMode(true);

//Image conversion loop
for (i=0; i<list.length; i++) 
{
  inFullname = dir + list[i];
  print("Converting", i+1, "of", list.length, list[i]);

  //open(inFullname);
  
  //Apply function to convert file type
  //run("Grouped Z Project...", "projection=[Max Intensity] group="+getSliceNumber());
  
  close(File.getName(inFullname)+".tif");

  print("...done.");
}


function getFromFileList(ext, fileList)
{
  selectedFileList = newArray(fileList.length);
  ext = toLowerCase(ext);
  j = 0;
  for (i=0; i<fileList.length; i++)
    {
      extHere = toLowerCase(getExtension(fileList[i]));
      if (extHere == ext)
        {
          selectedFileList[j] = fileList[i];
          j++;
        }
    }
  selectedFileList = Array.trim(selectedFileList, j);
  return selectedFileList;
}

// Print array items
function printArray(array)
{ 
  for (i=0; i<array.length; i++)
    print(array[i]);
}

function getExtension(filename)
{
  ext = substring( filename, lastIndexOf(filename, ".") + 1 );
  return ext;
}