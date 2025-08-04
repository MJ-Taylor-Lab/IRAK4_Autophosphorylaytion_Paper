run("Grouped Z Project...", "projection=[Max Intensity] group=8314");
run("Median...", "radius=3");
run("Maximum...", "radius=21");
run("Gaussian Blur...", "sigma=51");
run("Median...", "radius=75");
run("mpl-inferno");
