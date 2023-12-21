img = rayimage::generate_2d_gaussian(sd = 1, power = 1, dim = c(11,11), width = 3)
# 
ijtiff::display(img)



scale_matrix <- function(m, k) {
  rows = nrow(m)
  cols = ncol(m)
  
  rs <- rep(1:rows, each = k)
  cs <- rep(1:cols, each = k)
  
  return(m[rs, ][, cs])
}
scale_matrix <- cmpfun(scale_matrix)




img <- ijtiff::read_tif("/Users/u_deliz/Desktop/20200110 cl069 1R myd88gfp mScarlett-traf6/Combined_intensity_ref.tif", frame = 1)
img = img[,,,1]

r   = 2.5
x   = 44+.5
y   = 487+.5

h = NROW(img)
w = NCOL(img)


k = 250
d = r*2
new_scale = k*d
circle_mask = (rep(1:new_scale, new_scale)- (r*k))^2 + (rep(1:new_scale, each=new_scale) - (r*k))^2 <= (r*k)^2
circle_mask = matrix(circle_mask, nrow = new_scale)

IntensityFx <- function(img, x, y, r, k){
  # Get parameters
  x_min = floor(x - r)
  x_max = ceiling(x + r)
  y_min = floor(y - r)
  y_max = ceiling(y + r)
  # Crop image
  cropped_img = img[y_min:y_max, x_min:x_max]
  cropped_img = scale_matrix(cropped_img, k)
  # Blow up
  x_min = (x-r) - x_min + 1
  x_min = round(x_min*k)
  x_max = x_min + (r*2*k)-1
  y_min = (y-r) - y_min + 1
  y_min = round(y_min*k)
  y_max = y_min + (r*2*k)-1
  cropped_img = cropped_img[y_min:y_max, x_min:x_max]
  # Mask image
  cropped_img = cropped_img*circle_mask
  # Get intensity
  intensity = sum(cropped_img)/(k^2)
  # display(cropped_img)
  return(intensity)
}
IntensityFx <- cmpfun(IntensityFx)
IntensityFx(img, x, y, r, k)

(1733.462/1731.815)*100-100

library(raster)
library(plotrix)
r1 = raster(img)
width=600
height=600

x <- crop(r1, extent(0,width,0,height))
raster::plot(x)
circlex=44
circley=487
radius=2.5
draw.circle(circlex,circley,radius,border="black")

is_outside_circ = (rasterToPoints(x)[,1] - circlex)^2 + (rasterToPoints(x)[,2] - circley)^2 >= radius^2
x@data@values[ is_outside_circ,] <- 0
raster::plot(x)
draw.circle(circlex,circley,radius,border="blue")


x <- crop(r1, extent(0,25,0,25))












mask = matrix(0, nrow = w, ncol = h)

(x + r/2)
(x - r/2)

(w-x)^2+(h-y)^2 > r^2


mask = 

x = x/w
y = y/h

# img <- melt(img)
# names(img) <- c("x", "y", "z")

(1-x:h-x)^2 + (1-y:w-y)^2 <= r^2

# mask = ((1-x:h-x).' .^2 + (1-y:w-y) .^2) <= r^2


library(raster)
library(plotrix)
r1 <- brick(system.file("external/rlogo.grd", package="raster"))
width=50
height=40
x <- crop(r1, extent(0,width,0,height))
plotRGB(x)
circlex=20
circley=15
radius=5
draw.circle(circlex,circley,radius,border="blue")

is_outside_circ = (rasterToPoints(x)[,1] - circlex)^2 + (rasterToPoints(x)[,2] - circley)^2 >= radius^2
x@data@values[ is_outside_circ,] <- 0
plotRGB(x)
draw.circle(circlex,circley,radius,border="blue")





total_rows, total_cols, total_layers = image_data.shape
X, Y = np.ogrid[:total_rows, :total_cols]
center_row, center_col = total_rows/2, total_cols/2
dist_from_center = (X - total_rows)**2 + (Y - total_cols)**2
radius = (total_rows/2)**2
circular_mask = (dist_from_center > radius)



