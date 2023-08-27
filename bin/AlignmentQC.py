import tifffile as tif
from tifffile import TiffFile, TiffWriter, memmap
import numpy as np
import fire
import os
from csv import writer
import zarr
import csv
import sys
import concurrent.futures
import skimage.registration
import skimage.morphology
import skimage.measure
import skimage.util
import skimage.io
import itertools
import tqdm
import sklearn.linear_model
import os.path
import matplotlib.pyplot as plt
import warnings
import gc
warnings.filterwarnings("ignore")

def DrawRectangle(img, positions, intensity, thickness = 0):
    y0 = positions[0]
    x0 = positions[1]
    y1 = positions[2]
    x1 = positions[3]
    img[y0-thickness:y0, x0:x1] = intensity
    img[y1:y1+thickness, x0:x1] = intensity
    img[y0:y1, x0-thickness:x0] = intensity
    img[y0:y1, x1:x1+thickness] = intensity
    return img

def ColocCoef(img1, img2):
    X = img1.flatten()
    Y = img2.flatten()
    X_bar = np.average(X)
    Y_bar = np.average(Y)
    Pearson = np.sum((X-X_bar)*(Y-Y_bar))/(np.sqrt(np.sum((X-X_bar)**2))*(np.sqrt(np.sum((Y-Y_bar)**2))))
    Manders = np.sum((X)*(Y))/(np.sqrt(np.sum((X)**2))*(np.sqrt(np.sum((Y)**2))))
    del X, Y
    return Pearson, Manders



def GetImgs(filepath, channels):
    Imgs = []
    
    if channels:
        
        try:
            memmap_volume = memmap(filepath)
            img = np.zeros((len(channels), memmap_volume.shape[3], memmap_volume.shape[4]), dtype = 'uint16')
            i=0
            for nch in channels:
                img[i] = memmap_volume[0,nch,0,:,:].copy()
                i+=1
            gc.collect()
        except:
            tiff = tif.TiffFile(filepath)
            #store = tif.imread(filepath, aszarr=True)
            #z = zarr.open(store, mode='r')
            zstore = tiff.aszarr(level=0, key=0)
            immg = zarr.open(zstore)
            #print(z)
            img = np.zeros((len(channels), immg.shape[0], immg.shape[1]), dtype = 'uint16')
            i=0; del immg
            
            for nch in channels:
                zstore = tiff.aszarr(level=0, key=nch)
                img[i] = zarr.open(zstore)
                #print(store)
                
                #img[i] = z[0][nch,:,]
                i+=1
            gc.collect()
        
        
        
    else:
        print('If you use only one file, please specify numbers of channels to use!')
        
            
    return img

def TileAnalysis(img1, img2, TileSize):
    ResList = []
    PearsonC = np.zeros(int(np.floor(img1.shape[0]/TileSize))*int(np.floor(img1.shape[1]/TileSize)))
    for i in range(int(np.floor(img1.shape[0]/TileSize))):
        for j in range(int(np.floor(img1.shape[1]/TileSize))):
            subimg1 = img1[i*TileSize:(i+1)*TileSize, j*TileSize:(j+1)*TileSize]
            subimg2 = img2[i*TileSize:(i+1)*TileSize, j*TileSize:(j+1)*TileSize]
            P, M = ColocCoef(subimg1, subimg2)
            ResList.append([(j+0.5)*TileSize, (i+0.5)*TileSize, P])
            PearsonC[i*int(np.floor(img1.shape[1]/TileSize))+j] = P
    return ResList, PearsonC

def GenerateTileImage(img1, img2, ListPosPear, TileSize):
    border = int(TileSize/20)
    Thickness = int(TileSize/50)
    IMG = np.zeros((3, img1.shape[0], img1.shape[1]), dtype = 'uint16')
    IMG[1,:,:] = img1; IMG[2,:,:] = img2
    MaxIntensity = 65536
    img3 = np.zeros(img1.shape)
    #print(ListPosPear)
    for posP in ListPosPear:
        if not np.isnan(posP[2]):
            x0 = int(posP[0] - TileSize/2 + border); x1 = int(posP[0] + TileSize/2 - border);
            y0 = int(posP[1] - TileSize/2 + border); y1 = int(posP[1] + TileSize/2 - border);
            intensity = int(MaxIntensity*(1-posP[2]))
            if intensity < 1:
                intensity = 1
            img3 = DrawRectangle(img3, [y0, x0, y1, x1], intensity, thickness = Thickness)
    IMG[0,:,:] = img3
    IMG = np.uint16(IMG)
    return IMG
    
def SaveImgPyr(IMG, filepath_out, subresolutions = 2):
    pixelsize = 0.24
    with TiffWriter(filepath_out, bigtiff=True, ome = True) as tif:
        tif.write(
             IMG,
             subifds=subresolutions,
             resolution=(1e4 / pixelsize, 1e4 / pixelsize),
             #metadata=imagej_metadata,
             #**options
        )


        for level in range(subresolutions):
            mag = 2**(level + 1)
            tif.write(
            IMG[:,::mag, ::mag],
            subfiletype=1,
            resolution=(1e4 / mag / pixelsize, 1e4 / mag / pixelsize),
            )


def main_ColocCoef(Imgs, output_folder, PerformTileAnalysis, TileSize, downscale):
    main_csv = os.path.join(output_folder, 'ColocCoef_mean.csv')
    for i in range(len(Imgs)-1):
        for j in range(i+1, len(Imgs)):
            print('Colocalization QC - Cycle ' + str(i) + ' vs Cycle ' + str(j))
            string = 'Cycle' + str(i) + '_vs_Cycle' + str(j)
            if PerformTileAnalysis:
                XYP_List, P = TileAnalysis(Imgs[i], Imgs[j], TileSize)
                
                #saving csv with computed Pearson coefficents
                fields = ['Xpos', 'Ypos', 'Pearson'] 
                filename = string + '.csv'
                output_csv_tiles = os.path.join(output_folder, filename)
                with open(output_csv_tiles, 'w') as f:
                    csv_writer = csv.writer(f)
                    csv_writer.writerow(fields)
                    csv_writer.writerows(XYP_List)
                
                #save average vaues to main file
                if os.path.isfile(main_csv):
                    with open(main_csv, 'a') as f:
                        csv_writer = csv.writer(f) 
                        csv_writer.writerow([string, np.nanmean(P)])
                        f.close()
                else:
                    with open(main_csv, 'w') as f: 
                        csv_writer = csv.writer(f)
                        csv_writer.writerow(['Cycles_compared', 'PearsonCoef'])
                        csv_writer.writerow([string, np.nanmean(P)])
                        f.close()

                    
                #generating and saving image with areas where coef was calculating and intensity of the box ~ Pearson coef
                VisImg = GenerateTileImage(Imgs[i], Imgs[j], XYP_List, TileSize)
                #print(sys.getsizeof(VisImg))
                filename2 = string + '_area_vis.ome.tif'
                filepathout = os.path.join(output_folder, filename2)
                #SaveImgPyr(VisImg, filepathout)
                VisImg_downscaled = VisImg[:, ::downscale, ::downscale]
                tif.imwrite(filepathout, VisImg_downscaled, photometric='rgb')
                del VisImg    
            else:
                P, M = ColocCoef(Imgs[i], Imgs[j])
                Row = []
                Row.append(string)
                Row.append(P)
                Row.append(M)
                if os.path.isfile(output_csv):
                    with open(output_csv, 'a') as f_object:
                        writer_object = writer(f_object)
                        writer_object.writerow(Row)
                        f_object.close() 
                else: 
                    with open(output_csv, 'w') as f_object:
                        writer_object = writer(f_object)
                        writer_object.writerow(['Cycles_compared', 'Pearson', 'Manders'])
                        writer_object.writerow(Row)
                        f_object.close() 

                        
                        
                        
## from here there are functions which are used in OptFlow
def read_tiff_channel(path, channel):
    tiff = tif.TiffFile(path)
    try:
        zstore = tiff.aszarr(level=0, key=channel)
    except IndexError:
        raise ValueError(
            f"Selected channel ({channel}) out of range: {path}"
        ) from None
    img = zarr.open(zstore)
    if img.ndim != 2:
        raise ValueError(f"TIFF must be 2-D or 3-D: {path}")
    return img

def compose_block(base, img, brightness, in_range):
    img = skimage.exposure.rescale_intensity(
        img, in_range=in_range, out_range=float
    )
    block = skimage.img_as_float(base) + img[..., None] * brightness
    block = skimage.img_as_ubyte(np.clip(block, 0, 1))
    base[:] = block


def crop_to(img, shape, offset=None):
    begin = np.array(offset if offset is not None else [0, 0])
    end = begin + shape
    if any(end > img.shape):
        raise ValueError("offset+shape is larger than image size")
    out = img[begin[0]:end[0], begin[1]:end[1]]
    # Above check should have handled this, but let's be sure.
    assert out.shape == tuple(shape)
    return out

def pool_apply_blocked(pool, w, func, img1, img2, desc, *args, **kwargs):
    assert img1.shape[:2] == img2.shape[:2]
    range_y = range(0, img1.shape[0], w)
    range_x = range(0, img2.shape[1], w)
    futures = []
    for y, x in itertools.product(range_y, range_x):
        f = pool.submit(
            func, img1[y:y+w, x:x+w], img2[y:y+w, x:x+w], *args, **kwargs
        )
        futures.append(f)
    #show_progress(futures, desc)
    results = np.array([f.result() for f in futures], dtype=object)
    results = results.reshape(len(range_y), len(range_x), -1)
    #print('pool_apply_blocked')
    return results


def show_progress(futures, desc=None):
    n = len(futures)
    progress = tqdm.tqdm(
        concurrent.futures.as_completed(futures), total=n, desc=desc
    )
    for _ in progress:
        pass


def optical_flow(img1, img2, w, pool):
    assert img1.dtype == img2.dtype
    results = pool_apply_blocked(
        pool, w, skimage.registration.phase_cross_correlation, img1, img2,
        "    computing optical flow", upsample_factor=10,
    )
    shifts = results[..., 0]
    # Unwrap inner-most nested numpy arrays (register_translation shifts).
    shifts = np.stack(shifts.ravel()).reshape(shifts.shape + (-1,))
    # skimage phase correlation tells us how to shift the second image to match
    # the first, but we want the opposite i.e. how the second image has been
    # shifted with respect to the first. So we reverse the shifts.
    shifts = -shifts
    return shifts

def add_vector_plot(img, shifts, pos, output_path):
    my_dpi = 1200 #this defines quality and respectively size of the image
    plt.figure(frameon=False, dpi=my_dpi)
    plt.quiver(pos[:,:,1], pos[:,:,0], shifts[:,:,1], shifts[:,:,0], color='y', scale = 1000)
    plt.imshow(img, cmap='gray')
    plt.show()
    plt.savefig(output_path, dpi = my_dpi,bbox_inches='tight')
    

def build_panel2(
    img1, img2, bmask, w, out_scale, dmax, brightness, intensity_pct, pool, distance_range, folderpath, I, J, bsize):
    assert w % out_scale == 0
    assert np.all(np.mod(img1.shape, w) == 0)
    shifts = optical_flow(img1, img2, w, pool)
    p1 = np.dstack(np.meshgrid(
        range(shifts.shape[0]), range(shifts.shape[1]), indexing='ij'
    ))
    p1 *= w #p1 is like position of tiles
    p1_1 = p1 + bsize/2
    p2 = p1 + shifts #p2 is like updated position of tiles including computed shifts
    lr = sklearn.linear_model.LinearRegression() #just simply linear fit with minimizing square displacements
    lr.fit(p1[bmask], p2[bmask])
    det = np.linalg.det(lr.coef_) #determinant
    assert det != 0, "Degenerate matrix"
    (a, b), (c, d) = lr.coef_ #it is like transf matrix between stilled and moved set of tiles
    #T = [[a b], [c d]
    scale_y = np.linalg.norm([a, b])
    scale_x = det / scale_y
    shear = (a * c + b * d) / det
    rotation = np.arctan2(b, a)
    
    
        
    
    #print("    recovered affine transform: ")
    #print(f"      scale = [{scale_y:.3g} {scale_x:.3g}]")
    #print(f"      shear = {shear:.3g}")
    #print(f"      rotation = {rotation:.3g}")
    #if np.allclose(lr.coef_, np.eye(2), atol=1e-4, rtol=1e-4):
    #    print("      (affine correction is trivial)")
    #else:
    #    print("      (affine correction is non-trivial)")

    shifts = (p2 - lr.intercept_) @ np.linalg.inv(lr.coef_.T) - p1
    
    filepath = folderpath + '/OptFlow_QCpars.csv'
    pairnames = 'Cyc' + str(I) + '_vs_Cyc' + str(J)
    line_csv = [pairnames, np.mean(abs(shifts), axis=(0, 1)), np.median(abs(shifts), axis=(0, 1)), [scale_y, scale_x], shear, rotation]
    if os.path.isfile(filepath):
        with open(filepath,'a') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(line_csv)
            f.close()
            #outfile.writerow("\n")
    else:
        with open(filepath,'w') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(['Cycles_compared','mean_shift', 'median_shift', 'scale', 'shear', 'rotation'])
            csv_writer.writerow(line_csv)
            f.close()


    
    #with np.printoptions(precision=3):
    #    print("    mean shift:", np.mean(shifts, axis=(0, 1)))
    #    print("    median shift:", np.median(shifts, axis=(0,1)))

    angle = np.arctan2(shifts[..., 0], shifts[..., 1]) 
    distance = np.linalg.norm(shifts, axis=2)
    #print("    colorizing")
    heatmap_small = colorize2(distance, distance_range) * bmask[..., None] #colorization on level of tiles
    hs = w // out_scale
    panel = heatmap_small.repeat(hs, axis=1).repeat(hs, axis=0) #extend colorization on whole (downscaled) iamge
    img1_scaled = downscale(img1, out_scale)
    kwargs = {}
    if intensity_pct:
        kwargs["in_range"] = tuple(np.percentile(img1_scaled, intensity_pct))
    compose(panel, img1_scaled, pool, brightness=brightness, **kwargs)
    return panel, distance, shifts, p1_1

def colorize2(distance, distance_range):
    img_lch = np.zeros(distance.shape[:2] + (3,), np.float32)
    MaxInt = 255 #consider 8 bits image
    range_min = distance_range[0]; range_max = distance_range[1]
    #img_lch[..., 0] = 0
    img_lch[..., 1] = (((distance-range_min)/(range_max-range_min))*MaxInt)
    #img_lch[..., 2] = ((distance-range_min)/(range_max-range_min)*MaxInt)
    img = skimage.color.lab2rgb(skimage.color.lch2lab(img_lch))
    img = skimage.img_as_ubyte(img)
    return img

def compose(base, img, pool, brightness, in_range="image"):
    pool_apply_blocked(
        pool,
        1000,
        compose_block,
        base,
        img,
        "    composing image",
        brightness,
        in_range,
    )

def mean_round_sametype(a, axis=None):
    """Compute mean, round, and cast back to original dtype (typically int)."""
    return a.mean(axis=axis).round().astype(a.dtype)


def downscale(image, block_width):
    """Downscale 2D or 3D image using as little extra memory as possible."""
    if image.shape[0] % block_width or image.shape[1] % block_width:
        raise ValueError(
            "`image` width and height must be a multiple of `block_width`."
        )
    block_size = (block_width, block_width)
    if image.ndim == 2:
        pass
    elif image.ndim == 3:
        block_size = block_size + (1,)
    else:
        raise ValueError("`image` must be 2-D or 3-D")
    return block_reduce_nopad(image, block_size, func=mean_round_sametype)

def block_reduce_nopad(image, block_size, func=np.sum, cval=0):
    """Like block_reduce but requires image.shape is multiple of block_size"""
    if len(block_size) != image.ndim:
        raise ValueError(
            "`block_size` must have the same length as `image.shape`."
        )
    # This check lets us skip calling np.pad, which always returns a copy.
    if (np.array(image.shape) % block_size).any():
        raise ValueError(
            "`image.shape` must be an integer multiple of `block_size`."
        )
    blocked = skimage.util.view_as_blocks(image, block_size)
    return func(blocked, axis=tuple(range(image.ndim, blocked.ndim)))


def main_OptFlow(Imgs, output_folder, block_threshold, ga_downscale, bsize):
    area_threshold = 7 #Isolated regions of the image smaller than this number of blocks will be considered background and excluded from analysis 
    output_brightness_scale = 1 #Multiplier for reference image intensity in the output image
    output_contrast_percentile = None #Lower and upper brightness percentile for contrast rescaling. The reference image's brightness and contrast will be scaled to place these values at the bottom and top end of the brightness scale, respectively. If not specified, the image's actual dynamic range will be used."
    out_scale = 20 #Factor by which to downsample the final output image. Must be an integer factor of block-size "
    dmax = 6 #"Large flow magnitudes (outliers) will be clamped down to this value for the purpose of rendering the output image
    distance_range = [0, 1000] #Specify the range of shift distances for visualization Usually it is between 0 and the maximum expected shift value from all images" 
    
    
    for i in range(len(Imgs)-1):
        for j in range(i+1, len(Imgs)):
            print('Optical flow QC - Cycle ' + str(i) + ' vs Cycle ' + str(j))
            img1 = Imgs[i]; img2 = Imgs[j];
            its = np.minimum(img1.shape, img2.shape)
            its_round = its // ga_downscale * ga_downscale
            c1 = crop_to(img1, its_round)
            c2 = crop_to(img2, its_round)
            bmax = skimage.measure.block_reduce(c1, (bsize, bsize), np.max) #Downsample image by    applying function function 'np.max' to local blocks (but other functions like np.median can be used)
            bmask = bmax > block_threshold #here we filter out low intensity areas
            bmask = skimage.morphology.remove_small_objects(bmask, min_size=area_threshold)
            #print("Performing global image alignment")
            r1 = downscale(c1, ga_downscale)
            r2 = downscale(c2, ga_downscale)
            shift = skimage.registration.phase_cross_correlation(
                r1, r2, upsample_factor=ga_downscale, return_error=False,
                )
            shift = (shift * ga_downscale).astype(int)
            border = np.abs(shift)
            offset1 = np.zeros(2, int)
            offset2 = -shift
            for d in 0, 1:
                if offset2[d] < 0:
                    offset1[d] = shift[d]
                    offset2[d] = 0

            shape = (its - border) // bsize * bsize
            c1 = crop_to(img1, shape, offset1)
            c2 = crop_to(img2, shape, offset2)
            #print("Computing foreground image mask")
            bmax = skimage.measure.block_reduce(c1, (bsize, bsize), np.max) #Downsample image by applying function function 'np.max' to local blocks (but other functions like np.median can be used)
            bmask = bmax > block_threshold #here we filter out low intensity areas
            bmask = skimage.morphology.remove_small_objects(bmask, min_size=area_threshold)
            if np.sum(bmask)>3:
                pool = concurrent.futures.ThreadPoolExecutor(len(os.sched_getaffinity(0)))
                panel, dist, shifts, pos_shifts = build_panel2(
                        c1,
                        c2,
                        bmask,
                        bsize,
                        out_scale,
                        dmax,
                        output_brightness_scale,
                        output_contrast_percentile,
                        pool,
                        distance_range,
                        output_folder,
                        i,j,bsize,
                    )
                output_path = output_folder + '/OptFlow_Cyc' + str(i) + '_vs_Cyc' + str(j) + '.tif'
                add_vector_plot(panel, shifts/(2*ga_downscale), pos_shifts/(2*ga_downscale), output_path)
            
            else:
                print('Optical flow for Cycle ' + str(i) + ' vs Cycle ' + str(j) + ' doesnt work due to high OF_BlockThresh')
            #print(panel)
            #print(dist)
            #print(shifts)
            #print("    saving")
           
            #pos_shifts += int(bsize/(2*ga_downscale))

            #x = pos_shifts.reshape(360,2)
            #print(x.shape)
            #np.savetxt('pos_shifts.txt',x)
            
            
            
            
            #skimage.io.imsave(output_path, panel, check_contrast=False)
            
            
def main(filepath, output_folder, channels = None, PerformTileAnalysis = True, TileSize = 5000, ModeColoc = True, OF_BlockThresh = 500, Downscale = 10,  ModeOptFlow = True, OF_TileSize = 2000):
    Imgs = GetImgs(filepath, channels = channels)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if ModeColoc:
        main_ColocCoef(Imgs, output_folder, PerformTileAnalysis, TileSize, Downscale)
    if ModeOptFlow:
        main_OptFlow(Imgs, output_folder, OF_BlockThresh, Downscale, OF_TileSize)
    
        
        
if __name__ == "__main__":
    fire.Fire(main)  
