#!/usr/bin/env python3

############################################################
# Script to generate kmz from gdal readable binary files.  #
# Modified after save_kmz.py from  mintpy                  #
# https://github.com/insarlab/MintPy                       #
# Author: Talib Oliver, 2023 (Caltech/JPL)                 #
############################################################

import os
import glob
import shutil
import argparse
from zipfile import ZipFile
from osgeo import gdal
import numpy as np
from lxml import etree
from matplotlib import colorbar, colors, pyplot as plt, ticker
from pykml.factory import KML_ElementMaker as KML

def createParser():
    EXAMPLE = """example:
      gdal2kmz.py -i file_to_conver.dat -o out_file.kmz --zero-mask -v -10 10 --dpi 150 
    """

    parser = argparse.ArgumentParser(description='Script to generate kmz from gdal readable binary files.',
                                     formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)
    
    parser.add_argument('-i', '--input', dest='in_file', type=str, required=True,
            help='Input gdal readable raster.')
    parser.add_argument('-b', '--band', dest='band', metavar='NUM', type=int, default=1,
            help='Raster band to read if multiple bandas available, default band 1.')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')
    parser.add_argument('--kk','--keep-kml','--keep-kml-file', dest='keep_kml_file', action='store_true',
                        help='Do not remove KML and data/resource files after compressing into KMZ file.')
    parser.add_argument('--zero-mask', dest='zero_mask', action='store_true',
                        help='Mask pixels with zero value.')
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('MIN', 'MAX'), type=float,
                        help='Y/value limits for plotting.')
    parser.add_argument('-c', '--cm', '--colormap', dest='cmap_name', default='jet',
                        help='Colormap for plotting (default: %(default)s), such as jet, RdBu, etc.\n'
                             'More details at https://mintpy.readthedocs.io/en/latest/api/colormaps/')
    parser.add_argument('--dpi', dest='fig_dpi', metavar='NUM', type=int, default=600,
                     help='Figure DPI (dots per inch). Default: 600')


    return parser


def cmdLineParse(iargs=None):
    '''
    Command line parser.
    '''

    parser = createParser()
    return parser.parse_args(args = iargs)


min_figsize_single = 6.0       # default min size in inch, for single plot
max_figsize_single = 10.0      # default min size in inch, for single plot
max_figsize_height = 8.0       # max figure size in vertical direction in inch
def auto_figure_size(ds_shape, scale=1.0, disp_cbar=False, disp_slider=False,
                     cbar_ratio=0.25, slider_ratio=0.15, print_msg=True):
    """Get auto figure size based on input data shape
    Adjust if display colobar on the right and/or slider on the bottom

    Parameters: ds_shape          - tuple/list of 2 int for the 2D matrix shape in [length, width]
                scale             - floag, scale the final figure size
                disp_cbar/slider  - bool, plot colorbar on the right / slider on the bottom
                cbar/slider_ratio - float, size ratio of the additional colobar / slider
    Returns:    figsize           - list of 2 float for the figure size in [width, lenght] in inches
    """
    # figure shape
    fig_shape = list(ds_shape)[::-1]
    if disp_cbar:
        fig_shape[0] *= (1 + cbar_ratio)
    if disp_slider:
        fig_shape[1] *= (1 + slider_ratio)

    # get scale to meet the min/max figure size constrain
    fig_scale = min(min_figsize_single / min(fig_shape),
                    max_figsize_single / max(fig_shape),
                    max_figsize_height / fig_shape[1])

    # fig_shape/scale --> fig_size
    fig_size = [i*fig_scale*scale for i in fig_shape]
    if print_msg:
        print(f'figure size : [{fig_size[0]:.2f}, {fig_size[1]:.2f}]')

    return fig_size


def four_corners(atr):
    """Get the 4 corners coordinates from metadata dict in geo-coordinates.
    Parameters: atr - dict
    Returns:    south, north, west, east - float, in degrees or meters
    Examples:   S, N, W, E = ut.four_corners(atr)
                SNWE = ut.four_corners(atr)
    """
    width  = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    west  = float(atr['X_FIRST'])
    north = float(atr['Y_FIRST'])
    south = north + lat_step * length
    east  = west  + lon_step * width
    return south, north, west, east
    

def write_kmz_file(out_file_base, kml_doc, data_files=None, res_files=None, keep_kml_file=False):
    """Write KML and KMZ files.
    Parameters: out_file_base - str, output file name without extension
                kml_doc       - KML.Document() object
                data_files    - list of str, rel path of data files
                res_files     - list of str, rel path of resource files
                keep_kml_file - bool, do not remove KML files after zipping.
    Returns:    kmz_file      - str, zipped KMZ file.
    """
    # default values
    data_files = [] if data_files is None else data_files
    res_files  = [] if res_files  is None else res_files

    work_dir = os.path.dirname(out_file_base)
    kml_file = f'{out_file_base}.kml'
    kmz_file = f'{out_file_base}.kmz'

    # 1. Write KML file
    kml = KML.kml()
    kml.append(kml_doc)

    print('writing '+kml_file)
    with open(kml_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf8'))

    # 2. Copy resource files
    if res_files:
        res_dir = os.path.join(os.path.dirname(mintpy.__file__), "data")
        for fname in res_files:
            src_file = os.path.join(res_dir, os.path.basename(fname))
            shutil.copy2(src_file, work_dir)
            print(f"copy {src_file} to the local directory")

    # 3. Generate KMZ file, by
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(work_dir)

    # 2) zip all data files
    with ZipFile(kmz_file, 'w') as fz:
        for fname in [kml_file] + data_files + res_files:
            fz.write(os.path.relpath(fname))
            if not keep_kml_file:
                os.remove(fname)
                print(f'remove {fname}')

    # 3) go back to the running directory
    os.chdir(run_dir)
    print(f'merged all files to {kmz_file}')

    return kmz_file


def write_kmz_overlay(data, meta, out_file, inps):
    """Generate Google Earth Overlay KMZ file for data in GEO coordinates.
    Parameters: data     - 2D np.array in int/float, data matrix to write
                meta     - dict, containing the following attributes:
                           WIDTH/LENGTH      : required, file size
                           X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                           REF_X/Y           : optional, column/row number of reference pixel
                out_file - string, output file name
                inps     - Namespace, optional, input options for display
    Returns:    kmz_file - string, output KMZ filename
    """

    south, north, west, east = four_corners(meta)

    # 1. Make PNG file - Data
    print('plotting data ...')

    fig_size = None
    # Figure size
    if not fig_size:
        fig_size = auto_figure_size(ds_shape=[north-south, east-west], scale=2.0)
    fig = plt.figure(figsize=fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    # Plot - data matrix
    ax.imshow(data, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.cmap_name,
              aspect='auto', interpolation='nearest')

    width = int(meta['WIDTH'])
    length = int(meta['LENGTH'])
    ax.set_xlim([0, width])
    ax.set_ylim([length, 0])

    out_file_base = os.path.splitext(out_file)[0]
    data_png_file = out_file_base + '.png'
    fig_dpi = inps.fig_dpi
    print(f'writing {data_png_file} with dpi={fig_dpi}')
    plt.savefig(data_png_file, pad_inches=0.0, transparent=True, dpi=fig_dpi)

    # 2. Generate KML file
    kml_doc = KML.Document()

    # Add data png file
    img_name = os.path.splitext(os.path.basename(data_png_file))[0]
    img_overlay = KML.GroundOverlay(
        KML.name(img_name),
        KML.Icon(
            KML.href(os.path.basename(data_png_file))
        ),
        KML.altitudeMode('clampToGround'),
        KML.LatLonBox(
            KML.north(str(north)),
            KML.east(str(east)),
            KML.south(str(south)),
            KML.west(str(west)),
        ),
    )
    kml_doc.append(img_overlay)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[data_png_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


def write_kmz_placemark(data, meta, out_file, geom_file, inps):
    """Generate Google Earth Placemark KMZ file for data in RADAR coordinates.
    Parameters: data      - 2D np.array in int/float, data matrix to write
                meta      - dict, containing the following attributes:
                            WIDTH/LENGTH      : required, file size
                            X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                            REF_X/Y           : optional, column/row number of reference pixel
                geom_file - str, path of the geometry file with latitude/longitude datasets
                out_file  - string, output file name
                inps      - Namespace, optional, input options for display
    Returns:    kmz_file  - string, output KMZ filename
    """

    out_file_base = os.path.splitext(out_file)[0]
    dot_file = 'shaded_dot.png'
    star_file = 'star.png'

    # read latitude / longitude
    lats = readfile.read(geom_file, datasetName='latitude',  box=inps.pix_box)[0]
    lons = readfile.read(geom_file, datasetName='longitude', box=inps.pix_box)[0]

    ## Generate KML file
    kml_doc = KML.Document()

    # 1. colorbar png file
    print('plot and add colorbar as a ScreenOverlay element')
    cbar_file = f'{out_file_base}_cbar.png'
    cbar_overlay = generate_cbar_element(
        cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        unit=inps.disp_unit,
        loc=inps.cbar_loc,
        nbins=inps.cbar_bin_num,
        label=inps.cbar_label)
    kml_doc.append(cbar_overlay)

    # 2. reference point
    xmin = int(meta.get('SUBSET_XMIN', 0))
    ymin = int(meta.get('SUBSET_YMIN', 0))

    if 'REF_Y' in meta.keys():
        print('add reference point as a star icon')
        ry, rx = int(meta['REF_Y']), int(meta['REF_X'])
        rlat = lats[ry, rx]
        rlon = lons[ry, rx]
        ref_point = create_placemark_element(
            lat=rlat,
            lon=rlon,
            row=ry + ymin,
            col=rx + xmin,
            val=0.0,
            icon_file=star_file,
            inps=inps)
        ref_point.name = 'ReferencePoint'
        ref_point.Style.IconStyle.scale = 1.0
        kml_doc.append(ref_point)

        # do not plot reference point as data again
        data[ry, rx] = np.nan

    # 3. data folder for all points
    data_folder = KML.Folder(KML.name("Data"))

    print(f'generating point element with step size of {inps.step} pixels')
    length, width = data.shape
    prog_bar = ptime.progressBar(maxValue=length)
    for y in range(0, length, inps.step):
        for x in range(0, width, inps.step):
            value = data[y, x]
            if not np.isnan(value):
                lat = lats[y, x]
                lon = lons[y, x]

                # create KML icon element
                placemark = create_placemark_element(
                    lat=lat,
                    lon=lon,
                    row=y + ymin,
                    col=x + xmin,
                    val=value,
                    icon_file=dot_file,
                    inps=inps)
                data_folder.append(placemark)

        prog_bar.update(y+1, every=1, suffix=f'row={y+1}/{length}')
    prog_bar.close()
    kml_doc.append(data_folder)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[cbar_file],
        res_files=[dot_file, star_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


def gdal_load(path, band):
    """
    Load image file into array
    using gdal.
    
    Inputs
    ------
    path = path to image file, text.
    band = band to load, int.
    
    Returns
    -------
    loaded file as numpy array
    """
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    arr = ds.GetRasterBand(band).ReadAsArray()
    length, width = ds.RasterYSize, ds.RasterXSize
    # ulx, uly stands for upper left corner, lrx, lry for lower right corner
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (width * xres)
    lry = uly + (length * yres)
    atr = {"WIDTH": width, "LENGTH": length, 
           "X_FIRST":ulx, "X_STEP":xres,"X_UNIT": 'degrees',
           "Y_FIRST":uly,"Y_STEP":yres,"Y_UNIT": 'degrees'}
    ds = None
    
    return arr, atr


def main(iargs=None):
    # Main driver
    inps = cmdLineParse(iargs)
    in_file = os.path.abspath(inps.in_file)
    out_file = os.path.abspath(inps.outfile)
    # read array
    data, atr = gdal_load(in_file, inps.band)
    ref_date = None
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    print(f'data   coverage in y/x: {data_box}')
    # zero mask
    if inps.zero_mask:
        print('masking out pixels with zero value')
        data[data == 0] = np.nan

    # write KMZ file
    write_kmz_overlay(
        data,
        meta=atr,
        out_file=out_file,
        inps=inps) 


if __name__ == '__main__':

    main()

    '''
    Main driver.
    '''
