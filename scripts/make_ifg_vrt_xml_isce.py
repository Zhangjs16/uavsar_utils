#!/usr/bin/env python3

############################################################
# Script to generate UAVSAR ifgs xmls and vrts using isce. #
# Author: Talib Oliver, 2021 (Caltech/JPL)                 #
############################################################

import os
import glob
import argparse
import isce
import isceobj


def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser(description='Script to generate UAVSAR ifgs xmls and vrts using isce: int | amp | cor | unw | conncomp')
    parser.add_argument('-i', '--input', dest='ifgdir', type=str, required=True,
            help='Ifgram directory.')
    parser.add_argument('-s', '--samples', dest='samples', type=int, required=True,
            help='Image sample number. Found as Slant Range Data Range Samples in UAVSAR .ann file')
    return parser


def cmdLineParse(iargs=None):
    '''
    Command line parser.
    '''

    parser = createParser()
    return parser.parse_args(args = iargs)


def run_ifg(int_file, samples):
    outInt = isceobj.Image.createIntImage()
    outInt.setFilename(int_file)
    outInt.setWidth(samples)
    outInt.setAccessMode('read')
    outInt.renderHdr()
    outInt.renderVRT()

def run_unw(unw_file, samples):
    outCor = isceobj.Image.createImage()
    outCor.setFilename(unw_file)
    outCor.setWidth(samples)
    outCor.setAccessMode('read')
    outCor.setDataType('FLOAT')
    outCor.renderHdr()
    outCor.renderVRT()

def run_amp(amp_file, samples):
    outAmp = isceobj.Image.createAmpImage()
    outAmp.setFilename(amp_file)
    outAmp.setWidth(samples)
    outAmp.setAccessMode('read')
    outAmp.renderHdr()
    outAmp.renderVRT()

def run_cor(cor_file, samples):
    outCor = isceobj.Image.createImage()
    outCor.setFilename(cor_file)
    outCor.setWidth(samples)
    outCor.setAccessMode('read')
    outCor.setDataType('FLOAT')
    outCor.renderHdr()
    outCor.renderVRT()

def run_cc(cc_file, samples):
    outCc = isceobj.Image.createImage()
    outCc.setFilename(cc_file)
    outCc.setWidth(samples)
    outCc.setAccessMode('read')
    outCc.setDataType('BYTE') # BYTE | FLOAT | UINT
    outCc.renderHdr()
    outCc.renderVRT()


def main(iargs=None):

    inps = cmdLineParse(iargs)  
    ifgs = glob.glob(inps.ifgdir+'/*')
    print("Input files:", ifgs)
    for ifg_file in ifgs:
        img = os.path.basename(ifg_file)
        ftype = img.split('.')[1]
        if ftype == 'int': 
            run_ifg(img, inps.samples)
        elif ftype == 'amp':
            run_amp(img, inps.samples)
        elif ftype == 'cor':
            run_cor(img, inps.samples)
        elif ftype == 'unw':
            run_unw(img, inps.samples)
        elif ftype == 'conncomp':
            run_cc(img, inps.samples)
    print('Finished vrts and xmls')

if __name__ == '__main__':

    main() 
 
    '''
    Main driver.
    '''
