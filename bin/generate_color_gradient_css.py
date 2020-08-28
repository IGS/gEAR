#!/usr/bin/env python3

"""

Creates a css file with 200 color level classes, named like:

.level1_10
.level1_15
.level1_20

Where the range is 0.00 -> 10.00 in 0.05 increments.

http://www.perbang.dk/rgbgradient/

"""
import os
import re
import sys
from colour import Color

def main():
    purple_steps = [
        '#FFFFFF',
        '#EAE7EB', 
        '#D5CFD8',
        '#C0B7C5',
        '#AB9FB1',
        '#97879E', 
        '#826F8B', 
        '#6D5777', 
        '#583F64', 
        '#432751', 
        '#2F103E'
    ]

    dark_purple_hsv_steps = [
        '#FFFFFF',
        '#E9D8DA',
        '#D4B4BD',
        '#BF94A5',
        '#AA7792',
        '#955D82',
        '#7F4674',
        '#6A3367',
        '#512255',
        '#371540',
        '#1F0A2A'
    ]

    # red -> white -> green
    reverse_italian_steps = [
        '#A30004',
        '#B53336',
        '#C76668',
        '#DA999A',
        '#ECCCCC',
        '#FFFFFF',
        '#CCF1CC',
        '#99E49A',
        '#66D768',
        '#33CA36',
        '#00BD04',
    ]

    # blue -> white -> red
    netherlands_steps = [
        '#0121BC',
        '#334DC9',
        '#6679D6',
        '#99A6E4',
        '#CCD2F1',
        '#FFFFFF',
        '#ECCDCC',
        '#D99B99',
        '#C76967',
        '#B43734',
        '#A20602'
    ]

    netherlands_linear = [
        '#0048FA',
        '#3F75FB',
        '#7FA3FC',
        '#BFD1FD',
        '#FFFFFF',
        '#FFC0BF',
        '#FF817F',
        '#FF423F',
        '#FF0300'
    ]

    # white needs to be at 6_65
    netherlands_neg4_to_2_steps = [
        '#0121BC',
        '#2B46C7',
        '#556BD2',
        '#8090DD',
        '#AAB5E8',
        '#D4DAF3',
        '#FFFFFF',
        '#E7C0BF',
        '#D08280',
        '#B94441',
        '#A20602'
    ]

    # dataset 6
    netherlands_neg3_to_6_steps = [
        '#0121BC',
        '#556BD2',
        '#AAB5E8',
        '#FFFFFF',
        '#EFD5D4',
        '#E0ACAA',
        '#D08280',
        '#C15956',
        '#B12F2C',
        '#A20602',
        '#6B0300'
    ]

    blue_to_red_steps = [
        '#0121BC',
        '#191DA9',
        '#311A97',
        '#4A1785',
        '#621373',
        '#7B1061',
        '#930D4F',
        '#AB093D',
        '#C4062B',
        '#DC0319',
        '#F50007'
    ]

    custom_ronna_linear = [
        '#F9F9F9',
        '#EFE2E4',
        '#E5CDD1',
        '#DBB8C2',
        '#D1A4B5',
        '#C791AB',
        '#BD80A2',
        '#B36F9C',
        '#A86097',
        '#9E5193',
        '#944490',
        '#87388A',
        '#752D80',
        '#632376',
        '#511A6C',
        '#401362'
    ]

    steps = netherlands_linear

#    for level in range(0, 10):
#        step = 0
        
#        while step < 100:
            #print("svg#dataset6 .level{0}_{1} {{ fill: {2} }}".format(level, step, palette[idx]))
            # "7_95" : "rgb(255,95,95)",
            #print("            \"{0}_{1}\" : \"{2}\",".format(level, step, palette[idx]))
#            step += 5
#            idx += 1

    # red to white
    current_w = 0
    for i in range(0, 127):
        print("'rgb(255,{0},{0})',".format(current_w + (i * 2)), end="")

    print("'rgb(255,255,255)',", end="")

    current_w = 255
    for i in range(0, 127):
        print("'rgb({0},{0},255)',".format(current_w - (i * 2)), end="")

    sys.exit()

    for i in range(0, len(steps) - 1):
        minc = Color(steps[i])
        maxc = Color(steps[i+1])
        rangec = list(minc.range_to(maxc, 29))

        for x in range(0, len(rangec)):
            #print(".level{0}_{1} {{ fill: {2} }}".format(i, x * 5, rangec[x]))
            print("'{0}',".format(rangec[x]), end="")

        #print(".level{0}_100 {{ fill: {1} }}".format(i, rangec[-1]))
        print("'{0}',".format(rangec[-1]), end="")

    #print(".level10_0 {{ fill: {0} }}".format(steps[-1]))
    print("'{0}',".format(steps[-1]), end="")

if __name__ == '__main__':
    main()







