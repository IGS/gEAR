#!/usr/bin/env python3

"""
Amiel's software doesn't have any export tools which will set the <path> 
element class names properly.  This script fixes this, reading any id
attributes and setting that also as the class attribute.

Note:  Only does <path> elements

"""
import argparse
import os
import re

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def main():
    parser = argparse.ArgumentParser(description='SVG id -> class copier')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    tree = ET.ElementTree(file=args.input_file)
    root = tree.getroot()

    path_count = process_element(root)

    # write the file back out
    tree.write(args.output_file)

def process_element(parent):
    # remove the namespace from the tag
    parent.tag = parent.tag.split('}', 1)[1]
    
    # is this as path?  If so, modify it
    if parent.tag == 'path':
        # check/set attributes here
        if parent.get('id'):
            parent.set('class', parent.get('id'))

    # process all children
    for elem in list(parent):
        process_element(elem)

if __name__ == '__main__':
    main()



