#!/usr/bin/env python3

"""
Amiel's software doesn't have any export tools which will set the <path> 
element class names properly.  This script fixes this, reading any id
attributes and setting that also as the class attribute.

Note:  Does <path> element by default 

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
    parser.add_argument('-e', '--svg_element', type=str, required=False, default='path', help='SVG element type to be processed' )
    args = parser.parse_args()

    tree = ET.ElementTree(file=args.input_file)
    root = tree.getroot()

    path_count = process_element(root, args.svg_element)

    # write the file back out
    tree.write(args.output_file)

def process_element(parent, etype):
    # remove the namespace from the tag
    try:
        parent.tag = parent.tag.split('}', 1)[1]
    except IndexError as err:
        pass
    
    # is this as path?  If so, modify it
    if parent.tag == etype:
        # check/set attributes here
        if parent.get('id'):
            parent.set('class', parent.get('id'))

    # process all children
    for elem in list(parent):
        process_element(elem, etype)

if __name__ == '__main__':
    main()



