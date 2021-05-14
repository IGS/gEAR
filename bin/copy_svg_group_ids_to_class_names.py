#!/usr/bin/env python3

"""
Lance's software doesn't have any export tools which will set the <path> 
element class names properly.  This script fixes this, reading any <g element ids, 
processing them and then assigning them to the <circle child's class.

<g id="D1_x5F_SF18101610">
    <circle class="st103" cx="236.2" cy="366.20001" r="2.8" id="circle1443" style="fill:#f4f2f2;stroke:#666666;stroke-width:0.5;stroke-miterlimit:10" />
</g>

Becomes 

<g id="D1_x5F_SF18101610">
    <circle class="D1_SF18101610" cx="236.2" cy="366.20001" r="2.8" id="circle1443" style="fill:#f4f2f2;stroke:#666666;stroke-width:0.5;stroke-miterlimit:10" />
</g>




"""
import argparse
import os
import re

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def main():
    parser = argparse.ArgumentParser(description='SVG g id -> class copier')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    tree = ET.ElementTree(file=args.input_file)
    root = tree.getroot()

    path_count = process_element(root)

    # write the file back out
    tree.write(args.output_file)

def process_element(e):
    # remove the namespace from the tag
    if len(e.tag.split('}', 1)) > 1:
        e.tag = e.tag.split('}', 1)[1]
    
    # is this as path?  If so, modify it
    if e.tag == 'g':
        # check/set attributes here
        if e.get('id'):
            g_id = e.get('id')
            print("Got ID of: {0}".format(g_id))

            for child in e.getchildren():
                # remove the namespace from the tag
                if len(child.tag.split('}', 1)) > 1:
                    child.tag = child.tag.split('}', 1)[1]

                if child.tag == 'circle':
                    print("modifying child circle")
                    new_id = g_id.replace('x5F_', '')
                    child.set('class', new_id)
                else:
                    print("Skipping tag of {0}".format(child.tag))

    # process all children
    for elem in list(e):
        process_element(elem)

if __name__ == '__main__':
    main()



