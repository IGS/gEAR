#!/usr/bin/env python3

"""
Input: Excel file + Config (example below)

Take an Excel file and calculate:
    (1) averages
    (2) standard deviations

Output: Excel file containing original data and the new calculations

Example config file:
    [gene_symbol_column]
    #Index of gene symbol column (Zero-based index)
    column_header=1

    [calculations]
    #currently only averages and standard deviations can be performed
    averages=True
    standard_deviations=True
    p_values=False

    #A group is a list of replicates that average and standard deviation will be
    #  calculated from
    [replicates]
    replicates_1=IN1_CPM,IN2_CPM,IN3_CPM
    replicates_2=IP1_CPM,IP2_CPM,IP3_CPM

    [name_average_columns]
    average_1=IN--1
    average_2=IP--1
"""
import argparse
import configparser
import pandas as pd
import re

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to file containing data file' )
    parser.add_argument('-c', '--config_file', type=str, required=True, help='Path to file containing config file' )

    args = parser.parse_args()


    filename = args.input_file
    config_filename = args.config_file

    config = configparser.ConfigParser()
    config.read(config_filename)

    index_col=int(config['gene_symbol_column']['column_header'])

    columns_list = []
    for group in config['replicates']:
        columns_string = config.get('replicates', group)

        if ',' in columns_string:
            columns = columns_string.split(',')
        else:
            #copy and pasting column names from excel to config will tab delimit names
            columns = re.split(r'\t+', columns_string)
        #print(columns)
        #['IN1_CPM', 'IN2_CPM', 'IN3_CPM']
        #['IP1_CPM', 'IP2_CPM', 'IP3_CPM']
        columns_list.append(columns)

    avg_col_names = []
    for col in config['name_average_columns']:
        name = config.get('name_average_columns', col)
        avg_col_names.append(name)

    #what was_peformed
    calcs_performed = ''

    #Data file is excel
    if filename.lower().endswith( ('xlsx', 'xls') ):
        file_type = 'xlsx'
        df = pd.read_excel(filename, index_col=index_col)


    for group in columns_list:
        #Get a subset of the dataset that only contains 1 group of replicates
        df_group = df.loc[:, group].copy()

        #Get name of averaged column
        col_name = avg_col_names[columns_list.index(group)]

        #Calculate the averages.
        if config.getboolean('calculations', 'averages') is True:
            if '_avg' not in calcs_performed:
                calcs_performed += '_avg'

            df_group_avg = df_group.mean(axis=1).rename(col_name).to_frame()

            #Append averages col to main dataframe
            df = pd.concat([df, df_group_avg], axis=1)

        #Calculate the standard deviations.
        if config.getboolean('calculations', 'standard_deviations') is True:
            if '_sd' not in calcs_performed:
                calcs_performed += '_sd'
            df_group_std = df_group.std(axis=1).rename(col_name + '_sd').to_frame()

            #Append averages col to main dataframe
            df = pd.concat([df, df_group_std], axis=1)




    output_filename = filename.rsplit('.', 1)[0] + calcs_performed + '.' + file_type

    writer = pd.ExcelWriter(output_filename)
    df.to_excel(writer, 'Sheet1')
    writer.save()



if __name__ == '__main__':
    main()
