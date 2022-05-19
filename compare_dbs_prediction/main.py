import sys
import csv
from functools import *


def index_file(file_n):
    '''Receives a opened file
    and prints each column and
    its respective index'''
    first_line = file_n.readline() #getting header
    print(first_line)
    a = first_line.strip().split('\t')
    for i,ele in enumerate(a):
        print(i,ele,'\n')


def compare_tuples(list_unmatch_dbs):
    '''Receives a list of tuples (db,target)
    and returns a list of db with the 
    same target.'''

    res_ltup = []
    for i in range(len(list_unmatch_dbs) -1):
        for j in range(i + 1,len(list_unmatch_dbs)):
            if list_unmatch_dbs[i][1] == list_unmatch_dbs[j][1]:
                res_ltup = res_ltup + [list_unmatch_dbs[i][0], list_unmatch_dbs[j][0]]

    res_ltup = list(set(res_ltup))

    return res_ltup


def get_index_name(file_n):

    first_line = file_n.readline() #getting header
    new_first_line = '\t'.join([first_line.strip(),'number_dbs_targets', 'dbs_target_match_pred', 'number_dbs_target_match_pred', 'dbs_target_no_match_pred', 'dbs_agreement', 'number_dbs_agreement'])
    header = new_first_line.split('\t')

    with open(sys.argv[2], 'w') as f:
        
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        
        for line in file_n:
            samples = line.strip().split('\t')
            list_match_dbs = []
            list_unmatch_dbs = []

            #if add a column, change the indexes here!
            GEO = samples[21].lower()
            NGS = samples[39].lower()
            CA =  samples[50].lower()
            Cdb = samples[83].lower()
            EpiLaP = samples[85].lower()
            cols = [0,0,0,0,0,0] #solving else for EpiLaP not available
            samples = samples + cols
        
            #Filtering Pred values
            if EpiLaP != '----': #all test lines have pred values

                #GEO
                if GEO != '----' and GEO == EpiLaP:
                    list_match_dbs.append('GEO')

                elif GEO != '----' and GEO != EpiLaP:
                    list_unmatch_dbs.append(('GEO', GEO))

                #NGS-QC
                if NGS != '----' and NGS == EpiLaP:
                    list_match_dbs.append('NGS')

                elif NGS != '----' and NGS != EpiLaP:
                    list_unmatch_dbs.append(('NGS', NGS))

                #CA to ignore unclassified samples as well
                if CA != '----' and CA != 'unclassified' and CA == EpiLaP:
                    list_match_dbs.append('CA')

                elif CA != '----' and CA != 'unclassified' and CA != EpiLaP:
                    list_unmatch_dbs.append(('CA', CA))

                #Cdb
                if Cdb != '----' and Cdb == EpiLaP:
                    list_match_dbs.append('Cdb')

                elif Cdb != '----' and Cdb != EpiLaP:
                    list_unmatch_dbs.append(('Cdb', Cdb))

            else:
                #compare databases
                list_unmatch_dbs += [('GEO', GEO)] + [('NGS', NGS)] + [('CA', CA), ('Cdb', Cdb)]
                list_unmatch_dbs = list(filter(lambda dbs: dbs[1]!= '----' and dbs[1]!= 'unclassified', list_unmatch_dbs)) #filtering no target
               
            #Creating new columns
            #number of dbs with targets
            samples[87] = len(list_match_dbs) + len(list_unmatch_dbs)
            
            if len(list_match_dbs) > 0 : #list from comparison with EpiLaP pred
                #dbs with match match pred
                samples[88] = ','.join(list_match_dbs)

                #number of dbs with match prediction 
                samples[89] = len(list_match_dbs)

                # number of dbs no match with pred 
                samples[90] = len(list_unmatch_dbs)

            #dbs agreement no match pred; number dbs agreement no match pred (or just among them)
            if len(list_unmatch_dbs) >= 2:
                res_ltup = compare_tuples(list_unmatch_dbs)
                if len(res_ltup) == 0:
                    samples[91] = 'no match dbs'
                    
                else:
                    samples[91] = ",".join(res_ltup)
                  
                samples[92] = len(res_ltup)
                writer.writerow(samples)

            else: #if list len(unmatch) == 0 or 1
                list_dbs_nomatch = [i[0] for i in list_unmatch_dbs]
                samples[92] = ",".join(list_dbs_nomatch) # we do not have an agreement between databases 
                samples[92] = len(list_unmatch_dbs)

                #rewriting the file including the new columns
                writer.writerow(samples)


def main():

    print('Starting...')
    file_n = open(sys.argv[1], 'r') #file separated by tab (tsv) including Predicted and Max value columns from EpiLaP
    # index_file(file_n)
    # sys.exit()
    output_file = sys.argv[2] #path to output file
    get_index_name(file_n)
    print('File sucessfully saved!!')



if __name__ == "__main__":



    main()