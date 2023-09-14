import sys
import csv
from functools import *
from collections import Counter


def identical_col(dbs: list) -> str:
    """Receives a list of targets
    and retuns a string if the majority
    agree/disagree among them"""

    #cleaning list to check conditions
    to_filter = ['----', 'unclassified']
    clean_dbs = [ele for ele in dbs if ele not in to_filter]
    
    if len(set(clean_dbs)) == 1:
        return 'Identical'

    else:
        return 'Not_identical' 


def consensus_col(dbs: list) -> str:
    """Receives a list of targets
    and retuns a string if the majority
    agree/disagree among them"""

    #cleaning list to check conditions
    to_filter = ['----', 'unclassified']
    clean_dbs = [ele for ele in dbs if ele not in to_filter]
    threshold = len(clean_dbs)/2

    #counter dict for each term in the cleaning list
    dict_counter = Counter(clean_dbs)

    #checking conditions #to get the consensus target: get the higher v and print k (target)
    for k,v in dict_counter.items():
        if v > threshold:
            
            return 'Consensus'
            
    return 'No_consensus'



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
    new_first_line = '\t'.join([first_line.strip(),'Consensus_dbs_4',
                                'Consensus_dbs_3',  
                                'Identical_4', 
                                'Identical_3',
                                'number_dbs_targets', 
                                'dbs_target_match_pred', 
                                'number_dbs_target_match_pred', 
                                'number_dbs_target_no_match_pred', 
                                'dbs_agreement', 
                                'number_dbs_agreement', 
                                ])
    # new_first_line = '\t'.join([first_line.strip(),'Consensus_dbs', 'Identical_4', 'Identical_3'])
    
    header = new_first_line.split('\t') #tab

    with open(sys.argv[2], 'w') as f:
        
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        
        for line in file_n:
            samples = line.strip().split('\t')
            list_match_dbs = []
            list_unmatch_dbs = []
            list_dbs = []

            #if add a column, change the indexes here!
            GEO = samples[25].lower()
            NGS = samples[45].lower()
            CA =  samples[69].lower()
            Cdb = samples[102].lower()
            EpiLaP = samples[107].lower()
            cols = [0,0,0,0,0,0,0,0,0,0] #solving else for EpiLaP not available
            # cols = [0,0,0] #solving else for EpiLaP not available
            samples = samples + cols
        
            #creating consensus 4 column
            samples[152] = consensus_col([GEO, NGS, CA, Cdb])

            #creating consensus 3 column
            samples[153] = consensus_col([NGS, CA, Cdb])
           
            #creating identical_4
            samples[154] = identical_col([GEO, NGS, CA, Cdb])

            #creating identical_3 - excluding GEO
            samples[155] = identical_col([ NGS, CA, Cdb])
  
            # writer.writerow(samples)


            #append all dbs infor from metadata to create a complete dbs_agreement and number_dbs_agreement
            list_dbs += [('GEO', GEO)] + [('NGS', NGS)] + [('CA', CA), ('Cdb', Cdb)]
            # print('dbs no filter:',list_dbs)    

            list_dbs = list(filter(lambda dbs: dbs[1]!= '----' and dbs[1]!= 'unclassified', list_dbs))
            # print('dbs filter:',list_dbs)
     

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
            samples[156] = len(list_match_dbs) + len(list_unmatch_dbs) #before 87 index
            
            if len(list_match_dbs) > 0 : #list from comparison with EpiLaP pred
                #dbs with match match pred
                samples[157] = ','.join(list_match_dbs) 

                #number of dbs with match prediction 
                samples[158] = len(list_match_dbs)

                # number of dbs no match with pred 
                samples[159] = len(list_unmatch_dbs)

            #dbs agreement no match pred; number dbs agreement no match pred (or just among them) - using list dbs to fill the entire cols with metadata comparison
            if len(list_dbs) >= 2:
                res_ltup = compare_tuples(list_dbs)
                if len(res_ltup) == 0:
                    samples[160] = 'no match dbs'
                    # print(samples[160])
                    
                else:
                    samples[160] = ",".join(res_ltup)
                    # print(samples[160])

                samples[161] = len(res_ltup)
                # print(samples[161])
                writer.writerow(samples)

            else: #if list len(unmatch) == 0 or 1
                list_dbs_nomatch = [i[0] for i in list_dbs]
                samples[160] = ",".join(list_dbs_nomatch) # we do not have an agreement between databases 
                samples[161] = len(list_dbs)
                
                # print(samples[160])
                # print(samples[161])

            #     #rewriting the file including the new columns
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