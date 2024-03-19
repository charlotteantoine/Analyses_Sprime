##Check overlap between pop for Neanderthal segment to estimate the number admixture events

import csv

# Definir les populations de pop1 et pop2 pour la matrice
pop1_list = ['CHB', 'GBR', 'AKI', 'AKZ', 'BOU', 'HKS', 'KAZ', 'KIB', 'SHO', 'TAB', 'TJE', 'TLG', 'TUB', 'TUR', 'UZB', 'CBR', 'CJO', 'CKC', 'CMR', 'CPN', 'CRE', 'CSA', 'CST', 'CTA', 'LPR', 'LKM', 'LMT', 'LPO', 'LTA', 'KEK', 'KEM', 'MNG']
pop2_list = ['CHB', 'GBR', 'AKI', 'AKZ', 'BOU', 'HKS', 'KAZ', 'KIB', 'SHO', 'TAB', 'TJE', 'TLG', 'TUB', 'TUR', 'UZB', 'CBR', 'CJO', 'CKC', 'CMR', 'CPN', 'CRE', 'CSA', 'CST', 'CTA', 'LPR', 'LKM', 'LMT', 'LPO', 'LTA', 'KEK', 'KEM', 'MNG']

# Create une matrice de confusion
confusion_matrix = {}
for pop1 in pop1_list:
    confusion_matrix[pop1] = {}
    for pop2 in pop2_list:
        confusion_matrix[pop1][pop2] = 0

# Add fragment to dico for pop 1
for pop1 in pop1_list:
    CHB = open(f"/home/charlotte/Data/Sprime/step5/{pop1}/{pop1}.YRI.150000.summary.txt", 'r')
    next(CHB)
    length_vindija_CHB = 0
    fragment_vindija_CHB = {}
    key_vindija_CHB = 'chr0_0'
    fragment_vindija_CHB[key_vindija_CHB] = {'start': 0, 'end': 0}
    
    for line in CHB:
        line = line.split()
        chrom, seg, start, end, vindija, denisova = line[0], line[1], line[2], line[3], line[4], line[5]

        if vindija == 'NA':
            vindija = float(0)
        if denisova == "NA":
            denisova = float(0)
        vindija = float(vindija)
        denisova = float(denisova)

        if vindija > 0.6 and denisova < 0.4:
            
            if chrom == key_vindija_CHB.split('_')[0] and (int(start) > fragment_vindija_CHB[key_vindija_CHB]['start']) and int(end) > int(fragment_vindija_CHB[key_vindija_CHB]['end']) and int(start) < fragment_vindija_CHB[key_vindija_CHB]['end']:
                length_vindija_CHB += int(end) - int(fragment_vindija_CHB[key_vindija_CHB]['end'])
                key_vindija_CHB = f'{chrom}_{seg}'
                fragment_vindija_CHB[key_vindija_CHB] = {'start': int(start), 'end': int(end)}
            
            elif chrom == key_vindija_CHB.split('_')[0] and int(start) > fragment_vindija_CHB[key_vindija_CHB]['start'] and int(end) < fragment_vindija_CHB[key_vindija_CHB]['end']:
                length_vindija_CHB += 0
                key_vindija_CHB = f'{chrom}_{seg}'
                fragment_vindija_CHB[key_vindija_CHB] = {'start': int(start), 'end': int(end)}
            
            elif chrom == key_vindija_CHB.split('_')[0] and int(start) > fragment_vindija_CHB[key_vindija_CHB]['start'] and int(end) == fragment_vindija_CHB[key_vindija_CHB]['end']:
                length_vindija_CHB += 0
                key_vindija_CHB = f'{chrom}_{seg}'
                fragment_vindija_CHB[key_vindija_CHB] = {'start': int(start), 'end': int(end)}
            
            elif chrom < key_vindija_CHB.split('_')[0] :
                length_vindija_CHB += int(end) - int(start)
                key_vindija_CHB = f'{chrom}_{seg}'
                fragment_vindija_CHB[key_vindija_CHB] = {'start': int(start), 'end': int(end)}
            
            elif chrom == key_vindija_CHB.split('_')[0] and (int(start) > fragment_vindija_CHB[key_vindija_CHB]['start']) and int(end) > int(fragment_vindija_CHB[key_vindija_CHB]['end']) and int(start) > fragment_vindija_CHB[key_vindija_CHB]['end']:
                length_vindija_CHB += int(end) - int(start)
                key_vindija_CHB = f'{chrom}_{seg}'
                fragment_vindija_CHB[key_vindija_CHB] = {'start': int(start), 'end': int(end)}
    CHB.close()


##Add fragment to dico for pop2
    for pop2 in pop2_list:
        length_vindija_GBR = 0
        fragment_vindija_GBR = {}
        key_vindija_GBR = 'chr0_0'
        fragment_vindija_GBR[key_vindija_GBR] = {'start': 0, 'end': 0}
        
        GBR = open(f"/home/charlotte/Data/Sprime/step5/{pop2}/{pop2}.YRI.150000.summary.txt", 'r')
        next(GBR)
        
        for line in GBR:
            line = line.split()
            chrom, seg, start, end, vindija, denisova = line[0], line[1], line[2], line[3], line[4], line[5]
            
            if vindija == 'NA':
                vindija = float(0)
            if denisova == "NA":
                denisova = float(0)
            vindija = float(vindija)
            denisova = float(denisova)
            
            if vindija > 0.6 and denisova < 0.4:
                
                if (chrom == key_vindija_GBR.split('_')[0]) and (int(start) > fragment_vindija_GBR[key_vindija_GBR]['start']) and int(end) > int(fragment_vindija_GBR[key_vindija_GBR]['end']) and int(start) < fragment_vindija_GBR[key_vindija_GBR]['end']:
                    length_vindija_GBR += int(end) - int(fragment_vindija_GBR[key_vindija_GBR]['end'])
                    key_vindija_GBR = f'{chrom}_{seg}'
                    fragment_vindija_GBR[key_vindija_GBR] = {'start': int(start), 'end': int(end)}
                
                elif chrom == key_vindija_GBR.split('_')[0] and int(start) > fragment_vindija_GBR[key_vindija_GBR]['start'] and int(end) < fragment_vindija_GBR[key_vindija_GBR]['end']:
                    length_vindija_GBR += 0
                    key_vindija_GBR = f'{chrom}_{seg}'
                    fragment_vindija_GBR[key_vindija_GBR] = {'start': int(start), 'end': int(end)}
                
                elif chrom == key_vindija_GBR.split('_')[0] and int(start) > fragment_vindija_GBR[key_vindija_GBR]['start'] and int(end) == fragment_vindija_GBR[key_vindija_GBR]['end']:
                    length_vindija_GBR += 0
                    key_vindija_GBR = f'{chrom}_{seg}'
                    fragment_vindija_GBR[key_vindija_GBR] = {'start': int(start), 'end': int(end)}
                
                elif chrom < key_vindija_GBR.split('_')[0] :
                    length_vindija_GBR += int(end) - int(start)
                    key_vindija_GBR = f'{chrom}_{seg}'
                    fragment_vindija_GBR[key_vindija_GBR] = {'start': int(start), 'end': int(end)}
                
                elif chrom == key_vindija_GBR.split('_')[0] and (int(start) > fragment_vindija_GBR[key_vindija_GBR]['start']) and int(end) > int(fragment_vindija_GBR[key_vindija_GBR]['end']) and int(start) > fragment_vindija_GBR[key_vindija_GBR]['end']: 
                    length_vindija_GBR += int(end) - int(start)
                    key_vindija_GBR = f'{chrom}_{seg}'
                    fragment_vindija_GBR[key_vindija_GBR] = {'start': int(start), 'end': int(end)}
        GBR.close()

##Compute overlap between pop1 and pop2
        overlap = 0
        for key_vindija_CHB in fragment_vindija_CHB:
            frag_CHB = fragment_vindija_CHB[key_vindija_CHB]
            for key_vindija_GBR in fragment_vindija_GBR:
                frag_GBR = fragment_vindija_GBR[key_vindija_GBR]
                if key_vindija_GBR.split('_')[0] == key_vindija_CHB.split('_')[0]:
                    
                    if frag_GBR['start'] < frag_CHB['start'] and frag_GBR['end'] < frag_CHB['end'] and frag_CHB['start'] < frag_GBR['end']:
                        overlap += int(frag_GBR['end']) - int(frag_CHB['start'])
                    
                    elif frag_GBR['start'] > frag_CHB['start'] and frag_CHB['end'] > frag_GBR['start'] and frag_GBR['end'] > frag_CHB['end']:
                        overlap += int(frag_CHB['end']) - int(frag_GBR['start'])
                    
                    elif frag_GBR['start'] > frag_CHB['start'] and frag_GBR['end'] < frag_CHB['end']:
                        overlap += int(frag_GBR['end']) - int(frag_GBR['start'])
                    
                    elif frag_GBR['start'] < frag_CHB['start'] and frag_GBR['end'] > frag_CHB['end']:
                        overlap += int(frag_CHB['end']) - int(frag_CHB['start'])
                    
                    elif frag_GBR['start'] < frag_CHB['start'] and frag_GBR["end"] == frag_CHB['end']:
                        overlap += int(frag_CHB['end']) - int(frag_CHB['start'])
                    
                    elif frag_GBR['start'] > frag_CHB['start'] and frag_GBR["end"] == frag_CHB['end']:
                        overlap += int(frag_GBR['end']) - int(frag_GBR['start'])
                    
                    elif frag_GBR['start'] == frag_CHB['start'] and frag_GBR['end'] == frag_CHB['end']:
                        overlap += int(frag_GBR['end']) - int(frag_GBR['start'])
                    
                    elif frag_GBR['start'] == frag_CHB['start'] and frag_GBR['end'] < frag_CHB['end']:
                        overlap += int(frag_GBR['end']) - int(frag_GBR['start'])
                    
                    elif frag_GBR['start'] == frag_CHB['start'] and frag_GBR['end'] > frag_CHB['end']:
                        overlap += int(frag_CHB['end']) - int(frag_CHB['start'])
        
        confusion_matrix[pop1][pop2] = round((overlap*2)/(length_vindija_GBR + length_vindija_CHB) *100, 2)
        
with open('/home/charlotte/Data/Sprime/step5/vindija_overlap_between_pop.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write the header
        csv_writer.writerow([''] + pop2_list)
        
        # Write the data
        for pop1 in pop1_list:
            row = [pop1]
            for pop2 in pop2_list:
                row.append(str(confusion_matrix[pop1][pop2]))
            csv_writer.writerow(row)
