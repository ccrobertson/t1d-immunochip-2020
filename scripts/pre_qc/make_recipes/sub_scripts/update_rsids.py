import pandas as pd
import sys
#module load anaconda/5.2.0-py3.6

bimfilename=sys.argv[1]
print(bimfilename)

# define functions
flipStrand = {
    "A":"T",
    "T":"A",
    "G":"C",
    "C":"G",
    "0":"0"
    }

mapToXYM = {
    23:"X",
    24:"Y",
    25:"PAR",
    26:"M",
    }

def getChrString(chrNum):
    if chrNum > 22:
        chrStr = "chr" + mapToXYM[chrNum]
    else:
        chrStr = "chr" + str(chrNum)
    return chrStr


def pasteVals(chrpos, a1, a2):
    var = chrpos + ":" + a1 + ":" + a2
    return var


#load dbsnp data
dbmap = pd.read_table("ichip_snp_positions_to_dbsnp.txt", header=None, low_memory=False)
dbmap.columns = ["chr","pos","snp","ref","alt","u1","u2","info"]
dbmap['chrpos'] = "chr" + dbmap['chr'].map(str) + ":" + dbmap['pos'].map(str)
dbmap['var'] = dbmap['chrpos'] + ":" + dbmap['ref'] + ":" + dbmap['alt']
map_dict = dbmap.set_index('var').T.to_dict('dict')

#load ichip map file
bim = pd.read_table(bimfilename, header=None, low_memory=False)
bim.columns = ["chr","snp","map","pos","ref","alt"]
bim['chrStr'] = bim['chr'].apply(getChrString)
bim['chrpos'] = bim['chrStr'].map(str) + ":" + bim['pos'].map(str)
bim['var'] = bim['chrpos'] + ":" + bim['ref'] + ":" + bim['alt']
bim_dict = bim.set_index('var').T.to_dict('dict')

        
count=0
countIn=0
header = ["variant_mega","rsid_mega","variant_dbsnp","rsid_dbsnp","info_dbsnp"]
print("\t".join(str(x) for x in header))
for var in bim_dict:
    count += 1
    #if count > 100:
    #    break
    
    chrpos = bim_dict[var]['chrpos']
    s1_allele1 = bim_dict[var]['ref']
    s1_allele2 = bim_dict[var]['alt']
    s2_allele1 = flipStrand[s1_allele1]
    s2_allele2 = flipStrand[s1_allele2]
    
    #define all possible allele configurations
    var1 = pasteVals(chrpos, s1_allele1, s1_allele2)
    var2 = pasteVals(chrpos, s1_allele2, s1_allele1)
    var3 = pasteVals(chrpos, s2_allele1, s2_allele2)
    var4 = pasteVals(chrpos, s2_allele2, s2_allele1)
    
    if var1 in map_dict:
        countIn += 1
        record = [var, bim_dict[var]['snp'], var1, map_dict[var]['snp'], map_dict[var]['info']]
    elif var2 in map_dict:  
        countIn += 1
        record = [var, bim_dict[var]['snp'], var2, map_dict[var2]['snp'], map_dict[var2]['info']]
    elif var3 in map_dict:
        countIn += 1
        record = [var, bim_dict[var]['snp'], var3, map_dict[var3]['snp'], map_dict[var3]['info']]
    elif var4 in map_dict:
        countIn += 1
        record = [var, bim_dict[var]['snp'], var4, map_dict[var4]['snp'], map_dict[var4]['info']]
    else:
        record = [var, bim_dict[var]['snp'], "NA", "NA", "NA"]
    
    print("\t".join(str(x) for x in record))


#print(count)       
#print(countIn)        

        
        
        
        
        
        
        
        
        
        
        
        
        