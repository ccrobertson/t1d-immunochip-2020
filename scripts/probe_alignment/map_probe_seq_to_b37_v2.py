import inspect
import pandas as pd


### read in manifest (create dictionary -- #https://stackoverflow.com/questions/26716616/convert-a-pandas-dataframe-to-a-dictionary)
manifest = pd.read_csv("Immuno_BeadChip_11419691_Manifest_A.csv", header=7, low_memory=False)
manifest_dict = manifest.set_index('Name').T.to_dict('dict')


### read in blat results
blatA = pd.read_table("manifest_A_alleleA.psl", skiprows=5, header=None, low_memory=False)
blatA.columns = ["matches", "mismatches", "repmatches", "ncount", "qnuminsert", "qbaseinsert","tnuminsert","tbaseinsert","strand","qname","qsize","qstart","qend","tname","tsize","tstart","tend","blockcount","blocksizes","qstarts","tstarts"]
blatB = pd.read_table("manifest_A_alleleB.psl", skiprows=5, header=None, low_memory=False)
blatB.columns = ["matches", "mismatches", "repmatches", "ncount", "qnuminsert", "qbaseinsert","tnuminsert","tbaseinsert","strand","qname","qsize","qstart","qend","tname","tsize","tstart","tend","blockcount","blocksizes","qstarts","tstarts"]


### remove strange contigs
def removeContigs(results):
    bool_crit=["_" not in x for x in results['tname'].values] 
    return results[bool_crit]    
    
blatA_tmp = removeContigs(blatA)
blatB_tmp = removeContigs(blatB)


### restrict results to perfect matches
def choosePerfectMatch(results):
    perfect_matches = results[(results['mismatches']==0) & (results['qstart']==0) & (results['qend']==results['qsize'])]
    return perfect_matches
    
blatA_perfect = choosePerfectMatch(blatA_tmp)
blatB_perfect = choosePerfectMatch(blatB_tmp)


### define PAR1 and PAR2 
#https://www.ncbi.nlm.nih.gov/grc/human
#PAR1 = X:60001-2699520 & Y:10001-2649520
#PAR2 = X:154931044-155260560 & Y:59034050-59363566
PAR1 = {'chrX': {'start':60001, 'end':2699520},
        'chrY': {'start':10001, 'end':2649520}}            
PAR2 = {'chrX': {'start':154931044, 'end':155260560},
        'chrY': {'start':59034050, 'end':59363566}}

### define function to check if variant is in PAR
def isInPar(chr,pos, PAR): 
    if (chr in PAR):
        if pos >= PAR[chr]['start'] and pos <= PAR[chr]['end']:
            return True
        else:
            return False
    else:
        return False
            
            
### for multimapping probes from PAR keep X match only (drop Y matches)
def removeYPARMatches(results):
    multiple_matches = results['qname'].value_counts().index[results['qname'].value_counts()>1]
    results_mm = results[results['qname'].isin(multiple_matches)]
    xy_list = []
    for var in multiple_matches:
        df = results_mm[results_mm['qname']==var]    
        no_matches = df.shape[0]
        no_chrX = sum(df['tname'].values=="chrX")
        no_chrY = sum(df['tname'].values=="chrY")
        if (no_chrX==1) & (no_chrY==1) & (no_matches==2):
            if isInPar(df['tname'].values[0], df['tstart'].values[0], PAR1) and isInPar(df['tname'].values[1], df['tstart'].values[1], PAR1):
                xy_list.append(var)
            elif isInPar(df['tname'].values[0], df['tstart'].values[0], PAR2) and isInPar(df['tname'].values[1], df['tstart'].values[1], PAR2):
                xy_list.append(var)             
    results_dropped_ypar = results.loc[(~((results['qname'].isin(xy_list)) & (results['tname']=="chrY")))]
    results_dropped_ypar.loc[results_dropped_ypar['qname'].isin(xy_list), 'tname'] = "chrXY" 
    return results_dropped_ypar


blatA_perfect2 = removeYPARMatches(blatA_perfect)
blatB_perfect2 = removeYPARMatches(blatB_perfect)
                                                                
### remove any remaining variants with duplicate matches (can't trust these!)
def removeDupMatches(results):
    multiple_matches = results['qname'].value_counts().index[results['qname'].value_counts()>1]
    print("Number of probes excluded due to more than 1 perfect match:", len(multiple_matches))
    return results[~results['qname'].isin(multiple_matches)]

blatA_final = removeDupMatches(blatA_perfect2)
blatB_final = removeDupMatches(blatB_perfect2)



### create blat results dictionary
blatA_dict = blatA_final.set_index('qname').T.to_dict('dict')
blatB_dict = blatB_final.set_index('qname').T.to_dict('dict')

### read in quinlin map
qmap = pd.read_table("quinlin_map.txt", header=0, low_memory=False)
qmap_dict = qmap.set_index('name').T.to_dict('dict')


### function to extract snps position from blat output
def extract_snppos(varname, blat_dict):
    if varname in blat_dict:
        strand = blat_dict[varname]['strand']
        chr = blat_dict[varname]['tname']
        if strand=="+":
            start = blat_dict[varname]['tend']
            end = blat_dict[varname]['tend'] + 1
        elif strand=="-":
            start = blat_dict[varname]['tstart'] - 1
            end = blat_dict[varname]['tstart']
        else:
            start = "invalid_strand"
            end =  "invalid_strand"
    else:
        strand = 'NA'
        chr = 'NA'
        start = 'NA'
        end = 'NA'
    return [strand, chr, start, end]     


### define start and end position for each variant
for var in manifest_dict: 
        
    #get mapped snp position from allele A probeseq    
    resultsA = extract_snppos(var,blatA_dict)        
    
    #get mapped snp position from allele B probeseq
    resultsB = extract_snppos(var,blatB_dict)
        
    #get quinlin snp position    
    if var in qmap_dict:
        quin_start = qmap_dict[var]['hg19_start']
        quin_end =  qmap_dict[var]['hg19_end']
        quin_chr = qmap_dict[var]['hg19_chrom']
    else:
        quin_start = 'NA'
        quin_end = 'NA'
        quin_chr = 'NA'                    
    
    #get top allele
    snp =  manifest_dict[var]['SNP']
    
    #return results    
    record = [var, snp, resultsA, resultsB, quin_chr, quin_start, quin_end]
    print("\t".join(str(x) for x in record))
    
        
