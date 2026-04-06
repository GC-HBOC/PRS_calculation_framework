import sys
import argparse
import os

# implementation based exclusively on Python standard library


parser = argparse.ArgumentParser(epilog="See https://github.com/GC-HBOC/PRS_calculation_framework/ for further documentation.")
parser.add_argument("-o", "--output", help="Specification of CanRisk-compatible output VCF file. Default: Input VCF file name with suffix *.canrisk.vcf")
parser.add_argument("-a", "--anc", type=str, choices=['AFR', 'EAS', 'EUR', 'SAS'], help="Specification of presumed ancestry (AFR, EAS, EUR, or SAS). If set, ancestry check is omitted.")
parser.add_argument("-ap", "--anc_prefix", action = 'store_true', help="If set, a prefix specifying predicted or pre-defined ancestry will be added to the output file name.")
parser.add_argument("-dec", "--dec_places", type=int, choices=range(0,10), default=3, help="Maximum number of decimal places in terminal output (default: 3).")
parser.add_argument("-d", "--min_depth", type=int, default=10, help="Minimum read depth required for genotyping (default: 10).")
parser.add_argument("prs_template_file", help="PRS template file in TSV format.")
parser.add_argument("VCF_file", help="Sample VCF file, generated via BCFtools. File ending has to be *.vcf")
args = parser.parse_args()

TEMPLATE = args.prs_template_file
INPUT_VCF = args.VCF_file
if INPUT_VCF[-4:] != ".vcf": sys.exit("VCF input file name must end with .vcf\n")
OFNAME = INPUT_VCF[:-4] + '.canrisk.vcf' if args.output == None else args.output
MINDP = args.min_depth
ANC = args.anc
VERSION = '1.0'

CONVERT_DICT = dict()
AFS = dict()
REVERSE_DICT = dict()
PC1, PC2, PC3 = [], [], []
REF_VARS = []
WEIGHTS = []
AFR_MEAN, AFR_SD = None, None
EAS_MEAN, EAS_SD = None, None
EUR_MEAN, EUR_SD = None, None
SAS_MEAN, SAS_SD = None, None

### PARSE TEMPLATE FILE

with open(TEMPLATE) as infile:
    for line in infile:
        if line.startswith('#CHROM'):
            #CHROM	POS	REF	ALT	EFF	AF_AFR	AF_EAS	AF_EUR	AF_SAS	AF_MEAN	PC1	PC2	PC3	VCF_CHROM	VCF_POS	VCF_REF	VCF_ALT	REVERSE
            # get column indices, first 4 columns have to be PRS locus, remaining columns can be swapped
            ll = line.rstrip().split('\t')
            try:
                EFF_IND = ll.index('EFF')
                AF_AFR_IND, AF_EAS_IND, AF_EUR_IND, AF_SAS_IND, AF_MEAN_IND =  ll.index('AF_AFR'), ll.index('AF_EAS'), ll.index('AF_EUR'), ll.index('AF_SAS'), ll.index('AF_MEAN')
                PC1_IND, PC2_IND, PC3_IND = ll.index('PC1'), ll.index('PC2'), ll.index('PC3')
                VCF_CHROM_IND, VCF_POS_IND, VCF_REF_IND, VCF_ALT_IND, REVERSE_IND = ll.index('VCF_CHROM'), ll.index('VCF_POS'), ll.index('VCF_REF'), ll.index('VCF_ALT'), ll.index('REVERSE')
            except:
                sys.exit('Can not parse template file, is it in valid format??\n')
        
        elif line.startswith('#'):
            # parse means and standard deviations (sds)
            if line.startswith('#AFR_MEAN='):
                try:
                    AFR_MEAN = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse AFR_MEAN.\n')
            elif line.startswith('#AFR_SD='):
                try:
                    AFR_SD = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse AFR_SD.\n')
            if line.startswith('#EAS_MEAN='):
                try:
                    EAS_MEAN = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse EAS_MEAN.\n')
            elif line.startswith('#EAS_SD='):
                try:
                    EAS_SD = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse EAS_SD.\n')
            elif line.startswith('#EUR_MEAN='):
                try:
                    EUR_MEAN = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse EUR_MEAN.\n')
            elif line.startswith('#EUR_SD='):
                try:
                    EUR_SD = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse EUR_SD.\n')
            elif line.startswith('#SAS_MEAN='):
                try:
                    SAS_MEAN = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse SAS_MEAN.\n')
            elif line.startswith('#SAS_SD='):
                try:
                    SAS_SD = float(line.rstrip().split('=')[1])
                except:
                    sys.stderr.write('Unable to parse SAS_SD.\n')
            
        else:
            ll = line.rstrip().split()
            REF_VAR = (ll[0], ll[1], ll[2], ll[3])
            ### store every variant inclusively ALT allele with REF[0] and ALT = . 
            CONVERT_DICT[(ll[VCF_CHROM_IND], ll[VCF_POS_IND], ll[VCF_REF_IND])] = (REF_VAR, (ll[VCF_ALT_IND], '.')) # check ALT  alleles extra
            if len(ll[VCF_REF_IND]) > 1:
                CONVERT_DICT[(ll[VCF_CHROM_IND], ll[VCF_POS_IND], ll[VCF_REF_IND][0])] = (REF_VAR, ('.')) # check ALT  alleles extra
                  
            REF_VARS.append(REF_VAR)
            REVERSE_DICT[REF_VAR] = True if ll[REVERSE_IND] == "1" else False
            AFS[REF_VAR] = (float(ll[AF_AFR_IND]),  float(ll[AF_EAS_IND]), float(ll[AF_EUR_IND]), float(ll[AF_SAS_IND]), float(ll[AF_MEAN_IND]) ) # AF_AFR	AF_EAS	AF_EUR	AF_SAS	AF_MEAN
            # TODO catch typecast errors
            PC1.append(float(ll[PC1_IND]))
            PC2.append(float(ll[PC2_IND]))
            PC3.append(float(ll[PC3_IND]))
            WEIGHTS.append(float(ll[EFF_IND]))

## PARSE VCF            
FOUND = [0] * len(REF_VARS)
GT = [0] * len(REF_VARS)
DP = [0] * len(REF_VARS) ### extra DP list in case input VCFs contains already dosages

with open(INPUT_VCF) as infile:
    for line in infile:
        if not line.startswith('#'):
            ll = line.rstrip().split('\t')
            VAR = (ll[0], ll[1], ll[3])
            gt = ll[9].split(':')[0] # should be '0/0', '0/1', or '1/1' 
            if VAR in CONVERT_DICT and ll[4] in CONVERT_DICT[VAR][1]:
                vind = REF_VARS.index(CONVERT_DICT[VAR][0])
                FOUND[vind] = 1
                GT[vind] = gt
                try:
                    dp = int([_ for _ in ll[7].split(';') if _.startswith('DP=')][0].split('=')[1])
                    DP[vind] = dp
                except:
                    pass
        elif line.startswith('#CHROM'):
            ll = line.rstrip().split('\t')
            SAMPLE = ll[9]

sys.stderr.write('### Sample ' + SAMPLE + '\n' )


sys.stderr.write('Found ' + str(sum(FOUND)) + ' of ' + str(len(FOUND)) +' PRS variants!\n')
if sum(FOUND) < len(FOUND): sys.stderr.write('Could not find:\n' )
for i in range(len(FOUND)):
    if FOUND[i] == 0:
        REF_VAR = REF_VARS[i]
        sys.stderr.write('\t'.join([REF_VAR[0], REF_VAR[1], REF_VAR[2], REF_VAR[3]]) + '\n') 


ANCS = ["AFR", "EAS", "EUR", "SAS"]
if ANC == None:
    #### ANC CHECK
    ANC_GT = []
    c = 0
    for i in range(len(FOUND)):
        if FOUND[i] and GT[i] in ['0/0', '0/1', '1/1'] and DP[i] >= MINDP:
            if not REVERSE_DICT[REF_VARS[i]]:
                if GT[i] == '0/0':
                    ANC_GT.append(0)
                elif GT[i] == '0/1':
                    ANC_GT.append(0.5)
                else:
                     ANC_GT.append(1)
            else:
                if GT[i] == '0/0':
                    ANC_GT.append(1)
                elif GT[i] == '0/1':
                    ANC_GT.append(0.5)
                else:
                     ANC_GT.append(0)
        else:
            c +=1
            ANC_GT.append(float(AFS[REF_VARS[i]][-1]))
            
    sys.stderr.write('### ANCESTRY CHECK\n' )
    sys.stderr.write('Using genotypes from ' + str(len(REF_VARS)-c) + ' out of ' + str(len(REF_VARS)) + ' variants\n' )

    # AFR
    af = [AFS[_][0] for _ in REF_VARS]
    X, Y, Z = 0, 0, 0
    for pc1, pc2, pc3, gt in zip(PC1, PC2, PC3, af): X, Y, Z = X + pc1*gt, Y + pc2*gt, Z + pc3*gt
    AFR_X, AFR_Y, AFR_Z = X, Y, Z
    sys.stderr.write('AFR data point is ' + str(round(X, args.dec_places)) + ' ' + str(round(Y, args.dec_places)) + ' ' + str(round(Z, args.dec_places)) + '\n' )
    # EAS
    af = [AFS[_][1] for _ in REF_VARS]
    X, Y, Z = 0, 0, 0
    for pc1, pc2, pc3, gt in zip(PC1, PC2, PC3, af): X, Y, Z = X + pc1*gt, Y + pc2*gt, Z + pc3*gt
    EAS_X, EAS_Y, EAS_Z = X, Y, Z
    sys.stderr.write('EAS data point is ' + str(round(X, args.dec_places)) + ' ' + str(round(Y, args.dec_places)) + ' ' + str(round(Z, args.dec_places)) + '\n' )
    # EUR
    af = [AFS[_][2] for _ in REF_VARS]
    X, Y, Z = 0, 0, 0
    for pc1, pc2, pc3, gt in zip(PC1, PC2, PC3, af): X, Y, Z = X + pc1*gt, Y + pc2*gt, Z + pc3*gt
    EUR_X, EUR_Y, EUR_Z = X, Y, Z
    sys.stderr.write('EUR data point is ' + str(round(X, args.dec_places)) + ' ' + str(round(Y, args.dec_places)) + ' ' + str(round(Z, args.dec_places)) + '\n' )
    # SAS
    af = [AFS[_][3] for _ in REF_VARS]
    X, Y, Z = 0, 0, 0
    for pc1, pc2, pc3, gt in zip(PC1, PC2, PC3, af): X, Y, Z = X + pc1*gt, Y + pc2*gt, Z + pc3*gt
    SAS_X, SAS_Y, SAS_Z = X, Y, Z
    sys.stderr.write('SAS data point is ' + str(round(X, args.dec_places)) + ' ' + str(round(Y, args.dec_places)) + ' ' + str(round(Z, args.dec_places)) + '\n' )


    X = 0
    for pc, gt in zip(PC1, ANC_GT): X += pc*gt
    Y = 0
    for pc, gt in zip(PC2, ANC_GT): Y += pc*gt
    Z = 0
    for pc, gt in zip(PC3, ANC_GT): Z += pc*gt
    sys.stderr.write('Sample data point is ' + str(round(X, args.dec_places)) + ' ' + str(round(Y, args.dec_places)) + ' ' + str(round(Z, args.dec_places)) + '\n' )
    edafr = ( (X-AFR_X)**2 + (Y - AFR_Y)**2 + (Z - AFR_Z)**2 )**(1/2) 
    sys.stderr.write('Euclidean distance to AFR data point is ' + str(round(edafr, args.dec_places)) + '\n' )
    edeas = ( (X-EAS_X)**2 + (Y - EAS_Y)**2 + (Z - EAS_Z)**2 )**(1/2) 
    sys.stderr.write('Euclidean distance to EAS data point is ' + str(round(edeas, args.dec_places)) + '\n' )
    edeur = ( (X-EUR_X)**2 + (Y - EUR_Y)**2 + (Z- EUR_Z)**2 )**(1/2) 
    sys.stderr.write('Euclidean distance to EUR data point is ' + str(round(edeur, args.dec_places)) + '\n' )
    edsas = ( (X-SAS_X)**2 + (Y - SAS_Y)**2 + (Z - SAS_Z)**2 )**(1/2) 
    sys.stderr.write('Euclidean distance to SAS data point is ' + str(round(edsas, args.dec_places)) + '\n' )
    #print(edafr, edeas, edeur, edsas)

    minind = [edafr, edeas, edeur, edsas].index(min([edafr, edeas, edeur, edsas]))
    ANC = ANCS[minind]
sys.stderr.write('=> Sample is ' + ANC + '\n' )

# OFNAME was specified initially 
if args.anc_prefix: OFNAME = os.path.join(os.path.dirname(OFNAME), ANC + '_' + os.path.basename(OFNAME))

PRS_SUM = 0
sys.stderr.write('### WRITING OUTPUT VCF\n' ) 

with open(OFNAME, 'w') as outfile:
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("##source=GC-HBOC-CanRisk-Pipeline_v" + VERSION + "\n")
    outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    outfile.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">\n')
    outfile.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Imputed Alternate Allele Dosage">\n')
    outfile.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Imputed Alternate Allele Dosage">\n')
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +  SAMPLE + "\n")
    
    for i in range(len(REF_VARS)):
        outfile.write(REF_VARS[i][0] + '\t')
        outfile.write(REF_VARS[i][1] + '\t.\t')
        outfile.write(REF_VARS[i][2] + '\t') # REF
        outfile.write(REF_VARS[i][3] + '\t.\t') # ALT & QUAL
        
        if FOUND[i]:
            tags = []
            if DP[i] < MINDP:
                tags.append('DP<' + str(MINDP))
            if not len(tags):
                tags.append("PASS")
            outfile.write(';'.join(tags) + '\t') # FILTER & INFO
            if DP[i] < MINDP:
                anc_ind = ANCS.index(ANC)
                outfile.write('.\tGT:DS\t./.:' + str(AFS[REF_VARS[i]][anc_ind] * 2) )
                PRS_SUM += AFS[REF_VARS[i]][anc_ind] * 2 * WEIGHTS[i]
            else:
                outfile.write('DP=' + str(DP[i]) + '\t')
                if not REVERSE_DICT[REF_VARS[i]]:
                    outfile.write('GT\t' + GT[i])
                    if GT[i] == '0/1': PRS_SUM += WEIGHTS[i]
                    elif GT[i] == '1/1': PRS_SUM += (2 * WEIGHTS[i])
                else:
                    if GT[i] == '0/0':
                        outfile.write('GT\t1/1')
                        PRS_SUM += 2 * WEIGHTS[i]
                    elif GT[i] == '0/1':
                        outfile.write('GT\t0/1')
                        PRS_SUM += WEIGHTS[i]
                    else:
                         outfile.write('GT\t0/0')          
                    
        else:
            outfile.write("NotGenotyped\t.\t") # FILTER & INFO
            anc_ind = ANCS.index(ANC)
            outfile.write('GT:DS\t./.:' + str(AFS[REF_VARS[i]][anc_ind] * 2) )
            PRS_SUM += AFS[REF_VARS[i]][anc_ind] * 2 * WEIGHTS[i]
        outfile.write('\n')



sys.stderr.write('=> OUTPUT VCF written to ' + OFNAME + '\n' )

sys.stderr.write('### PRS\n' )
sys.stderr.write('=> Raw PRS is ' + str(round(PRS_SUM, args.dec_places)) + '\n')
ZSCORE = None
if ANC == "AFR":
    if AFR_MEAN == None and AFR_SD == None:
        sys.stderr.write("Mean and standard deviation for AFR PRS unkonwn\n")
    elif AFR_MEAN == None:
        sys.stderr.write("Mean for AFR PRS unkonwn\n")
    elif AFR_SD == None:
        sys.stderr.write("Standard deviation for AFR PRS unkonwn\n")
    else: ZSCORE = (PRS_SUM - AFR_MEAN)/AFR_SD
elif ANC == "EAS":
    if EAS_MEAN == None and EAS_SD == None:
        sys.stderr.write("Mean and standard deviation for EAS PRS unkonwn\n")
    elif EAS_MEAN == None:
        sys.stderr.write("Mean for EAS PRS unkonwn\n")
    elif EAS_SD == None:
        sys.stderr.write("Standard deviation for EAS PRS unkonwn\n")
    else: ZSCORE = (PRS_SUM - EAS_MEAN)/EAS_SD
elif ANC == "EUR":
    if EUR_MEAN == None and EUR_SD == None:
        sys.stderr.write("Mean and standard deviation for EUR PRS unkonwn\n")
    elif EUR_MEAN == None:
        sys.stderr.write("Mean for EUR PRS unkonwn\n")
    elif EUR_SD == None:
        sys.stderr.write("Standard deviation for EUR PRS unkonwn\n")
    else: ZSCORE = (PRS_SUM - EUR_MEAN)/EUR_SD
elif ANC == "SAS":
    if SAS_MEAN == None and SAS_SD == None:
        sys.stderr.write("Mean and standard deviation for SAS PRS unkonwn\n")
    elif SAS_MEAN == None:
        sys.stderr.write("Mean for SAS PRS unkonwn\n")
    elif SAS_SD == None:
        sys.stderr.write("Standard deviation for SAS PRS unkonwn\n")
    else: ZSCORE = (PRS_SUM - SAS_MEAN)/SAS_SD   
if ZSCORE: sys.stderr.write("=> Normalized z-score is " + str(round(ZSCORE, args.dec_places)) +"\n")


