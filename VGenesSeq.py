__author__ = 'PCW-MacBookProRet'
import os, sys
import IsoelectricPoint, IUPACData, ProtParamData
import sqlite3 as db
import VGenesSQL
from PyQt5.QtWidgets import QApplication

global working_prefix
global temp_folder
global bin_folder

working_prefix = os.path.dirname(os.path.realpath(sys.argv[0]))
temp_folder = os.path.join(working_prefix, 'Temp')
bin_folder = os.path.join(working_prefix, 'Tools')

def Translator(Sequence, frame):

        # Translate sequence into a list of codons
    CodonList = [ ]
    for x in range(frame, len(Sequence), 3):
            CodonList.append(Sequence[x:x+3])
    # For each codon, translate it into amino acid, and create a list
    ProteinSeq = [ ]
    for codon in CodonList:
        if codon in CodonDict:
            ProteinSeq.append(CodonDict[codon])
        else:
            ProteinSeq.append('~')

    AASeq = ''.join(ProteinSeq)

    # print("Translated in frame %d: %s (%.1f Da)" % ((frame+1), ''.join(ProteinSeq), sum(ProteinWeight)))
    # Check position of stop codon, making sure it's at the end, and the only one
    XCount = 0
    UCount = 0
    for acid in ProteinSeq:
        if acid == "*":
            XCount += 1
    for acid in ProteinSeq:
        if acid == "~":
            UCount += 1
    ErMessage = []
    # ErMessage.append()
    if XCount > 0:
        if XCount == 1:
            ErMes =  'WARNING: '+ str(XCount) + ' stop codon was found (marked as "*")!'
        else:
            ErMes =  'WARNING: '+ str(XCount) + ' stop codons found (marked as "*")!'
        ErMessage.append(ErMes)
    if UCount > 0:
        # todo this doesn't label errors properly
        AASeq2 = AASeq.replace ('.', '')
        ErMes = 'Codon errors (marked as "~"): '
        if len(Sequence) % 3 != 0 and UCount == 1:
            ErMes += 'Incomplete codon at end.'
            ErMessage.append(ErMes)
            return AASeq, ErMessage

        elif UCount == 1:

            if AASeq2[0] == '~':
                ErMes += 'The first codon is incomplete.'
                ErMessage.append(ErMes)
                return AASeq, ErMessage

            else:
                ErMes += '1 codon error internally.'
                ErMessage.append(ErMes)
                return AASeq, ErMessage

        elif UCount > 1:


            if AASeq2[0] == '~':
                ErMes += 'The first codon is incomplete. '
                ErMessage.append(ErMes)
                UCount -= 1

            if len(Sequence) % 3 != 0:
                if UCount > 1:
                    ErMes += '1 incomplete on end and '
                    if UCount-1 > 1:
                        ErMes += str(UCount-1) + ' others with errors internally.'
                    elif UCount - 1 == 1:
                        ErMes += '1 other with errors internally.'
                else:
                    ErMes += '1 incomplete on end.'

            else:
                ErMes += str(UCount) + ' errors within the sequence.'
        ErMessage.append(ErMes)

    return AASeq, ErMessage
# AminoDict={'A':89.09,   'R':174.20, 'N':132.12, 'D':133.10, 'C':121.15,
# 'Q':146.15, 'E':147.13, 'G':75.07,  'H':155.16, 'I':131.17, 'L':131.17,
# 'K':146.19, 'M':149.21, 'F':165.19, 'P':115.13, 'S':105.09, 'T':119.12,
# 'W':204.23, 'Y':181.19, 'V':117.15, 'X':0.0,    '-':0.0,    '*':0.0,
# '?':0.0}

CodonDict={'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',
'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',
'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',
'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',
'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',
'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',
'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',
'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',
'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',
'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',
'AGA':'R',  'AGG':'R',  'TAA':'*',  'TAG':'*',  'TGA':'*',  '...':'.',
'NNN':'.'}

def isoelectric_point(AASeq):
    """Calculate the isoelectric point.

    Uses the module IsoelectricPoint to calculate the pI of a protein.
    """
    aa_content = count_amino_acids(AASeq) #dictionary:  {AA:number present}

    # ie_point = IsoelectricPoint.IsoelectricPoint(AASeq.sequence, aa_content)
    ie_point = IsoelectricPoint.IsoelectricPoint(AASeq, aa_content)
    return ie_point.pi()

def count_amino_acids(AASeq):
    """Count standard amino acids, returns a dict.

    Counts the number times each amino acid is in the protein
    sequence. Returns a dictionary {AminoAcid:Number}.

    The return value is cached in self.amino_acids_content.
    It is not recalculated upon subsequent calls.
    """
    # if AASeq.amino_acids_content is None:
    prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
    for aa in prot_dic:
        prot_dic[aa] = AASeq.count(aa)

    amino_acids_content = prot_dic

    return amino_acids_content

def molecular_weight(AAseq):
    """Calculate MW from Protein sequence"""
    # make local dictionary for speed
    # if AAseq.monoisotopic:
    #     water = 18.01
    #     iupac_weights = IUPACData.monoisotopic_protein_weights
    # else:
    iupac_weights = IUPACData.protein_weights
    water = 18.02

    aa_weights = {}
    for i in iupac_weights:
        # remove a molecule of water from the amino acid weight
        aa_weights[i] = iupac_weights[i] - water

    total_weight = water  # add just one water molecule for the whole sequence
    Failed = False
    for aa in AAseq:   #.sequence:
        try:
            total_weight += aa_weights[aa]
        except:
            failed = True
    if Failed == True: total_weight = 0
    return total_weight

def get_amino_acids_percent(AAseq):
        """Calculate the amino acid content in percentages.

        The same as count_amino_acids only returns the Number in percentage of
        entire sequence. Returns a dictionary of {AminoAcid:percentage}.

        The return value is cached in self.amino_acids_percent.

        input is the dictionary self.amino_acids_content.
        output is a dictionary with amino acids as keys.
        """
        # if AAseq.amino_acids_percent is None:
        aa_counts = count_amino_acids(AAseq)

        percentages = {}
        for aa in aa_counts:
            percentages[aa] = aa_counts[aa] / float(AAseq.length)

        amino_acids_percent = percentages

        return amino_acids_percent


def OtherParam(AASequence, Param, WindowSize, SupressErrors):
    import ProtParamData
    from PyQt5 import QtWidgets

    if AASequence.count('*') > 0:
        if SupressErrors != True:
            QtWidgets.QMessageBox.critical(None, "Sequence error",
            "The sequence contains a stop codon.\n",QtWidgets.QMessageBox.Cancel)
        return 0

    # todo if windo is 0 then whole seqeunce can use for report


    if WindowSize <= 0 or WindowSize > (len(AASequence)-1):
        WindowSize = len(AASequence)

    # todo code to put in output: AACounts, AAPercents, AAMW, Instability, flex, Hydropobicity, Hydrophilicity, Surface (last four are for graphs)
    if Param == 'AACounts':
        RetValue = count_amino_acids(AASequence) #returns a dictionary in form {AA: number in seq}
    elif Param == 'AAPercents':
        RetValue = get_amino_acids_percent(AASequence)  #returns a dictionary in form {AA: percentage in seq}
    elif Param == 'AAMW': #molecular weight
        RetValue = molecular_weight(AASequence)

    elif Param == 'AApI': #isoelectric point

        RetValue = isoelectric_point(AASequence)

    elif Param == 'MapAApI': #isoelectric point

        preMAP= []
        for i in range(len(AASequence) - WindowSize + 1):
            subsequence = AASequence[i:i+WindowSize]
            score = 0.0

            preMAP.append(isoelectric_point(subsequence))
        #     this is the pI for each peptide at this window length,
        # now need to adjust based on surrounding

        weights = _weight_list(WindowSize, 0.5)
        scores = []

        # the score in each Window is divided by the sum of weights
        # (* 2 + 1) since the weight list is one sided:
        sum_of_weights = sum(weights) * 2 + 1

        for i in range(WindowSize//2): #first add unweighted from end as won't be positoned otherwise
            scores.append(preMAP[i])

        # now add all weighted ones that can be
        i = WindowSize // 2
        for i in range(WindowSize//2,len(preMAP) - WindowSize + 1):
            # subsequence = AASequence[i:i+WindowSize]
            score = 0.0

            for j in range(WindowSize // 2):
                # walk from the outside of the Window towards the middle.
                # Iddo: try/except clauses added to avoid raising an exception on a non-standard amino acid

                try:
                    front = preMAP[(i-WindowSize//2)+j]       #  param_dict[subsequence[j]]
                    back =  preMAP [(i+WindowSize//2)-j]                   #  param_dict[subsequence[WindowSize - j - 1]]
                    score += weights[j] * front + weights[j] * back
                except KeyError:
                    sys.stderr.write('warning: problem with sequence')
                    #          (subsequence[j], subsequence[WindowSize - j - 1]))

            # Now add the middle value, which always has a weight of 1.
            score += preMAP[i]




            scores.append(score / sum_of_weights)
        difference = len(preMAP)-len(scores)
        for i in range(len(preMAP)-difference, len(preMAP)):
            scores.append(preMAP[i])


        return scores


    elif Param == 'Instability':  # can use to see if antibody might not express well
        """Calculate the instability index according to Guruprasad et al 1990.

        Implementation of the method of Guruprasad et al. 1990 to test a
        protein for stability. Any value above 40 means the protein is unstable
        (has a short half life).

        See: Guruprasad K., Reddy B.V.B., Pandit M.W.
        Protein Engineering 4:155-161(1990).
        """
        RetValue = instability_index(AASequence)

    elif Param == 'MapInstability':  # can use to see if antibody might not express well
        scores = []
        score = 0.0
        # """Calculate the flexibility according to Vihinen, 1994.
        if WindowSize <4: WindowSize = 4
        for i in range(len(AASequence) - WindowSize + 1):
            subsequence = AASequence[i:i+WindowSize]

            score = instability_index(subsequence)
            scores.append(score)
        RetValue = scores

    elif Param == 'Flexibility':
        scores = []
        score = 0.0
        # """Calculate the flexibility according to Vihinen, 1994.
        for i in range(len(AASequence) - 9 + 1):
            subsequence = AASequence[i:i+9]

            score = flexibility(subsequence)
            scores.append(score)
        RetValue = scores

    elif Param == 'Hydrophobicity':     # Kyte & Doolittle index of hydrophobicity

        RetValue = protein_scale(AASequence, param_dict=ProtParamData.kd,window=WindowSize,edge=0.5)
        #  returns a a scaled list of hydrophobicity that can be graphed

    elif Param == 'Hydrophilicity':     # 1 Hopp & Wood, Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).

        RetValue = protein_scale(AASequence, param_dict=ProtParamData.hw,window=WindowSize,edge=0.5)

        # returns a a scaled list of Hydrophilicity that can be graphed
        # 1 Hopp & Wood, Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).

    elif Param == 'Surface': # Emini Surface fractional probability...based on average effects on if antibody can target

        RetValue = protein_scale(AASequence, param_dict=ProtParamData.em,window=WindowSize,edge=0.5)
        # returns a a scaled list of Hydrophilicity that can be graphed

    # todo make mutability indices...that is C, G, A, and T in hotspots and/or based on silent frequencies from large pool
    # todo can predict if selected mutations based on large pool silent for particular V gene and HS/CS..use on clones

    return RetValue

def protein_scale(AAseq, param_dict, window, edge):
    """Compute a profile by any amino acid scale.

    An amino acid scale is defined by a numerical value assigned to each type of
    amino acid. The most frequently used scales are the hydrophobicity or
    hydrophilicity scales and the secondary structure conformational parameters
    scales, but many other scales exist which are based on different chemical and
    physical properties of the amino acids.  You can set several parameters that
    control the computation  of a scale profile, such as the window size and the
    window edge relative weight value.

    WindowSize: The window size is the length
    of the interval to use for the profile computation. For a window size n, we
    use the i-(n-1)/2 neighboring residues on each side to compute
    the score for residue i. The score for residue i is the sum of the scaled values
    for these amino acids, optionally weighted according to their position in the
    window.

    Edge: The central amino acid of the window always has a weight of 1.
    By default, the amino acids at the remaining window positions have the same
    weight, but you can make the residue at the center of the window  have a
    larger weight than the others by setting the edge value for the  residues at
    the beginning and end of the interval to a value between 0 and 1. For
    instance, for Edge=0.4 and a window size of 5 the weights will be: 0.4, 0.7,
    1.0, 0.7, 0.4.

    The method returns a list of values which can be plotted to
    view the change along a protein sequence.  Many scales exist. Just add your
    favorites to the ProtParamData modules.

    Similar to expasy's ProtScale: http://www.expasy.org/cgi-bin/protscale.pl
    """
    # generate the weights
    #   _weight_list returns only one tail. If the list should be [0.4,0.7,1.0,0.7,0.4]
    #   what you actually get from _weights_list is [0.4,0.7]. The correct calculation is done
    #   in the loop.
    weights = _weight_list(window, edge)
    scores = []

    # the score in each Window is divided by the sum of weights
    # (* 2 + 1) since the weight list is one sided:
    sum_of_weights = sum(weights) * 2 + 1

    for i in range(len(AAseq) - window + 1):
        subsequence = AAseq[i:i+window]
        score = 0.0

        for j in range(window // 2):
            # walk from the outside of the Window towards the middle.
            # Iddo: try/except clauses added to avoid raising an exception on a non-standard amino acid
            try:
                front = param_dict[subsequence[j]]
                back = param_dict[subsequence[window - j - 1]]
                score += weights[j] * front + weights[j] * back
            except KeyError:
                sys.stderr.write('warning: %s or %s is not a standard amino acid.\n' %
                         (subsequence[j], subsequence[window - j - 1]))

        # Now add the middle value, which always has a weight of 1.
        middle = subsequence[window // 2]
        if middle in param_dict:
            score += param_dict[middle]
        else:
            sys.stderr.write('warning: %s  is not a standard amino acid.\n' % (middle))

        scores.append(score / sum_of_weights)

    return scores

def _weight_list(window, edge):
    """Makes a list of relative weight of the
    window edges compared to the window center. The weights are linear.
    it actually generates half a list. For a window of size 9 and edge 0.4
    you get a list of [0.4, 0.55, 0.7, 0.85].
    """
    unit = 2 * (1.0 - edge) / (window - 1)
    weights = [0.0] * (window // 2)

    for i in range(window // 2):
        weights[i] = edge + unit * i

    return weights

def instability_index(AAseq):
    """Calculate the instability index according to Guruprasad et al 1990.

    Implementation of the method of Guruprasad et al. 1990 to test a
    protein for stability. Any value above 40 means the protein is unstable
    (has a short half life).

    See: Guruprasad K., Reddy B.V.B., Pandit M.W.
    Protein Engineering 4:155-161(1990).
    """
    index = ProtParamData.DIWV
    score = 0.0

    for i in range(len(AAseq) - 1):
        this, next = AAseq[i:i+2]
        dipeptide_value = index[this][next]
        score += dipeptide_value

    return (10.0 / len(AAseq)) * score

def flexibility(subsequence):
    """Calculate the flexibility according to Vihinen, 1994.

    No argument to change window size because parameters are specific for a
    window=9. The parameters used are optimized for determining the flexibility.
    """
    flexibilities = ProtParamData.Flex
    window_size = 9
    weights = [0.25, 0.4375, 0.625, 0.8125, 1]
    # scores = []


    score = 0.0

    for j in range(window_size // 2):
        front = subsequence[j]
        back = subsequence[window_size - j - 1]
        score += (flexibilities[front] + flexibilities[back]) * weights[j]

    middle = subsequence[window_size // 2 + 1]
    score += flexibilities[middle]

    score = score / 5.25

    return score




def ClustalO(SeqDict, wrapLength, ordered):
    # input is a list of lists containing filename and sequence: ((fielname1, seq1),(fielname2, seq2))
    # import time
    # clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v
    CurDir = os.getcwdb()
    # workingfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'ClustalOmega', 'my-in-seqs.fa')
    # workingfilename = '/Applications/ClustalOmega/my-in-seqs.fa'
    # workingdir, NameBase = os.path.split(DBname)
    #
    #
    #
    # NameBase = NameBase[:(len(NameBase) - 4)]

    import uuid

    NameBase  =str(uuid.uuid4())
    # NameBase = NameBase[:12]
    NameBase = NameBase.replace('-', '')

    NameBase  = NameBase.replace(' ', '')


    MyInFiles = NameBase + 'In.fa'
    MyOutFiles = NameBase + 'Out.fa'

    # workingfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-in-seqs.fa')
    # savefilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-out-seqs.fa')

    workingfilename = os.path.join(temp_folder, MyInFiles)
    savefilename = os.path.join(temp_folder, MyOutFiles)


    workingdir, filename = os.path.split(workingfilename)




    os.chdir(workingdir)
    i = 1
    FASTAfile = ''
    for seq in SeqDict:
        for item in seq:
            if i % 2 != 0:
                FASTAfile += '>' + item + '\n'
                if i == 1: SeqName = FASTAfile
            else:
                if item == '': item = 'N'
                if SeqName == "Germline":
                    if len(item) < 50:
                        return False

                FASTAfile += item + '\n'
            i += 1

    with open(workingfilename,
              'w') as currentFile:  # using with for this automatically closes the file even if you crash
        currentFile.write(FASTAfile)
    # todo with the ./ in front of clustalO all I need is the 'stand alone MAC binary' from http://www.clustal.org/omega/ no other installation
    # Mac users should rename the downloaded file to clustalo and
    #  place in the location of their choice. This file may need to
    # be made executable e.g.: chmod u+x clustalo
    # ClustalOCommandLine = "clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v --full --force"  # --full uses better alignment but slower
    # with open('/Applications/ClustalOmega/my-out-seqs.fa', 'w') as currentFile:  #to set a test phrase so that it only progresses after clustal is doen
    with open(savefilename, 'w') as currentFile:
        currentFile.write('The clustal output is not done yet')

    if ordered == True:
        ClustalOCommandLine = bin_folder + '/clustalo -i '+ MyInFiles + ' -o '+ MyOutFiles + ' -v --force --output-order=tree-order --outfmt=vie --resno --wrap=' + str(wrapLength)  #
    else:
        ClustalOCommandLine = bin_folder + '/clustalo -i '+ MyInFiles + ' -o '+ MyOutFiles + ' -v --force --output-order=tree-order --outfmt=vie --resno --wrap=1000'
    ClustalOut = os.popen(ClustalOCommandLine)

    # if ordered == True:
    #     ClustalOCommandLine = "./clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v --force --output-order=tree-order --outfmt=vie --resno --wrap=" + str(
    #         wrapLength)  #
    # else:
    #
    #     ClustalOCommandLine = './clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v --force --output-order=tree-order --outfmt=vie --resno --wrap=1000'
    # ClustalOut = os.popen(ClustalOCommandLine)







    ClustalOutLines = []
    for line in ClustalOut:
        ClustalOutLines.append(line)

    send = False

    while send == False:
        # with open('/Applications/ClustalOmega/my-out-seqs.fa', 'r') as currentFile:  #using with for this automatically closes the file even if you crash
        with open(savefilename, 'r') as currentFile:
            testName = currentFile.readline()
            if testName != 'The clustal output is not done yet': send = True
            if send == False and len(ClustalOutLines) < 2:
                return False
                # print('first')



    os.remove(workingfilename)

    return savefilename

#
# Alignment Output ClustalO:
#   -o, --out, --outfile={file,-} Multiple sequence alignment output file (default: stdout)
#   --outfmt={a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]} MSA output file format (default: fasta)
#   --residuenumber, --resno  in Clustal format print residue numbers (default no)
#   --wrap=<n>                number of residues before line-wrap in output
#   --output-order={input-order,tree-order} MSA output order like in input/guide-tree
def ReverseComp(Sequence):
    LenSeq = len(Sequence)
    Sequence = Sequence.upper()
    Sequence = Sequence.strip()
    NewSeq = ''
    for i in range(LenSeq-1,0, -1):
        Ch = Sequence[i]
        if Ch == 'A':
            NewSeq += 'T'
        elif Ch == 'T':
            NewSeq += 'A'
        elif Ch == 'G':
            NewSeq += 'C'
        elif Ch == 'C':
            NewSeq += 'G'
        else:
            NewSeq += 'N'
    return NewSeq



def readClustalOutput(outfilename):
    # reads a FASTA file and returns a list of tupels with seqname and then seqeunce
    # useful for building alignments, etc.

    ReadFile = []
    # ReadFile.clear
    # SeqRead = []
    Readline = ''
    currentFile2 = ''
    with open(outfilename, 'r') as currentFile2:  #using with for this automatically closes the file even if you crash
        # currentFile.write(FASTAfile)
        Seq  = ''
        SeqName = ''
        Readline = ''
        for line in currentFile2:
            Readline = line.replace('\n', '').replace('\r', '')

            if Readline[0] == '>':
                if Seq != '':  #saves the previous except on first round
                    SeqRead = (SeqName, Seq)
                    ReadFile.append(SeqRead)
                    Seq = ''
                SeqName = Readline[1:]

            else:
                Seq += Readline
        SeqRead = (SeqName, Seq) # must save last one at end
        ReadFile.append(SeqRead)


    # os.remove(workingfilename)
    # os.chdir(CurDir)

    # Returns a list of seqname and sequences, but now aligned
    return ReadFile



def Mutations(Alignment, GVend):
    global StartID
    global IDLength
    # todo junctions should be identical so won't matter
    GermSeq = ''
    TestSeq = ''
    if Alignment is None: return ''
    for sequence in Alignment:

        if sequence[0] == 'Germline':
            GermSeq = sequence[1]
        else:
            TestSeq  = sequence[1]

    GermSeq = GermSeq.upper()
    TestSeq = TestSeq.upper()
    AlignedSeq = TestSeq
    MutList = []
    IDLength = 0
    try:
        while GermSeq[len(GermSeq)-1] == '-':
            GermSeq = GermSeq[:len(GermSeq)-1]  #  code to remove dashes from end where seq is matched to C region
    except:
        return MutList, TestSeq, IDLength


    while TestSeq[len(TestSeq)-1] == '-':
        TestSeq = TestSeq[:len(TestSeq)-1]  #  code to remove dashes from end where seq is matched to C region

    #     todo need remove '-' from front also and remove that length from other seq
    i = 0
    TrimFront = 0
    for i in range(len(GermSeq)):
        if GermSeq[i] == '-' or TestSeq[i] == '-':
            TrimFront += 1
        else:
            break
    if TrimFront > 0:
        GermSeq = GermSeq[TrimFront:]
        TestSeq = TestSeq[TrimFront:]

    if TestSeq.count('-') > 0 or GermSeq.count('-') > 0:
        # TestSeq, GermSeq = IDPlacer(TestSeq, GermSeq)
        TestSeq = IDPlacer(TestSeq, GermSeq)
    i = 0
    # GermSeq = ''
    # TestSeq = ''
    GermNuc = ''
    TestNuc = ''
    event = ''
    IDSeq = ''
    IDnuc = ''

    MutIs = []

    if len(GermSeq) < len(TestSeq):
        SeqLen = len(GermSeq)
    else:
        SeqLen = len(TestSeq)

    StartID = True
    IDLength = 0
    for i in range(i,SeqLen-1):
        GermNuc = GermSeq[i]
        TestNuc = TestSeq[i]
        if GermNuc != TestNuc and GermNuc != 'N' and TestNuc != 'N':
            if GermNuc == '-':
                if GermSeq[i-1] == '-':
                    StartID = False
                else:
                    StartID = True
                if StartID == True:
                    event  = 'Insertion'
                    pos = i
                    IDSeq = ''
                    IDnuc = GermSeq[i]
                    j = i
                    while IDnuc == '-':
                        IDSeq += TestSeq[j]
                        j +=1
                        IDnuc = GermSeq[j]
                        if pos < GVend:
                            IDLength += 1

                    MutIs = (event, (pos+1), IDSeq)
                    MutList.append(MutIs)

            elif TestNuc == '-':
                if TestSeq[i-1] == '-':
                    StartID = False
                else:
                    StartID = True
                if StartID == True:
                    event  = 'Deletion'
                    pos = i
                    IDSeq = ''
                    IDnuc = TestSeq[i]
                    j = i
                    while IDnuc == '-':
                        IDSeq += GermSeq[j]
                        j +=1
                        IDnuc = TestSeq[j]
                        # if pos < GVend:
                        #     IDLength -= 1

                    MutIs = (event, (pos+1), IDSeq)
                    MutList.append(MutIs)

            else:
                StartID = False

                baseFrom = GermNuc
                baseTo = TestNuc
                pos = i


                MutIs = (baseFrom, (pos+1), baseTo)


                MutList.append(MutIs)
                #

        # else:
            # i += 1
        # EditedSeq = 'a'
    return MutList, TestSeq, IDLength, AlignedSeq

def RScaller(MutList, Vgene, species, DBpathname):
    #     code to determine if mutations are replacements or silent
    from math import ceil
    GVseq = []
    try:
        DBpathname = os.path.join(working_prefix, 'VDJGenes.db')

        (dirname, filename) = os.path.split(DBpathname)
        os.chdir(dirname)

        GetProductive = False
        conn = db.connect(DBpathname)


    except:
        DBpathname = '/Volumes/Promise Pegasus/Dropbox/VGenes/VDJGenes.db'
        (dirname, filename) = os.path.split(DBpathname)
        os.chdir(dirname)

        GetProductive = False
        conn = db.connect(DBpathname)

    # then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()

    if species == 'Human':
        SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Vgene + '"'
        # SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'

    elif species == 'Mouse':
        SqlStatementV = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Vgene + '"'
        # SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'

    GVseq.clear()

    cursor.execute(SqlStatementV)
    for row in cursor:
        for column in row:
            GVseq.append(column)  # 0 SeqName, 1 Allele, 2 Species, 3 CodingStartNucleotide, 4 Sequence, 5 IMGTSequence
    VGeneLocus = GVseq[0]
    GVgene = GVseq[4]
    GStart = int(GVseq[3]) - 1
    GVgene = GVgene[GStart:]
    Codons  = split_str(GVgene, 3, skip_tail=True)  # give back list of string split by three


    SeqMuts = MutList.split(',')

    NumMuts = len(SeqMuts)
    Codon = ''
    MutCodon = ''
    if NumMuts > 0:
        for j in range(0, NumMuts-1):
            Mutation = SeqMuts[j]
            MutDet = Mutation.split('-')
            MutNuc = int(MutDet[1])
            CodonNum = ceil(MutNuc/3)-1 #with base 0 so need -1
            MutPos = MutNuc - ((CodonNum)*3)
            try:
                Codon = Codons[CodonNum]

                # MutPos = (MutNuc - ((CodonNum)*3))#-1

                base = MutDet[2]
                newCod = list(Codon)
                newCod[MutPos-1] = base #need -1 for base 0
                MutCodon = ''.join(newCod)
                # MutCodon[MutPos] = base
            except:
                print('nope')
            if Codon in CodonDict:
                GermAA = CodonDict[Codon]
            else:
                GermAA = '?'


            if MutCodon in CodonDict:
                MutAA = CodonDict[MutCodon]
            else:
                MutAA = '?'

            if MutAA == '?' or GermAA == '?':
                RS = '?'
            elif MutAA == GermAA:
                RS = 'S'
            else:
                RS = 'R' + GermAA + MutAA

            NewMut = Mutation + '-' + RS
            SeqMuts[j] = NewMut



    return SeqMuts



def split_str(seq, chunk, skip_tail=False):
    lst = []
    if chunk <= len(seq):
        lst.extend([seq[:chunk]])
        lst.extend(split_str(seq[chunk:], chunk, skip_tail))
    elif not skip_tail and seq:
        lst.extend([seq])
    return lst

def Intraclonal(DataIn, DBFilename):
    # DataIn: list of tuples: (SeqName, Sequence, ClonalPool, GermlineSequence, Mutations, Vbeg, Vend, species, Vgene)
    # produce a CSV file: Name1, Name2, #same, #different, Sequence-skelton, sequence

    #			fields = ['SeqName', 'Sequence', 'ClonalPool', 'GermlineSequence', 'Mutations', 'GVbeg', 'GVend', 'Species',
	#		          'V1', 'Project', 'Grouping', 'Isotype', 'Quality' , 'Subspecificity']
    from operator import itemgetter
    # from collections import Counter
    import itertools
    # import csv

    DataIn.sort(key=itemgetter(2,4))
    ClonalPool = []
    ClonalPools = []
    SimDif = ()


    for k,v in itertools.groupby(DataIn, key=itemgetter(2)):   #first split out seperate clonal pools
        i = int(k)

        if i != 0:
            for item in v:
                 ClonalPool.append(item)
            CurrentPool = tuple(ClonalPool)
            ClonalPools.append(CurrentPool)
            ClonalPool.clear()

    MutRep = []

    compMutsR = []
    compMutsS = []
    SeqMutsR = []
    SeqMutsS = []
    RDifferences2 = []
    SDifferences2 = []
    RMatches2 = []
    SMatches2 = []
    Warn = ' '
    for pool in ClonalPools:
        lenPool = len(pool)
        for i in range(0, lenPool-1):
            seq = pool[i]
            SeqMutstr = seq[4]



            SeqMuts = SeqMutstr.split(',')
            Vbeg = int(seq[5])
            Vend = int(seq[6])
            NumMuts = len(SeqMuts)
            # if NumMuts > 0:
            #     for j in range(0, NumMuts-1):
            #         Mutation = SeqMuts[j]
            #         MutDet = Mutation.split('-')
            #         MutPos = MutDet[1]
            #         AdjGLPos = int(MutPos) + (Vbeg)
            #         AdjGL = MutDet[0] + '-' + str(AdjGLPos) + '-' + MutDet[2]
            #         SeqMuts[j] = AdjGL
            #
            #     SeqMutstr = ','.join(SeqMuts)

            SeqMuts = RScaller(SeqMutstr, seq[8], seq[7], DBFilename)
            # SeqMuts = RScaller(SeqMutstr, [VH gene], [species], DBFilename)





            for j in range(1, lenPool):

                compSeq = pool[j]
                compMutstr = compSeq[4]

                compMuts = compMutstr.split(',')
                CVbeg = int(compSeq[5])
                CVend = int(compSeq[6])
                NumMuts = len(compMuts)
                # if NumMuts > 0:
                #     for k in range(0, NumMuts - 1):
                #         Mutation = compMuts[k]
                #         MutDet = Mutation.split('-')
                #         MutPos = MutDet[1]
                #         AdjGLPos = int(MutPos) + (CVbeg)
                #         AdjGL = MutDet[0] + '-' + str(AdjGLPos) + '-' + MutDet[2]
                #         compMuts[k] = AdjGL
                #
                #     compMutstr = ','.join(compMuts)

                compMuts = RScaller(compMutstr, compSeq[8], compSeq[7], DBFilename)
                Warn = ' '

                if seq[8] != compSeq[8]:
	                Warn = 'V genes mismatched: ' + seq[8] + ' x ' + compSeq[8]

                if Vbeg >= CVbeg:
                    beg = Vbeg
                else:
                    beg = CVbeg

                if Vend <= CVend:
                    end = Vend
                else:
                    end = CVend

                compMutsS.clear()
                SeqMutsS.clear()

                compMutsR.clear()
                SeqMutsR.clear()
                for k in range(0, NumMuts - 1):
                    Mutation = compMuts[k]
                    MutDet = Mutation.split('-')
                    RS = MutDet[3]
                    MutPos = int(MutDet[1])
                    if MutPos >= beg and MutPos <= end:
                        if RS[0] == 'R':
                            compMutsR.append(Mutation)
                        elif RS[0] == 'S':
                            compMutsS.append(Mutation)

                NumMuts = len(SeqMuts)
                for k in range(0, NumMuts - 1):
                    Mutation = SeqMuts[k]
                    MutDet = Mutation.split('-')
                    RS = MutDet[3]
                    MutPos = int(MutDet[1])
                    if MutPos >= beg and MutPos <= end:
                        if RS[0] == 'R':
                            SeqMutsR.append(Mutation)
                        elif RS[0] == 'S':
                            SeqMutsS.append(Mutation)

                RMatches = [element for element in SeqMutsR if element in compMutsR]
                RDifferences = [element for element in SeqMutsR if element not in compMutsR]
                RDifferences += [element for element in compMutsR if element not in SeqMutsR]

                #			fields = ['SeqName', 'Sequence', 'ClonalPool', 'GermlineSequence', 'Mutations', 'GVbeg', 'GVend', 'Species',
                #		          'V1', 'Project', 'Grouping', 'Isotype', 'Quality' , 'Subspecificity']

                SMatches = [element for element in SeqMutsS if element in compMutsS]
                SDifferences = [element for element in SeqMutsS if element not in compMutsS]
                SDifferences += [element for element in compMutsS if element not in SeqMutsS]
                # Differences.sort()
                Name1 = seq[0]
                Name2 = compSeq[0]
                Project = seq[9]
                try:
                    Subject = seq[10]
                except:
                    print('n')
                PreQuality = seq[13]
                Age = seq[14]
                if Project == 'Heavy-Aged':
                    Quality  = 'Aged'
                elif Project == 'Heavy-Young':
                    if PreQuality == 'pH1N1':
                        Quality = 'pH1N1'
                    else:
                        Quality = 'prePandemic'
                else:
                    Quality = 'None'


                Act  = seq[12]
                CompActivity = compSeq[12]
                if Act == 'HAI' or CompActivity == 'HAI':
                    Activity  = 'HAI'
                elif Act == 'MN' or CompActivity == 'MN':
                    Activity  = 'MN'
                else:
                    Activity  = 'Binding'

                # if Delin == True:
                #     CP = Subject + '_' + seq[2]
                # else:
                CP = seq[2]

                lengthCompared = (end - beg)+1
                Adjustment = (end/lengthCompared)
                RAdjMatch = Adjustment * len(RMatches)
                RAdjDif = Adjustment * len(RDifferences)
                SAdjMatch = Adjustment * len(SMatches)
                SAdjDif = Adjustment * len(SDifferences)

                TotDif = len(RDifferences) + len(SDifferences)
                TotAdjDif = Adjustment * TotDif

                TotMatch = len(SMatches) + len(RMatches)
                TotAdjMatch = Adjustment * TotMatch
                RDifferences2.clear()
                SDifferences2.clear()
                RMatches2.clear()
                SMatches2.clear()
                # RDifferences2.append('R-Differences: ')
                # SDifferences2.append('S-Differences: ')
                # RMatches2.append('R-RMatches: ')
                # SMatches2.append('S-Matches: ')


                for Muta in RDifferences:
	                RDifferences2.append(Muta + '|')
                for Muta in SDifferences:
	                SDifferences2.append(Muta + '|')
                for Muta in RMatches:
	                RMatches2.append(Muta + '|')
                for Muta in SMatches:
	                SMatches2.append(Muta + '|')
                RDif = ''.join(RDifferences2)
                RMatch = ''.join(RMatches2)
                SDif = ''.join(SDifferences2)
                SMatch = ''.join(SMatches2)
                Comparison = Name1 + '_x_' + Name2
                SimDif = (
                Comparison, Project, Subject, Quality, CP, Name1, Name2, Activity, TotDif, len(RDifferences), len(SDifferences), beg, end, lengthCompared, TotMatch, int(TotAdjMatch), int(TotAdjDif), len(RMatches),
                int(RAdjMatch), int(RAdjDif), len(SMatches), int(SAdjMatch),
                int(SAdjDif), Warn, RDif, RMatch, SDif, SMatch, Age)

                if Name1 != Name2:
                    MutRep.append(SimDif)
    return MutRep



def MutMap(Sequence):

#     goal to identify RGYW, WRCY, WA, TW, coldspots, and Pol-eta hs and cs, plus C-to-T targeting
# to optimally remove targetting under constraints of using codons at physiological frequencies
# can decinstruct a sequence and make every possible variant then choose least and most evolved
# can use to look at effects of mutation on changing substrate
# for aged study, what are the consequences of relying on memory cells to much?
    features = []
    Sequence = Sequence.upper()
    i=0
    j=1
    scores = ()
    for base in range(0, len(Sequence)-2):

        MutType = ''
        CTType = ''
        RScore = ''

        if j == 1:
            codon = Sequence[base:base+3]
            AA = Translator(codon, 0)
        Nuc = Sequence[base]
        if Nuc == 'C' and base >1 and base < len(Sequence)-1:
            if (Sequence[base-2] == 'A' or Sequence[base-2] == 'T') and (Sequence[base-1] == 'A' or Sequence[base-1] == 'G') and (Sequence[base+1] == 'C' or Sequence[base+1] == 'T'):
                MutType = 'WRCY'
            elif (Sequence[base-2] == 'A' or Sequence[base-2] == 'T') and (Sequence[base-1] == 'A' or Sequence[base-1] == 'G')and not (Sequence[base+1] == 'C' or Sequence[base+1] == 'T'):
                MutType = 'WRC'
            elif ((Sequence[base-2] == 'G' or Sequence[base-2] == 'C') and (Sequence[base-1] == 'C' or Sequence[base-1] == 'T')) or Sequence[base-2:base] == 'TTC' or Sequence[base-2:base] == 'CAC' or Sequence[base-2:base] == 'GGC' or Sequence[base-2:base] == 'GAC':
                MutType = 'AID-CS'
            else:
                MutType = "Neutral"
        if Nuc == 'C':
            if j == 1:
                CkCodon = 'T' + codon[1:]
            elif j ==2:
                try:
                    CkCodon = codon[0] + 'T' + codon[2]
                except:
                    print('stop')
            elif j == 3:
                CkCodon = codon[0:1] + 'T'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'G' and base >0 and base < len(Sequence)-2:
            if (Sequence[base-1] == 'G' or Sequence[base-1] == 'A') and (Sequence[base+1] == 'C' or Sequence[base+1] == 'T') and (Sequence[base+2] == 'A' or Sequence[base+2] == 'T'):
                MutType = 'WRCY'
            elif (Sequence[base+1] == 'C' or Sequence[base+1] == 'T') and (Sequence[base+2] == 'A' or Sequence[base+2] == 'T'):
                MutType = 'WRC'
            elif ((Sequence[base+1] == 'G' or Sequence[base+1] == 'A') and (Sequence[base+2] == 'C' or Sequence[base+2] == 'G')) or Sequence[base:base+2] == 'GAA' or Sequence[base:base+2] == 'GTG' or Sequence[base:base+2] == 'GCC' or Sequence[base:base+2] == 'GTC':
                MutType = 'AID-CS'
            else:
                MutType = "Neutral"
        if Nuc == 'G':
            if j == 1:
                CkCodon = 'A' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'A' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'A'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'A'and base >0:  #WA
            if (Sequence[base-1] == 'A' or Sequence[base-1] == 'T'):
                MutType = 'WA'
            else:
                MutType = 'Neutral'

        if Nuc == 'A':
            if j == 1:
                CkCodon = 'G' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'G' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'G'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'T'and base < len(Sequence)-1:  #WA
            if (Sequence[base+1] == 'A' or Sequence[base+1] == 'T'):
                MutType = 'TW'
            else:
                MutType = 'Neutral'

        if Nuc == 'T':
            if j == 1:
                CkCodon = 'C' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'C' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'C'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        NMut = Nuc+MutType
        TENDA, TENDG, TENDC, TENDT, MutIndex = 0,0,0,0,0

        if NMut in Mutability:
            MutIndex = Mutability[NMut]

        if NMut == 'CWRCY':
            TENDA = MutTendency['CAWRCY']
            TENDG = MutTendency['CGWRCY']
            TENDT = MutTendency['CTWRCY']
            TENDC = 0
        elif NMut == 'CWRC':
            TENDA = MutTendency['CAWRC']
            TENDG = MutTendency['CGWRC']
            TENDT = MutTendency['CTWRC']
            TENDC = 0
        elif NMut == 'CAID-CS':
            TENDA = MutTendency['CAAID-CS']
            TENDG = MutTendency['CGAID-CS']
            TENDT = MutTendency['CTAID-CS']
            TENDC = 0
        elif NMut == 'CNeutral':
            TENDA = MutTendency['CANeutral']
            TENDG = MutTendency['CGNeutral']
            TENDT = MutTendency['CTNeutral']
            TENDC = 0
        elif NMut == 'GWRCY':
            TENDA = MutTendency['GAWRCY']
            TENDG = 0
            TENDT = MutTendency['GTWRCY']
            TENDC = MutTendency['GCWRCY']
        elif NMut == 'GWRC':
            TENDA = MutTendency['GAWRC']
            TENDG = 0
            TENDT = MutTendency['GTWRC']
            TENDC = MutTendency['GCWRC']
        elif NMut == 'GAID-CS':
            TENDA = MutTendency['GAAID-CS']
            TENDG = 0
            TENDT = MutTendency['GTAID-CS']
            TENDC = MutTendency['GCAID-CS']
        elif NMut == 'GWRC':
            TENDA = MutTendency['GANeutral']
            TENDG = 0
            TENDT = MutTendency['GTNeutral']
            TENDC = MutTendency['GCNeutral']
        elif NMut == 'AWA':
            TENDA = 0
            TENDG = MutTendency['AGWA']
            TENDT = MutTendency['ATWA']
            TENDC = MutTendency['ACWA']
        elif NMut == 'ANeutral':
            TENDA = 0
            TENDG = MutTendency['AGNeutral']
            TENDT = MutTendency['ATNeutral']
            TENDC = MutTendency['ACNeutral']
        elif NMut == 'TTW':
            TENDA = MutTendency['TATW']
            TENDG = MutTendency['TGTW']
            TENDT = 0
            TENDC = MutTendency['TCTW']
        elif NMut == 'TNeutral':
            TENDA = MutTendency['TANeutral']
            TENDG = MutTendency['TGNeutral']
            TENDT = 0
            TENDC = MutTendency['TCNeutral']

        scores = (Nuc, MutType, CTType, RScore, MutIndex, TENDA, TENDG, TENDC, TENDT)
        #Nuc, obvious
        # MutType, whether in a hot, cold, or neutral for AID or Pol-eta
        # CTType, whether if C to T or G to A mutation would be silent (not verified code yet)
        # RScore, when scoring mutations oif it is a conservative amino acid change based on Blosum62
        # MutIndex, this is the empirical frequency of mutations determined from analyzing a million sequences if in hotspot, coldspot, etc
        # TENDA, similar to MutIndex but it's tendency for the current base to become an A based on empirical
        # TENDG, similar to MutIndex but it's tendency for the current base to become a G based on empirical
        # TENDC, similar to MutIndex but it's tendency for the current base to become a C based on empirical
        # TENDT similar to MutIndex but it's tendency for the current base to become a T based on empirical

        features.append(scores)
        if j <3:
            j+=1
        else:
            j=1
        # i+=1

    return features

def GCMutMap(Sequence):

#     goal to identify RGYW, WRCY, WA, TW, coldspots, and Pol-eta hs and cs, plus C-to-T targeting
# to optimally remove targetting under constraints of using codons at physiological frequencies
# can decinstruct a sequence and make every possible variant then choose least and most evolved
# can use to look at effects of mutation on changing substrate
# for aged study, what are the consequences of relying on memory cells to much?
    features = []
    Sequence = Sequence.upper()
    i=0
    j=1
    scores = ()
    for base in range(0, len(Sequence)-2):

        MutType = ''
        CTType = ''
        RScore = ''

        if j == 1:
            codon = Sequence[base:base+3]
            AA = Translator(codon, 0)
        Nuc = Sequence[base]
        if Nuc == 'C' and base >1 and base < len(Sequence)-1:
            if (Sequence[base-2] == 'A' or Sequence[base-2] == 'T') and (Sequence[base-1] == 'A' or Sequence[base-1] == 'G') and (Sequence[base+1] == 'C' or Sequence[base+1] == 'T'):
                MutType = 'WRCY'
            elif (Sequence[base-2] == 'A' or Sequence[base-2] == 'T') and (Sequence[base-1] == 'A' or Sequence[base-1] == 'G')and not (Sequence[base+1] == 'C' or Sequence[base+1] == 'T'):
                MutType = 'WRC'
            elif ((Sequence[base-2] == 'G' or Sequence[base-2] == 'C') and (Sequence[base-1] == 'C' or Sequence[base-1] == 'T')) or Sequence[base-2:base] == 'TTC' or Sequence[base-2:base] == 'CAC' or Sequence[base-2:base] == 'GGC' or Sequence[base-2:base] == 'GAC':
                MutType = 'AID-CS'
            else:
                MutType = "Neutral"
        if Nuc == 'C':
            if j == 1:
                CkCodon = 'T' + codon[1:]
            elif j ==2:
                try:
                    CkCodon = codon[0] + 'T' + codon[2]
                except:
                    print('stop')
            elif j == 3:
                CkCodon = codon[0:1] + 'T'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'G' and base >0 and base < len(Sequence)-2:
            if (Sequence[base-1] == 'G' or Sequence[base-1] == 'A') and (Sequence[base+1] == 'C' or Sequence[base+1] == 'T') and (Sequence[base+2] == 'A' or Sequence[base+2] == 'T'):
                MutType = 'WRCY'
            elif (Sequence[base+1] == 'C' or Sequence[base+1] == 'T') and (Sequence[base+2] == 'A' or Sequence[base+2] == 'T'):
                MutType = 'WRC'
            elif ((Sequence[base+1] == 'G' or Sequence[base+1] == 'A') and (Sequence[base+2] == 'C' or Sequence[base+2] == 'G')) or Sequence[base:base+2] == 'GAA' or Sequence[base:base+2] == 'GTG' or Sequence[base:base+2] == 'GCC' or Sequence[base:base+2] == 'GTC':
                MutType = 'AID-CS'
            else:
                MutType = "Neutral"
        if Nuc == 'G':
            if j == 1:
                CkCodon = 'A' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'A' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'A'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'A'and base >0:  #WA
            if (Sequence[base-1] == 'A' or Sequence[base-1] == 'T'):
                MutType = 'WA'
            else:
                MutType = 'Neutral'

        if Nuc == 'A':
            if j == 1:
                CkCodon = 'G' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'G' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'G'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        if Nuc == 'T'and base < len(Sequence)-1:  #WA
            if (Sequence[base+1] == 'A' or Sequence[base+1] == 'T'):
                MutType = 'TW'
            else:
                MutType = 'Neutral'

        if Nuc == 'T':
            if j == 1:
                CkCodon = 'C' + codon[1:]
            elif j ==2:
                CkCodon = codon[0] + 'C' + codon[2]
            elif j == 3:
                CkCodon = codon[0:1] + 'C'
            CkAA = Translator(CkCodon, 0)
            if AA[0] != CkAA[0]:
                if CkAA[0] != '*':
                    CTType = 'missense'
                    RScore = AASimilarity(AA[0]+CkAA[0])  #gets the BLOSUM62 score for conservation liklihood, negative is less, + is likly
                else:
                    CTType = 'nonsense'
            else:
                CTType = 'silent'

        NMut = Nuc+MutType
        TENDA, TENDG, TENDC, TENDT, MutIndex = 0,0,0,0,0

        if NMut in GCMutability:
            MutIndex = GCMutability[NMut]

        if NMut == 'CWRCY':
            TENDA = MutTendency['CAWRCY']
            TENDG = MutTendency['CGWRCY']
            TENDT = MutTendency['CTWRCY']
            TENDC = 0
        elif NMut == 'CWRC':
            TENDA = MutTendency['CAWRC']
            TENDG = MutTendency['CGWRC']
            TENDT = MutTendency['CTWRC']
            TENDC = 0
        elif NMut == 'CAID-CS':
            TENDA = MutTendency['CAAID-CS']
            TENDG = MutTendency['CGAID-CS']
            TENDT = MutTendency['CTAID-CS']
            TENDC = 0
        elif NMut == 'CNeutral':
            TENDA = MutTendency['CANeutral']
            TENDG = MutTendency['CGNeutral']
            TENDT = MutTendency['CTNeutral']
            TENDC = 0
        elif NMut == 'GWRCY':
            TENDA = MutTendency['GAWRCY']
            TENDG = 0
            TENDT = MutTendency['GTWRCY']
            TENDC = MutTendency['GCWRCY']
        elif NMut == 'GWRC':
            TENDA = MutTendency['GAWRC']
            TENDG = 0
            TENDT = MutTendency['GTWRC']
            TENDC = MutTendency['GCWRC']
        elif NMut == 'GAID-CS':
            TENDA = MutTendency['GAAID-CS']
            TENDG = 0
            TENDT = MutTendency['GTAID-CS']
            TENDC = MutTendency['GCAID-CS']
        elif NMut == 'GWRC':
            TENDA = MutTendency['GANeutral']
            TENDG = 0
            TENDT = MutTendency['GTNeutral']
            TENDC = MutTendency['GCNeutral']
        elif NMut == 'AWA':
            TENDA = 0
            TENDG = MutTendency['AGWA']
            TENDT = MutTendency['ATWA']
            TENDC = MutTendency['ACWA']
        elif NMut == 'ANeutral':
            TENDA = 0
            TENDG = MutTendency['AGNeutral']
            TENDT = MutTendency['ATNeutral']
            TENDC = MutTendency['ACNeutral']
        elif NMut == 'TTW':
            TENDA = MutTendency['TATW']
            TENDG = MutTendency['TGTW']
            TENDT = 0
            TENDC = MutTendency['TCTW']
        elif NMut == 'TNeutral':
            TENDA = MutTendency['TANeutral']
            TENDG = MutTendency['TGNeutral']
            TENDT = 0
            TENDC = MutTendency['TCNeutral']

        scores = (Nuc, MutType, CTType, RScore, MutIndex, TENDA, TENDG, TENDC, TENDT)
        #Nuc, obvious
        # MutType, whether in a hot, cold, or neutral for AID or Pol-eta
        # CTType, whether if C to T or G to A mutation would be silent (not verified code yet)
        # RScore, when scoring mutations oif it is a conservative amino acid change based on Blosum62
        # MutIndex, this is the empirical frequency of mutations determined from analyzing a million sequences if in hotspot, coldspot, etc
        # TENDA, similar to MutIndex but it's tendency for the current base to become an A based on empirical
        # TENDG, similar to MutIndex but it's tendency for the current base to become a G based on empirical
        # TENDC, similar to MutIndex but it's tendency for the current base to become a C based on empirical
        # TENDT similar to MutIndex but it's tendency for the current base to become a T based on empirical

        features.append(scores)
        if j <3:
            j+=1
        else:
            j=1
        # i+=1

    return features

Mutability = {'CWRCY':5.328476052, 'CWRC':1.615976588, 'CAID-CS':0.878356004, 'CNeutral':1.080514212,
                'GWRCY':3.75544955, 'GWRC':4.626183008, 'GAID-CS':0.725579452, 'GNeutral':1.512267645,
                'AWA':3.897581991, 'ANeutral':2.518229922, 'TTW':2.646667088, 'TNeutral':0.889080601}
MutTendency = {'CTWRCY':3.419071313, 'CTWRC':0.870149852, 'CTAID-CS':0.350768348, 'CTNeutral':0.72123155,
               'CGWRCY':0.022854755, 'CGWRC':0.031806369, 'CGAID-CS':0.201546946, 'CGNeutral':0.225848265,
               'CAWRCY':0.00152365, 'CAWRC':0.048097436, 'CAAID-CS':0.152787026, 'CANeutral':0.157784404,
                'GTWRCY':0.092268681, 'GTWRC':2.371327585, 'GTAID-CS':0.145110857, 'GTNeutral':0.188034032,
                'GCWRCY':0.241181924, 'GCWRC':0.573265113, 'GCAID-CS':0.073120796, 'GCNeutral':0.178289758,
                'GAWRCY':4.103190461, 'GAWRC':3.620687502, 'GAAID-CS':0.216648625, 'GANeutral':0.801415762,
                'ATWA':0.730736321, 'ATNeutral':0.273391661, 'TATW':0.275989829, 'TANeutral':0.163184205,
                'AGWA':2.383833337, 'AGNeutral':1.502258807, 'TGTW':0.425101597, 'TGNeutral':0.194067827,
                'ACWA':1.090225932, 'ACNeutral':0.608828997, 'TCTW':1.549320198, 'TCNeutral':0.409848594}
GCMutability = {'CWRCY':5.328476052, 'CWRC':1.615976588, 'CAID-CS':0.878356004, 'CNeutral':1.080514212,
                'GWRCY':3.75544955, 'GWRC':4.626183008, 'GAID-CS':0.725579452, 'GNeutral':1.512267645,
                'AWA':1, 'ANeutral':1, 'TTW':1, 'TNeutral':1}

def AASimilarity(AA):
    # Proc Natl Acad Sci U S A. 1992 Nov 15;89(22):10915-9.
    # Amino acid substitution matrices from protein blocks.
    # Henikoff S1, Henikoff JG.
    # Based on above reference and the BLOSUM62 similarity indices
    # BLOSUM62 = {}
    ExScore = ''
    try:
        if AA in BLOSUM62:
            ExScore = BLOSUM62[AA]
        else:
            AA = AA[1]+ AA[0]
            ExScore = BLOSUM62[AA]
    except:
            ExScore = ''
    return ExScore

BLOSUM62 = {'AA':4,'AR':-1, 'AN':-2, 'AD':-2, 'AC':0,'AQ':-1, 'AE':-1,
'AG':0, 'AH':-2,'AI':-1, 'AL':-1, 'AK':-1, 'AM':-1,'AF':-2,
'AP':-1, 'AS':1,'AT':0,'AW':-3, 'AY':-2, 'AV':0,'RR':4,
'RN':0, 'RD':-2, 'RC':-3,'RQ':1, 'RE':0,'RG':-2, 'RH':0,
'RI':-3, 'RL':-2, 'RK':2, 'RM':-1,'RF':-3,'RP':-2, 'RS':-1,
'RT':-1,'RW':-3, 'RY':-2, 'RV':-3,'NN':6, 'ND':1, 'NC':-3,
'NQ':0, 'NE':0,'NG':0, 'NH':1,'NI':-3, 'NL':-3, 'NK':0, 'NM':-2,
'NF':-3,'NP':-2, 'NS':1,'NT':0,'NW':-4, 'NY':-2, 'NV':-3,'DD':6,
'DC':-3,'DQ':0, 'DE':2,'DG':-1, 'DH':-1,'DI':-3, 'DL':-4, 'DK':-1,
'DM':-3,'DF':-3,'DP':-1, 'DS':0,'DT':-1,'DW':-4, 'DY':-3, 'DV':-3,
'CC':9,'CQ':-3, 'CE':-4,'CG':-3, 'CH':-3,'CI':-1, 'CL':-1, 'CK':-3,
'CM':-1,'CF':-2,'CP':-3, 'CS':-1,'CT':-1,'CW':-2, 'CY':-2, 'CV':-1,
'QQ':5, 'QE':2,'QG':-2, 'QH':0,'QI':-3, 'QL':-2, 'QK':1, 'QM':0,
'QF':-3,'QP':-1, 'QS':0,'QT':-1,'QW':-2, 'QY':-1, 'QV':-2,'EE':5,
'EG':-2, 'EH':0,'EI':-3, 'EL':-3, 'EK':1, 'EM':-2,'EF':-3,'EP':-1,
'ES':0,'ET':-1,'EW':-3, 'EY':-2, 'EV':-2,'GG':6, 'GH':-2,'GI':-4, 'GL':-4,
'GK':-2, 'GM':-3,'GF':-3,'GP':-2, 'GS':0,'GT':-2,'GW':-2, 'GY':-3, 'GV':-3,
'HH':8,'HI':-3, 'HL':-3, 'HK':-1, 'HM':-2,'HF':-1,'HP':-2, 'HS':-1,'HT':-2,
'HW':-2, 'HY':2, 'HV':-3,'II':4, 'IL':2, 'IK':-3, 'IM':1,'IF':0,'IP':-3,
'IS':-2,'IT':-1,'IW':-3, 'IY':-1, 'IV':3,'LL':4, 'LK':-2, 'LM':2,'LF':0,
'LP':-3, 'LS':-2,'LT':-1,'LW':-2, 'LY':-1, 'LV':1,'KK':5, 'KM':-1,'KF':-3,
'KP':-1, 'KS':0,'KT':-1,'KW':-3, 'KY':-2, 'KV':-2,'MM':5,'MF':0,'MP':-2,
'MS':-1,'MT':-1,'MW':-1, 'MY':-1, 'MV':1,'FF':6,'FP':-4, 'FS':-2,'FT':-2,
'FW':1, 'FY':3, 'FV':-1,'PP':7, 'PS':-1,'PT':-1,'PW':-4, 'PY':-3, 'PV':-2,
'SS':4,'ST':1,'SW':-3, 'SY':-2, 'SV':-2,'TT':5,'TW':-2, 'TY':-2, 'TV':0,
'WW':11, 'WY':2, 'WV':3,'YY':7, 'YV':-1, 'VV':5}

# template sequences for this function is downloaded from IMGT
# http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=human&latin=Homo%20sapiens&group=IGHC
def CallIsotype(IsoSeq):
    # import sys
    # sys.path.append('/Users/PCW-MacBookProRet/anaconda/lib/python3.5/site-packages/editdistance/__init__.py')
    # import editdistance
    # SetPath = os.path.join(os.path.expanduser('~'), 'anaconda', 'lib', 'python3.5', 'site-packages', 'editdistance', '__init__.py')
    # SetPath = os.path.join(os.path.expanduser('~'), 'anaconda', 'lib', 'python3.5', 'site-packages', 'editdistance')
    # from runpy import run_path
    # settings = run_path(SetPath)
    # os.chdir(SetPath)
    # import pkgutil
    # import importlib
    # import editdistance
    # packages = pkgutil.walk_packages(path='.')
    # for importer, name, is_package in packages:
    #     mod = importlib.import_module(name)
    #
    # SetPath = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Output.txt')
    # Distance  = 0
    # with open(SetPath, 'w') as currentFile:  #using with for this automatically closes the file even if you crash
    #     currentFile.write(str(Distance))

    # input is a list of lists containing filename and sequence: ((fielname1, seq1),(fielname2, seq2))
    import difflib

    IsoSeq = IsoSeq.replace('-', '')

    i = 0
    Ns = []
    SeqLen = len(IsoSeq)
    IsoSeq = IsoSeq.upper()

    for base in IsoSeq:
        i += 1
        if base == 'N':
            Count = (i,1)
            Ns.append(Count)
        else:
            Count = (i, 0)
            Ns.append(Count)
        if i > 9:
            Ncount = 0
            for j in range((i-10), i):
                Ncount += Ns[j][1]
            if Ncount >1:
                SeqLen = i
                IsoSeq = IsoSeq[:SeqLen]
                break


    # if SeqLen > 63:
    # IsoSeq = IsoSeq[:63]
    #     SeqLen = 63


    CS = len(IsoSeq)
    if CS < 5:
        return 'Unknown'

    ToClustalO = []
    # Seq = ('TestSeq', IsoSeq)
    # ToClustalO.append(Seq)


    IgMSeq = 'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGGCCGTTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTGTCCTGGAAATACAAGAACAACTCTGACATCAGCAGTACCCGGGGCTTCCCATCAGTCCTGAGAGGGGGCAAGTACGCAGCCACCTCACAGGTGCTGCTGCCTTCCAAGGACGTCATGCAGGGCACAGACGAACACGTGGTGTGCAAAGTCCAGCACCCCAACGGCAACAAAGAAAAGAACGTGCCTCTTCCAG'
    # IgMSeq = 'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGG'

    IgMSeq = IgMSeq[:SeqLen]
    Seq = ('IgM', IgMSeq)

    ToClustalO.append(Seq)



    IgG1Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGCGTGGTGGTGGACGTGAGCCACGAAGACCCTGAGGTCAAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTACCGGGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAATGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGATGAGCTGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTATCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
    # IgG1Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCC'
    IgG1Seq = IgG1Seq[:SeqLen]
    Seq = ('IgG1', IgG1Seq)
    ToClustalO.append(Seq)


    IgG2Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCAGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGAGCGCAAATGTTGTGTCGAGTGCCCACCGTGCCCAGCACCACCTGTGGCAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCACGGGAGGAGCAGTTCAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTTGTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACACCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
    # IgG2Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGG'
    IgG2Seq = IgG2Seq[:SeqLen]
    Seq = ('IgG2', IgG2Seq)
    ToClustalO.append(Seq)



    IgG3Seq = 'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACACCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGCTCAAAACCCCACTTGGTGACACAACTCACACATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCAAGGTGCCCAGCACCTGAACTCCTGGGAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGATACCCTTATGATTTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAAGTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGACAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAGCGGGCAGCCGGAGAACAACTACAACACCACGCCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACATCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCGCTTCACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
    # IgG3Seq = 'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGG'
    IgG3Seq = IgG3Seq[:SeqLen]
    Seq = ('IgG3', IgG3Seq)
    ToClustalO.append(Seq)



    IgG4Seq = 'CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACGAAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGTCCAAATATGGTCCCCCATGCCCATCATGCCCAGCACCTGAGTTCCTGGGGGGACCATCAGTCTTCCTGTTCCCCCCAAAACCCAAGGACACTCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCAGGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGATGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTTCAACAGCACGTACCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCGTCCTCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAGCCACAGGTGTACACCCTGCCCCCATCCCAGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAGGCTAACCGTGGACAAGAGCAGGTGGCAGGAGGGGAATGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACACAGAAGAGCCTCTCCCTGTCTCTGGGTAAATGA'
    # IgG4Seq = 'CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGG'
    IgG4Seq = IgG4Seq[:SeqLen]
    Seq = ('IgG4', IgG4Seq)
    ToClustalO.append(Seq)


    IgA1seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGGGCGTGACCGCCAGAAACTTCCCACCCAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCTAGCCGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCTCAACTCCACCTACCCCATCTCCCTCAACTCCACCTACCCCATCTCCCTCATGCTGCCACCCCCGACTGTCACTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCAGGTGTCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCGGGCTGTGCCGAGCCATGGAACCATGGGAAGACCTTCACTTGCACTGCTGCCTACCCCGAGTCCAAGACCCCGCTAACCGCCACCCTCTCAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGCGGCTTCAGCCCCAAGGACGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'
    # IgA1seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTG'
    IgA1seq = IgA1seq[:SeqLen]
    Seq = ('IgA1', IgA1seq)
    ToClustalO.append(Seq)



    IgA2seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGAACGTGACCGCCAGAAACTTCCCACCTAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCCAGACGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCCCACCTCCCCCATGCTGCCACCCCCGACTGTCGCTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCTGGTGCCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCTGGCTGTGCCCAGCCATGGAACCATGGGGAGACCTTCACCTGCACTGCTGCCCACCCCGAGTTGAAGACCCCACTAACCGCCAACATCACAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGTGGCTTCAGCCCCAAGGATGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'
    # IgA2seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGCCTG'
    IgA2seq = IgA2seq[:SeqLen]
    Seq = ('IgA2', IgA2seq)
    ToClustalO.append(Seq)



    IgEseq = 'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCTGCAAAAACATTCCCTCCAATGCCACCTCCGTGACTCTGGGCTGCCTGGCCACGGGCTACTTCCCGGAGCCGGTGATGGTGACCTGCGACACAGGCTCCCTCAACGGGACAACTATGACCTTACCAGCCACCACCCTCACGCTCTCTGGTCACTATGCCACCATCAGCTTGCTGACCGTCTCGGGTGCGTGGGCCAAGCAGATGTTCACCTGCCGTGTGGCACACACTCCATCGTCCACAGACTGGGTCGACAACAAAACCTTCAGCGTCTGCTCCAGGGACTTCACCCCGCCCACCGTGAAGATCTTACAGTCGTCCTGCGACGGCGGCGGGCACTTCCCCCCGACCATCCAGCTCCTGTGCCTCGTCTCTGGGTACACCCCAGGGACTATCAACATCACCTGGCTGGAGGACGGGCAGGTCATGGACGTGGACTTGTCCACCGCCTCTACCACGCAGGAGGGTGAGCTGGCCTCCACACAAAGCGAGCTCACCCTCAGCCAGAAGCACTGGCTGTCAGACCGCACCTACACCTGCCAGGTCACCTATCAAGGTCACACCTTTGAGGACAGCACCAAGAAGTGTGCAGATTCCAACCCGAGAGGGGTGAGCGCCTACCTAAGCCGGCCCAGCCCGTTCGACCTGTTCATCCGCAAGTCGCCCACGATCACCTGTCTGGTGGTGGACCTGGCACCCAGCAAGGGGACCGTGAACCTGACCTGGTCCCGGGCCAGTGGGAAGCCTGTGAACCACTCCACCAGAAAGGAGGAGAAGCAGCGCAATGGCACGTTAACCGTCACGTCCACCCTGCCGGTGGGCACCCGAGACTGGATCGAGGGGGAGACCTACCAGTGCAGGGTGACCCACCCCCACCTGCCCAGGGCCCTCATGCGGTCCACGACCAAGACCAGCGGCCCGCGTGCTGCCCCGGAAGTCTATGCGTTTGCGACGCCGGAGTGGCCGGGGAGCCGGGACAAGCGCACCCTCGCCTGCCTGATCCAGAACTTCATGCCTGAGGACATCTCGGTGCAGTGGCTGCACAACGAGGTGCAGCTCCCGGACGCCCGGCACAGCACGACGCAGCCCCGCAAGACCAAGGGCTCCGGCTTCTTCGTCTTCAGCCGCCTGGAGGTGACCAGGGCCGAATGGGAGCAGAAAGATGAGTTCATCTGCCGTGCAGTCCATGAGGCAGCGAGCCCCTCACAGACCGTCCAGCGAGCGGTGTCTGTAAATCCCGGTAAATGA'
    # IgEseq = 'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCTGCAAAAACATTCCCTCCAATGCCACCTCCGTGACTCTGGG'
    IgEseq = IgEseq[:SeqLen]
    # IgEseq = IgEseq[:30]
    Seq = ('IgE', IgEseq)
    ToClustalO.append(Seq)



    IgDseq = 'CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGGCATGCTTGATAACTGGGTACCACCCAACGTCCGTGACTGTCACCTGGTACATGGGGACACAGAGCCAGCCCCAGAGAACCTTCCCTGAGATACAAAGACGGGACAGCTACTACATGACAAGCAGCCAGCTCTCCACCCCCCTCCAGCAGTGGCGCCAAGGCGAGTACAAATGCGTGGTCCAGCACACCGCCAGCAAGAGTAAGAAGGAGATCTTCCGCTGGCCAGAGTCTCCAAAGGCACAGGCCTCCTCCGTGCCCACTGCACAACCCCAAGCAGAGGGCAGCCTCGCCAAGGCAACCACAGCCCCAGCCACCACCCGTAACACAGGAAGAGGAGGAGAAGAGAAGAAGAAGGAGAAGGAGAAAGAGGAACAAGAAGAGAGAGAGACAAAGACACCAGAGTGTCCGAGCCACACCCAGCCTCTTGGCGTCTACCTGCTAACCCCTGCAGTGCAGGACCTGTGGCTCCGGGACAAAGCCACCTTCACCTGCTTCGTGGTGGGCAGTGACCTGAAGGATGCTCACCTGACCTGGGAGGTGGCTGGGAAGGTCCCCACAGGGGGCGTGGAGGAAGGGCTGCTGGAGCGGCACAGCAACGGCTCCCAGAGCCAGCACAGCCGTCTGACCCTGCCCAGGTCCTTGTGGAACGCGGGGACCTCCGTCACCTGCACACTGAACCATCCCAGCCTCCCACCCCAGAGGTTGATGGCGCTGAGAGAACCCGCTGCGCAGGCACCCGTCAAGCTTTCTCTGAACCTGCTGGCCTCGTCTGACCCTCCCGAGGCGGCCTCGTGGCTCCTGTGTGAGGTGTCTGGCTTCTCGCCCCCCAACATCCTCCTGATGTGGCTGGAGGACCAGCGTGAGGTGAACACTTCTGGGTTTGCCCCCGCACGCCCCCCTCCACAGCCCAGGAGCACCACGTTCTGGGCCTGGAGTGTGCTGCGTGTCCCAGCCCCGCCCAGCCCTCAGCCAGCCACCTACACGTGTGTGGTCAGCCACGAGGACTCCCGGACTCTGCTCAACGCCAGCCGGAGCCTAGAAGTCAGCTACCTGGCCATGACCCCCCTGATCCCTCAGAGCAAGGATGAGAACAGCGATGACTACACGACCTTTGATGATGTGGGCAGCCTGTGGACCACCCTGTCCACGTTTGTGGCCCTCTTCATCCTCACCCTCCTCTACAGCGGCATTGTCACTTTCATCAAGGTGAAGTAG'
    # IgDseq = 'CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGG'
    IgDseq = IgDseq[:SeqLen]
    Seq = ('IgD', IgDseq)
    ToClustalO.append(Seq)

    ABVectorSeq = 'CGTCGACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCTGTGACGGTCTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGC'
    ABVectorSeq = ABVectorSeq[:SeqLen]
    Seq = ('Ab-vector', ABVectorSeq)
    ToClustalO.append(Seq)


    hit = "Unknown"

    IgSeqs = [IgMSeq,IgA1seq,IgA2seq,IgDseq,IgEseq,IgG1Seq,IgG2Seq,IgG3Seq,IgG4Seq,ABVectorSeq]

    best = difflib.get_close_matches(IsoSeq, IgSeqs, 1, 0.6)


    if len(best)>0:
        for item in ToClustalO:
            # test = item[1]
            if item[1] == best[0]:
                hit = item[0]

                if hit == 'IgG1' or hit == 'IgG2':
                    if SeqLen < 35:
                        hit = 'IgG'
                if hit == 'IgG3' or hit == 'IgG4':
                    if SeqLen < 20:
                        hit = 'IgG'
                if hit == 'IgA1' or hit == 'IgA2':
                    if SeqLen < 39:
                        hit = 'IgA'
                break
    else:
        if SeqLen >39:
            SeqLen = 40
        elif SeqLen >19:
            SeqLen = 20
        else:
            return 'Unknown'

        IsoSeq = IsoSeq[:SeqLen]

        CS = len(IsoSeq)
        if CS < 5:
            return 'Unknown'

        ToClustalO = []
        # Seq = ('TestSeq', IsoSeq)
        # ToClustalO.append(Seq)


        IgMSeq = 'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGGCCGTTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTGTCCTGGAAATACAAGAACAACTCTGACATCAGCAGTACCCGGGGCTTCCCATCAGTCCTGAGAGGGGGCAAGTACGCAGCCACCTCACAGGTGCTGCTGCCTTCCAAGGACGTCATGCAGGGCACAGACGAACACGTGGTGTGCAAAGTCCAGCACCCCAACGGCAACAAAGAAAAGAACGTGCCTCTTCCAG'
        # IgMSeq = 'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGG'

        IgMSeq = IgMSeq[:SeqLen]
        Seq = ('IgM', IgMSeq)
        ToClustalO.append(Seq)

        IgG1Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGCGTGGTGGTGGACGTGAGCCACGAAGACCCTGAGGTCAAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTACCGGGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAATGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGATGAGCTGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTATCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
        # IgG1Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCC'
        IgG1Seq = IgG1Seq[:SeqLen]
        Seq = ('IgG1', IgG1Seq)
        ToClustalO.append(Seq)

        IgG2Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCAGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGAGCGCAAATGTTGTGTCGAGTGCCCACCGTGCCCAGCACCACCTGTGGCAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCACGGGAGGAGCAGTTCAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTTGTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACACCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
        # IgG2Seq = 'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGG'
        IgG2Seq = IgG2Seq[:SeqLen]
        Seq = ('IgG2', IgG2Seq)
        ToClustalO.append(Seq)

        IgG3Seq = 'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACACCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGCTCAAAACCCCACTTGGTGACACAACTCACACATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCAAGGTGCCCAGCACCTGAACTCCTGGGAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGATACCCTTATGATTTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAAGTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGACAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAGCGGGCAGCCGGAGAACAACTACAACACCACGCCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACATCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCGCTTCACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'
        # IgG3Seq = 'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGG'
        IgG3Seq = IgG3Seq[:SeqLen]
        Seq = ('IgG3', IgG3Seq)
        ToClustalO.append(Seq)

        IgG4Seq = 'CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACGAAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGTCCAAATATGGTCCCCCATGCCCATCATGCCCAGCACCTGAGTTCCTGGGGGGACCATCAGTCTTCCTGTTCCCCCCAAAACCCAAGGACACTCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCAGGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGATGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTTCAACAGCACGTACCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCGTCCTCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAGCCACAGGTGTACACCCTGCCCCCATCCCAGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAGGCTAACCGTGGACAAGAGCAGGTGGCAGGAGGGGAATGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACACAGAAGAGCCTCTCCCTGTCTCTGGGTAAATGA'
        # IgG4Seq = 'CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGG'
        IgG4Seq = IgG4Seq[:SeqLen]
        Seq = ('IgG4', IgG4Seq)
        ToClustalO.append(Seq)

        IgA1seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGGGCGTGACCGCCAGAAACTTCCCACCCAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCTAGCCGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCTCAACTCCACCTACCCCATCTCCCTCAACTCCACCTACCCCATCTCCCTCATGCTGCCACCCCCGACTGTCACTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCAGGTGTCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCGGGCTGTGCCGAGCCATGGAACCATGGGAAGACCTTCACTTGCACTGCTGCCTACCCCGAGTCCAAGACCCCGCTAACCGCCACCCTCTCAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGCGGCTTCAGCCCCAAGGACGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'
        # IgA1seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTG'
        IgA1seq = IgA1seq[:SeqLen]
        Seq = ('IgA1', IgA1seq)
        ToClustalO.append(Seq)

        IgA2seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGAACGTGACCGCCAGAAACTTCCCACCTAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCCAGACGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCCCACCTCCCCCATGCTGCCACCCCCGACTGTCGCTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCTGGTGCCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCTGGCTGTGCCCAGCCATGGAACCATGGGGAGACCTTCACCTGCACTGCTGCCCACCCCGAGTTGAAGACCCCACTAACCGCCAACATCACAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGTGGCTTCAGCCCCAAGGATGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'
        # IgA2seq = 'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGCCTG'
        IgA2seq = IgA2seq[:SeqLen]
        Seq = ('IgA2', IgA2seq)
        ToClustalO.append(Seq)

        IgEseq = 'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCTGCAAAAACATTCCCTCCAATGCCACCTCCGTGACTCTGGGCTGCCTGGCCACGGGCTACTTCCCGGAGCCGGTGATGGTGACCTGCGACACAGGCTCCCTCAACGGGACAACTATGACCTTACCAGCCACCACCCTCACGCTCTCTGGTCACTATGCCACCATCAGCTTGCTGACCGTCTCGGGTGCGTGGGCCAAGCAGATGTTCACCTGCCGTGTGGCACACACTCCATCGTCCACAGACTGGGTCGACAACAAAACCTTCAGCGTCTGCTCCAGGGACTTCACCCCGCCCACCGTGAAGATCTTACAGTCGTCCTGCGACGGCGGCGGGCACTTCCCCCCGACCATCCAGCTCCTGTGCCTCGTCTCTGGGTACACCCCAGGGACTATCAACATCACCTGGCTGGAGGACGGGCAGGTCATGGACGTGGACTTGTCCACCGCCTCTACCACGCAGGAGGGTGAGCTGGCCTCCACACAAAGCGAGCTCACCCTCAGCCAGAAGCACTGGCTGTCAGACCGCACCTACACCTGCCAGGTCACCTATCAAGGTCACACCTTTGAGGACAGCACCAAGAAGTGTGCAGATTCCAACCCGAGAGGGGTGAGCGCCTACCTAAGCCGGCCCAGCCCGTTCGACCTGTTCATCCGCAAGTCGCCCACGATCACCTGTCTGGTGGTGGACCTGGCACCCAGCAAGGGGACCGTGAACCTGACCTGGTCCCGGGCCAGTGGGAAGCCTGTGAACCACTCCACCAGAAAGGAGGAGAAGCAGCGCAATGGCACGTTAACCGTCACGTCCACCCTGCCGGTGGGCACCCGAGACTGGATCGAGGGGGAGACCTACCAGTGCAGGGTGACCCACCCCCACCTGCCCAGGGCCCTCATGCGGTCCACGACCAAGACCAGCGGCCCGCGTGCTGCCCCGGAAGTCTATGCGTTTGCGACGCCGGAGTGGCCGGGGAGCCGGGACAAGCGCACCCTCGCCTGCCTGATCCAGAACTTCATGCCTGAGGACATCTCGGTGCAGTGGCTGCACAACGAGGTGCAGCTCCCGGACGCCCGGCACAGCACGACGCAGCCCCGCAAGACCAAGGGCTCCGGCTTCTTCGTCTTCAGCCGCCTGGAGGTGACCAGGGCCGAATGGGAGCAGAAAGATGAGTTCATCTGCCGTGCAGTCCATGAGGCAGCGAGCCCCTCACAGACCGTCCAGCGAGCGGTGTCTGTAAATCCCGGTAAATGA'
        # IgEseq = 'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCTGCAAAAACATTCCCTCCAATGCCACCTCCGTGACTCTGGG'
        IgEseq = IgEseq[:SeqLen]
        # IgEseq = IgEseq[:30]
        Seq = ('IgE', IgEseq)
        ToClustalO.append(Seq)

        IgDseq = 'CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGGCATGCTTGATAACTGGGTACCACCCAACGTCCGTGACTGTCACCTGGTACATGGGGACACAGAGCCAGCCCCAGAGAACCTTCCCTGAGATACAAAGACGGGACAGCTACTACATGACAAGCAGCCAGCTCTCCACCCCCCTCCAGCAGTGGCGCCAAGGCGAGTACAAATGCGTGGTCCAGCACACCGCCAGCAAGAGTAAGAAGGAGATCTTCCGCTGGCCAGAGTCTCCAAAGGCACAGGCCTCCTCCGTGCCCACTGCACAACCCCAAGCAGAGGGCAGCCTCGCCAAGGCAACCACAGCCCCAGCCACCACCCGTAACACAGGAAGAGGAGGAGAAGAGAAGAAGAAGGAGAAGGAGAAAGAGGAACAAGAAGAGAGAGAGACAAAGACACCAGAGTGTCCGAGCCACACCCAGCCTCTTGGCGTCTACCTGCTAACCCCTGCAGTGCAGGACCTGTGGCTCCGGGACAAAGCCACCTTCACCTGCTTCGTGGTGGGCAGTGACCTGAAGGATGCTCACCTGACCTGGGAGGTGGCTGGGAAGGTCCCCACAGGGGGCGTGGAGGAAGGGCTGCTGGAGCGGCACAGCAACGGCTCCCAGAGCCAGCACAGCCGTCTGACCCTGCCCAGGTCCTTGTGGAACGCGGGGACCTCCGTCACCTGCACACTGAACCATCCCAGCCTCCCACCCCAGAGGTTGATGGCGCTGAGAGAACCCGCTGCGCAGGCACCCGTCAAGCTTTCTCTGAACCTGCTGGCCTCGTCTGACCCTCCCGAGGCGGCCTCGTGGCTCCTGTGTGAGGTGTCTGGCTTCTCGCCCCCCAACATCCTCCTGATGTGGCTGGAGGACCAGCGTGAGGTGAACACTTCTGGGTTTGCCCCCGCACGCCCCCCTCCACAGCCCAGGAGCACCACGTTCTGGGCCTGGAGTGTGCTGCGTGTCCCAGCCCCGCCCAGCCCTCAGCCAGCCACCTACACGTGTGTGGTCAGCCACGAGGACTCCCGGACTCTGCTCAACGCCAGCCGGAGCCTAGAAGTCAGCTACCTGGCCATGACCCCCCTGATCCCTCAGAGCAAGGATGAGAACAGCGATGACTACACGACCTTTGATGATGTGGGCAGCCTGTGGACCACCCTGTCCACGTTTGTGGCCCTCTTCATCCTCACCCTCCTCTACAGCGGCATTGTCACTTTCATCAAGGTGAAGTAG'
        # IgDseq = 'CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGG'
        IgDseq = IgDseq[:SeqLen]
        Seq = ('IgD', IgDseq)
        ToClustalO.append(Seq)

        ABVectorSeq = 'CGTCGACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCTGTGACGGTCTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGC'
        ABVectorSeq = ABVectorSeq[:SeqLen]
        Seq = ('Ab-vector', ABVectorSeq)
        ToClustalO.append(Seq)

        hit = "Unknown"

        IgSeqs = [IgMSeq, IgA1seq, IgA2seq, IgDseq, IgEseq, IgG1Seq, IgG2Seq, IgG3Seq, IgG4Seq, ABVectorSeq]

        best = difflib.get_close_matches(IsoSeq, IgSeqs, 1, 0.6)

        if len(best) > 0:
            for item in ToClustalO:
                # test = item[1]
                if item[1] == best[0]:
                    hit = item[0]

                    if hit == 'IgG1' or hit == 'IgG2':
                        if SeqLen < 35:
                            hit = 'IgG'
                    if hit == 'IgG3' or hit == 'IgG4':
                        if SeqLen < 20:
                            hit = 'IgG'
                    if hit == 'IgA1' or hit == 'IgA2':
                        if SeqLen < 39:
                            hit = 'IgA'
                    break



    return hit




    # SeqLen += 10
    #
    #
    # IsoOrderred = ClustalO(ToClustalO, wrapLength=SeqLen, ordered=True)
    #
    # outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-out-seqs.fa') # '/Applications/ClustalOmega/my-out-seqs.fa'
    # NoneFound = ''
    # NextUp = False
    # if os.path.isfile(outfilename):
    #     with open(outfilename, 'r') as currentfile:
    #         for line in currentfile:
    #             Readline = line.replace('\n', '').replace('\r', '').replace('-', '.')
    #             Readline = Readline.strip()
    #             if Readline[0] == '>':
    #                 # PrevSeq = SeqName
    #                 SeqName = Readline[1:]
    #                 if NextUp == True:
    #                     if SeqName == 'IgE':
    #                         currentfile.close()
    #                         ToClustalO.clear()
    #
    #                         IgG1Seq = IgG1Seq[:40]
    #                         Seq = ('IgG1', IgG1Seq)
    #                         ToClustalO.append(Seq)
    #
    #                         IgG2Seq = IgG2Seq[:40]
    #                         Seq = ('IgG2', IgG2Seq)
    #                         ToClustalO.append(Seq)
    #
    #                         IgG3Seq = IgG3Seq[:40]
    #                         Seq = ('IgG3', IgG3Seq)
    #                         ToClustalO.append(Seq)
    #
    #                         IgG4Seq = IgG4Seq[:40]
    #                         Seq = ('IgG4', IgG4Seq)
    #                         ToClustalO.append(Seq)
    #
    #                         IgEseq = IgEseq[:40]
    #                         Seq = ('IgE', IgEseq)
    #                         ToClustalO.append(Seq)
    #
    #                         Seq = ('TestSeq', IsoSeq)
    #                         ToClustalO.append(Seq)
    #
    #                         IsoOrderred = ClustalO(ToClustalO, wrapLength=SeqLen, ordered=True)
    #
    #                         outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-out-seqs.fa') #'/Applications/ClustalOmega/my-out-seqs.fa'
    #
    #                         NoneFound = ''
    #                         NextUp = False
    #                         if os.path.isfile(outfilename):
    #                             with open(outfilename, 'r') as currentfile:
    #                                 for line in currentfile:
    #                                     Readline = line.replace('\n', '').replace('\r', '').replace('-', '.')
    #                                     Readline = Readline.strip()
    #                                     if Readline[0] == '>':
    #                                         # PrevSeq = SeqName
    #                                         SeqName = Readline[1:]
    #                                         if NextUp == True:
    #                                             SeqName = SeqName[0:3]
    #                                             return SeqName
    #                                         if SeqName == 'TestSeq':
    #                                             NextUp = True
    #
    #
    #                     return SeqName
    #                 if SeqName == 'TestSeq':
    #                     NextUp = True
    #
    #             else:
    #                 NoneFound  = Readline
    #                 # if SeqName != 'TestSeq:':
    #                 #     StartAll = True
    #                 # else:
    #                 #
    #                 #     StartAll = False
    #     if NextUp == False:
    #         return NoneFound
    #
    # else:
    #     return "Unknown"

# template sequences for this function is downloaded from IMGT
# http://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?species=house%20mouse&latin=Mus%20musculus&group=IGHC
def CallIsotypeMouse(IsoSeq):
    import difflib

    IsoSeq = IsoSeq.replace('-', '')

    i = 0
    Ns = []
    SeqLen = len(IsoSeq)
    IsoSeq = IsoSeq.upper()

    for base in IsoSeq:
        i += 1
        if base == 'N':
            Count = (i,1)
            Ns.append(Count)
        else:
            Count = (i, 0)
            Ns.append(Count)
        if i > 9:
            Ncount = 0
            for j in range((i-10), i):
                Ncount += Ns[j][1]
            if Ncount >1:
                SeqLen = i
                IsoSeq = IsoSeq[:SeqLen]
                break

    CS = len(IsoSeq)
    if CS < 5:
        return 'Unknown'

    ToClustalO = []

    IgMSeq = 'AGAGTCAGTCCTTCCCAAATGTCTTCCCCCTCGTCTCCTGCGAGAGCCCCCTGTCTGATAAGAATCTGGTGGCCATGGGCTGCCTGGCCCGGGACTTCCTGCCCAGCACCATTTCCTTCACCTGGAACTACCAGAACAACACTGAAGTCATCCAGGGTATCAGAACCTTCCCAACACTGAGGACAGGGGGCAAGTACCTAGCCACCTCGCAGGTCTTGCTGTCTCCCAAGAGCATCCTTGAAGGTTCAGATGAATACCTTGTATGCAAAATCCACTACGGAGGCAAAAACAGAGATCTGCATGTGCCCATTCCAGCTGTCGCAGAGATGAATCCCAATGTAAATGTGTTCGTCCCACCACGGGATGGCTTCTCTGGCCCTGCACCACGCAAGTCTAAACTCATCTGCGAGGCCACGAACTTCACTCCAAAACCGATCACAGTATCCTGGCTAAAGGATGGGAAGCTCGTGGAATCTGGCTTCACCACAGATCCGGTGACCATCGAGAACAAAGGATCCACACCCCAAACCTACAAGGTCATAAGCACACTTACCATCTCTGAAATCGACTGGCTGAACCTGAATGTGTACACCTGCCGTGTGGATCACAGGGGTCTCACCTTCTTGAAGAACGTGTCCTCCACATGTGCTGCCAGTCCCTCCACAGACATCCTAACCTTCACCATCCCCCCCTCCTTTGCCGACATCTTCCTCAGCAAGTCCGCTAACCTGACCTGTCTGGTCTCAAACCTGGCAACCTATGAAACCCTGAATATCTCCTGGGCTTCTCAAAGTGGTGAACCACTGGAAACCAAAATTAAAATCATGGAAAGCCATCCCAATGGCACCTTCAGTGCTAAGGGTGTGGCTAGTGTTTGTGTGGAAGACTGGAATAACAGGAAGGAATTTGTGTGTACTGTGACTCACAGGGATCTGCCTTCGCCACAGAAGAAATTCATCTCAAAACCCAATG'
    IgMSeq = IgMSeq[:SeqLen]
    Seq = ('IgM', IgMSeq)
    ToClustalO.append(Seq)

    IgG1Seq = 'CCAAAACGACACCCCCATCTGTCTATCCACTGGCCCCTGGATCTGCTGCCCAAACTAACTCCATGGTGACCCTGGGATGCCTGGTCAAGGGCTATTTCCCTGAGCCAGTGACAGTGACCTGGAACTCTGGATCCCTGTCCAGCGGTGTGCACACCTTCCCAGCTGTCCTGGAGTCTGACCTCTACACTCTGAGCAGCTCAGTGACTGTCCCCTCCAGCCCTCGGCCCAGCGAGACCGTCACCTGCAACGTTGCCCACCCGGCCAGCAGCACCAAGGTGGACAAGAAAATTGTGCCCAGGGATTGTGGTTGTAAGCCTTGCATATGTACAGTCCCAGAAGTATCATCTGTCTTCATCTTCCCCCCAAAGCCCAAGGATGTGCTCACCATTACTCTGACTCCTAAGGTCACGTGTGTTGTGGTAGACATCAGCAAGGATGATCCCGAGGTCCAGTTCAGCTGGTTTGTAGATGATGTGGAGGTGCACACAGCTCAGACGCAACCCCGGGAGGAGCAGTTCAACAGCACTTTCCGCTCAGTCAGTGAACTTCCCATCATGCACCAGGACTGGCTCAATGGCAAGGAGTTCAAATGCAGGGTCAACAGTGCAGCTTTCCCTGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGCAGACCGAAGGCTCCACAGGTGTACACCATTCCACCTCCCAAGGAGCAGATGGCCAAGGATAAAGTCAGTCTGACCTGCATGATAACAGACTTCTTCCCTGAAGACATTACTGTGGAGTGGCAGTGGAATGGGCAGCCAGCGGAGAACTACAAGAACACTCAGCCCATCATGAACACGAATGGCTCTTACTTCGTCTACAGCAAGCTCAATGTGCAGAAGAGCAACTGGGAGGCAGGAAATACTTTCACCTGCTCTGTGTTACATGAGGGCCTGCACAACCACCATACTGAGAAGAGCCTCTCCCACTCTCCTGGTAAA'
    IgG1Seq = IgG1Seq[:SeqLen]
    Seq = ('IgG1', IgG1Seq)
    ToClustalO.append(Seq)

    IgG2ASeq = 'CCAAAACAACAGCCCCATCGGTCTATCCACTGGCCCCTGTGTGTGGAGATACAACTGGCTCCTCGGTGACTCTAGGATGCCTGGTCAAGGGTTATTTCCCTGAGCCAGTGACCTTGACCTGGAACTCTGGATCCCTGTCCAGTGGTGTGCACACCTTCCCAGCTGTCCTGCAGTCTGACCTCTACACCCTCAGCAGCTCAGTGACTGTAACCTCGAGCACCTGGCCCAGCCAGTCCATCACCTGCAATGTGGCCCACCCGGCAAGCAGCACCAAGGTGGACAAGAAAATTGAGCCCAGAGGGCCCACAATCAAGCCCTGTCCTCCATGCAAATGCCCAGCACCTAACCTCTTGGGTGGACCATCCGTCTTCATCTTCCCTCCAAAGATCAAGGATGTACTCATGATCTCCCTGAGCCCCATAGTCACATGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGATGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTCTCCGGGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTCAACAACAAAGACCTCCCAGCGCCCATCGAGAGAACCATCTCAAAACCCAAAGGGTCAGTAAGAGCTCCACAGGTATATGTCTTGCCTCCACCAGAAGAAGAGATGACTAAGAAACAGGTCACTCTGACCTGCATGGTCACAGACTTCATGCCTGAAGACATTTACGTGGAGTGGACCAACAACGGGAAAACAGAGCTAAACTACAAGAACACTGAACCAGTCCTGGACTCTGATGGTTCTTACTTCATGTACAGCAAGCTGAGAGTGGAAAAGAAGAACTGGGTGGAAAGAAATAGCTACTCCTGTTCAGTGGTCCACGAGGGTCTGCACAATCACCACACGACTAAGAGCTTCTCCCGGACTCCGGGTAAA'
    IgG2ASeq = IgG2ASeq[:SeqLen]
    Seq = ('IgG2A', IgG2ASeq)
    ToClustalO.append(Seq)

    IgG2BSeq = 'CCAAAACAACACCCCCATCAGTCTATCCACTGGCCCCTGGGTGTGGAGATACAACTGGTTCCTCCGTGACCTCTGGGTGCCTGGTCAAGGGGTACTTCCCTGAGCCAGTGACTGTGACTTGGAACTCTGGATCCCTGTCCAGCAGTGTGCACACCTTCCCAGCTCTCCTGCAGTCTGGACTCTACACTATGAGCAGCTCAGTGACTGTCCCCTCCAGCACCTGGCCAAGTCAGACCGTCACCTGCAGCGTTGCTCACCCAGCCAGCAGCACCACGGTGGACAAAAAACTTGAGCCCAGCGGGCCCATTTCAACAATCAACCCCTGTCCTCCATGCAAGGAGTGTCACAAATGCCCAGCTCCTAACCTCGAGGGTGGACCATCCGTCTTCATCTTCCCTCCAAATATCAAGGATGTACTCATGATCTCCCTGACACCCAAGGTCACGTGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGACGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTATCCGGGTGGTCAGCACCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTGAACAACAAAGACCTCCCATCACCCATCGAGAGAACCATCTCAAAAATTAAAGGGCTAGTCAGAGCTCCACAAGTATACACTTTGCCGCCACCAGCAGAGCAGTTGTCCAGGAAAGATGTCAGTCTCACTTGCCTGGTCGTGGGCTTCAACCCTGGAGACATCAGTGTGGAGTGGACCAGCAATGGGCATACAGAGGAGAACTACAAGGACACCGCACCAGTTCTTGACTCTGACGGTTCTTACTTCATATATAGCAAGCTCAATATGAAAACAAGCAAGTGGGAGAAAACAGATTCCTTCTCATGCAACGTGAGACACGAGGGTCTGAAAAATTACTACCTGAAGAAGACCATCTCCCGGTCTCCGGGTAAA'
    IgG2BSeq = IgG2BSeq[:SeqLen]
    Seq = ('IgG2B', IgG2BSeq)
    ToClustalO.append(Seq)

    IgG2CSeq = 'CCAAAACAACAGCCCCATCGGTCTATCCACTGGCCCCTGTGTGTGGAGGTACAACTGGCTCCTCGGTGACTCTAGGATGCCTGGTCAAGGGTTATTTCCCTGAGCCAGTGACCTTGACCTGGAACTCTGGATCCCTGTCCAGTGGTGTGCACACCTTCCCAGCTCTCCTGCAGTCTGGCCTCTACACCCTCAGCAGCTCAGTGACTGTAACCTCGAACACCTGGCCCAGCCAGACCATCACCTGCAATGTGGCCCACCCGGCAAGCAGCACCAAAGTGGACAAGAAAATTGAGCCCAGAGTGCCCATAACACAGAACCCCTGTCCTCCACTCAAAGAGTGTCCCCCATGCGCAGCTCCAGACCTCTTGGGTGGACCATCCGTCTTCATCTTCCCTCCAAAGATCAAGGATGTACTCATGATCTCCCTGAGCCCCATGGTCACATGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGACGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTCTCCGGGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTCAACAACAGAGCCCTCCCATCCCCCATCGAGAAAACCATCTCAAAACCCAGAGGGCCAGTAAGAGCTCCACAGGTATATGTCTTGCCTCCACCAGCAGAAGAGATGACTAAGAAAGAGTTCAGTCTGACCTGCATGATCACAGGCTTCTTACCTGCCGAAATTGCTGTGGACTGGACCAGCAATGGGCGTACAGAGCAAAACTACAAGAACACCGCAACAGTCCTGGACTCTGATGGTTCTTACTTCATGTACAGCAAGCTCAGAGTACAAAAGAGCACTTGGGAAAGAGGAAGTCTTTTCGCCTGCTCAGTGGTCCACGAGGTGCTGCACAATCACCTTACGACTAAGACCATCTCCCGGTCTCTGGGTAAA'
    IgG2CSeq = IgG2CSeq[:SeqLen]
    Seq = ('IgG2C', IgG2CSeq)
    ToClustalO.append(Seq)

    IgG3seq = 'CTACAACAACAGCCCCATCTGTCTATCCCTTGGTCCCTGGCTGCAGTGACACATCTGGATCCTCGGTGACACTGGGATGCCTTGTCAAAGGCTACTTCCCTGAGCCGGTAACTGTAAAATGGAACTATGGAGCCCTGTCCAGCGGTGTGCGCACAGTCTCATCTGTCCTGCAGTCTGGGTTCTATTCCCTCAGCAGCTTGGTGACTGTACCCTCCAGCACCTGGCCCAGCCAGACTGTCATCTGCAACGTAGCCCACCCAGCCAGCAAGACTGAGTTGATCAAGAGAATCGAGCCTAGAATACCCAAGCCCAGTACCCCCCCAGGTTCTTCATGCCCACCTGGTAACATCTTGGGTGGACCATCCGTCTTCATCTTCCCCCCAAAGCCCAAGGATGCACTCATGATCTCCCTAACCCCCAAGGTTACGTGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGATGTCCATGTCAGCTGGTTTGTGGACAACAAAGAAGTACACACAGCCTGGACACAGCCCCGTGAAGCTCAGTACAACAGTACCTTCCGAGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGGGGCAAGGAGTTCAAATGCAAGGTCAACAACAAAGCCCTCCCAGCCCCCATCGAGAGAACCATCTCAAAACCCAAAGGAAGAGCCCAGACACCTCAAGTATACACCATACCCCCACCTCGTGAACAAATGTCCAAGAAGAAGGTTAGTCTGACCTGCCTGGTCACCAACTTCTTCTCTGAAGCCATCAGTGTGGAGTGGGAAAGGAACGGAGAACTGGAGCAGGATTACAAGAACACTCCACCCATCCTGGACTCAGATGGGACCTACTTCCTCTACAGCAAGCTCACTGTGGATACAGACAGTTGGTTGCAAGGAGAAATTTTTACCTGCTCCGTGGTGCATGAGGCTCTCCATAACCACCACACACAGAAGAACCTGTCTCGCTCCCCTGGTAAA'
    IgG3seq = IgG3seq[:SeqLen]
    Seq = ('IgG3', IgG3seq)
    ToClustalO.append(Seq)

    IgAseq = 'AGTCTGCGAGAAATCCCACCATCTACCCACTGACACTCCCACCAGTCCTGTGCAGTGATCCCGTGATAATCGGCTGCCTGATTCACGATTACTTCCCTTTCGGCACGATGAATGTGACCTGGGGAAAGAGTGGGAAGGATATAACCACCGTGAACTTTCCACCTGCCCTCGCCTCTGGGGGACGGTACACCATGAGCAGCCAGTTAACCCTGCCAGCTGTCGAGTGCCCAGAAGGAGAGTCCGTGAAATGTTCCGTGCAACATGACTCTAACCCCGTCCAAGAATTGGATGTGAATTGCTCTGGTCCTACTCCTCCTCCTCCTATTACTATTCCTTCCTGCCAGCCCAGCCTGTCACTGCAGCGGCCAGCTCTTGAGGACCTGCTCCTGGGTTCAGATGCCAGCATCACATGTACTCTGAATGGCCTGAGAAATCCTGAGGGAGCTGCTTTCACCTGGGAGCCCTCCACTGGGAAGGATGCAGTGCAGAAGAAAGCTGCGCAGAATTCCTGCGGCTGCTACAGTGTGTCCAGCGTCCTGCCTGGCTGTGCTGAGCGCTGGAACAGTGGCGCATCATTCAAGTGCACAGTTACCCATCCTGAGTCTGGCACCTTAACTGGCACAATTGCCAAAGTCACAGTGAACACCTTCCCACCCCAGGTCCACCTGCTACCGCCGCCGTCGGAGGAGCTGGCCCTGAATGAGCTCTTGTCCCTGACATGCCTGGTGCGAGCTTTCAACCCTAAAGAAGTGCTGGTGCGATGGCTGCATGGAAATGAGGAGCTGTCCCCAGAAAGCTACCTAGTGTTTGAGCCCCTAAAGGAGCCAGGCGAGGGAGCCACCACCTACCTGGTGACAAGCGTGTTGCGTGTATCAGCTGAAACCTGGAAACAGGGTGACCAGTACTCCTGCATGGTGGGCCACGAGGCCTTGCCCATGAACTTCACCCAGAAGACCATCGACCGTCTGTCGGGTAAACCCACCAATGTCAGCGTGTCTGTGATCATGTCAGAGGGAGATGGCATCTGCTAC'
    IgAseq = IgAseq[:SeqLen]
    Seq = ('IgA', IgAseq)
    ToClustalO.append(Seq)

    IgEseq = 'CCTCTATCAGGAACCCTCAGCTCTACCCCTTAAAGCCCTGTAAAGGCACTGCTTCCATGACCCTAGGCTGCCTAGTAAAGGACTACTTCCCTAATCCTGTGACTGTGACCTGGTATTCAGACTCCCTGAACATGAGCACTGTGAACTTCCCTGCCCTCGGTTCTGAACTCAAGGTCACCACCAGCCAAGTGACCAGCTGGGGCAAGTCAGCCAAGAACTTCACATGCCACGTGACACATCCTCCATCATTCAACGAAAGTAGGACTATCCTAGTTCGACCTGTCAACATCACTGAGCCCACCTTGGAGCTACTCCATTCATCCTGCGACCCCAATGCATTCCACTCCACCATCCAGCTGTACTGCTTCATTTATGGCCACATCCTAAATGATGTCTCTGTCAGCTGGCTAATGGACGATCGGGAGATAACTGATACACTTGCACAAACTGTTCTAATCAAGGAGGAAGGCAAACTAGCCTCTACCTGCAGTAAACTCAACATCACTGAGCAGCAATGGATGTCTGAAAGCACCTTCACCTGCAAGGTCACCTCCCAAGGCGTAGACTATTTGGCCCACACTCGGAGATGCCCAGATCATGAGCCACGGGGTGTGATTACCTACCTGATCCCACCCAGCCCCCTGGACCTGTATCAAAACGGTGCTCCCAAGCTTACCTGTCTGGTGGTGGACCTGGAAAGCGAGAAGAATGTCAATGTGACGTGGAACCAAGAGAAGAAGACTTCAGTCTCAGCATCCCAGTGGTACACTAAGCACCACAATAACGCCACAACTAGTATCACCTCCATCCTGCCTGTAGTTGCCAAGGACTGGATTGAAGGCTACGGCTATCAGTGCATAGTGGACCACCCTGATTTTCCCAAGCCCATTGTGCGTTCCATCACCAAGACCCCAGGCCAGCGCTCAGCCCCCGAGGTATATGTGTTCCCACCACCAGAGGAGGAGAGCGAGGACAAACGCACACTCACCTGTTTGATCCAGAACTTCTTCCCTGAGGATATCTCTGTGCAGTGGCTGGGGGATGGCAAACTGATCTCAAACAGCCAGCACAGTACCACAACACCCCTGAAATCCAATGGCTCCAATCAAGGCTTCTTCATCTTCAGTCGCCTAGAGGTCGCCAAGACACTCTGGACACAGAGAAAACAGTTCACCTGCCAAGTGATCCATGAGGCACTTCAGAAACCCAGGAAACTGGAGAAAACAATATCCACAAGCCTTGGTAACACCTCCCTCCGTCCCTCC'
    IgEseq = IgEseq[:SeqLen]
    Seq = ('IgE', IgEseq)
    ToClustalO.append(Seq)

    IgDseq = 'GTGATAAAAAGGAACCTGACATGTTCCTCCTCTCAGAGTGCAAAGCCCCAGAGGAAAATGAAAAGATAAACCTGGGCTGTTTAGTAATTGGAAGTCAGCCACTGAAAATCAGCTGGGAGCCAAAGAAGTCAAGTATAGTTGAACATGTCTTCCCCTCTGAAATGAGAAATGGCAATTATACAATGGTCCTCCAGGTCACTGTGCTGGCCTCAGAACTGAACCTCAACCACACTTGCACCATAAATAAACCCAAAAGGAAAGAAAAACCTTTCAAGTTTCCTGAGTCATGGGATTCCCAGTCCTCTAAGAGAGTCACTCCAACTCTCCAAGCAAAGAATCACTCCACAGAAGCCACCAAAGCTATTACCACCAAAAAGGACATAGAAGGGGCCATGGCACCCAGCAACCTCACTGTGAACATCCTGACCACATCCACCCATCCTGAGATGTCATCTTGGCTCCTGTGTGAAGTATCTGGCTTCTTCCCGGAAAATATCCACCTCATGTGGCTGGGTGTCCACAGTAAAATGAAGTCTACAAACTTTGTCACTGCAAACCCCACCGCCCAGCCTGGGGGCACATTCCAGACCTGGAGTGTCCTGAGACTACCAGTCGCTCTGAGCTCATCACTTGACACTTACACATGTGTGGTGGAACATGAGGCCTCAAAGACAAAGCTTAATGCCAGCAAGAGCCTAGCAATTAGTG'
    IgDseq = IgDseq[:SeqLen]
    Seq = ('IgD', IgDseq)
    ToClustalO.append(Seq)

    ABVectorSeq = 'CGTCGACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCTGTGACGGTCTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGC'
    ABVectorSeq = ABVectorSeq[:SeqLen]
    Seq = ('Ab-vector', ABVectorSeq)
    ToClustalO.append(Seq)


    hit = "Unknown-default"

    IgSeqs = [IgMSeq, IgG1Seq, IgG2ASeq, IgG2BSeq, IgG2CSeq, IgG3seq, IgAseq, IgDseq,IgEseq,ABVectorSeq]

    best = difflib.get_close_matches(IsoSeq, IgSeqs, 1, 0.6)

    if len(best)>0:
        for item in ToClustalO:
            if item[1] == best[0]:
                hit = item[0]
                break
    else:
        if SeqLen >39:
            SeqLen = 40
        elif SeqLen >19:
            SeqLen = 20
        else:
            return 'Unknown'

        IsoSeq = IsoSeq[:SeqLen]
        CS = len(IsoSeq)
        if CS < 5:
            return 'Unknown'

        ToClustalO = []

        IgMSeq = 'AGAGTCAGTCCTTCCCAAATGTCTTCCCCCTCGTCTCCTGCGAGAGCCCCCTGTCTGATAAGAATCTGGTGGCCATGGGCTGCCTGGCCCGGGACTTCCTGCCCAGCACCATTTCCTTCACCTGGAACTACCAGAACAACACTGAAGTCATCCAGGGTATCAGAACCTTCCCAACACTGAGGACAGGGGGCAAGTACCTAGCCACCTCGCAGGTCTTGCTGTCTCCCAAGAGCATCCTTGAAGGTTCAGATGAATACCTTGTATGCAAAATCCACTACGGAGGCAAAAACAGAGATCTGCATGTGCCCATTCCAGCTGTCGCAGAGATGAATCCCAATGTAAATGTGTTCGTCCCACCACGGGATGGCTTCTCTGGCCCTGCACCACGCAAGTCTAAACTCATCTGCGAGGCCACGAACTTCACTCCAAAACCGATCACAGTATCCTGGCTAAAGGATGGGAAGCTCGTGGAATCTGGCTTCACCACAGATCCGGTGACCATCGAGAACAAAGGATCCACACCCCAAACCTACAAGGTCATAAGCACACTTACCATCTCTGAAATCGACTGGCTGAACCTGAATGTGTACACCTGCCGTGTGGATCACAGGGGTCTCACCTTCTTGAAGAACGTGTCCTCCACATGTGCTGCCAGTCCCTCCACAGACATCCTAACCTTCACCATCCCCCCCTCCTTTGCCGACATCTTCCTCAGCAAGTCCGCTAACCTGACCTGTCTGGTCTCAAACCTGGCAACCTATGAAACCCTGAATATCTCCTGGGCTTCTCAAAGTGGTGAACCACTGGAAACCAAAATTAAAATCATGGAAAGCCATCCCAATGGCACCTTCAGTGCTAAGGGTGTGGCTAGTGTTTGTGTGGAAGACTGGAATAACAGGAAGGAATTTGTGTGTACTGTGACTCACAGGGATCTGCCTTCGCCACAGAAGAAATTCATCTCAAAACCCAATG'
        IgMSeq = IgMSeq[:SeqLen]
        Seq = ('IgM', IgMSeq)
        ToClustalO.append(Seq)

        IgG1Seq = 'CCAAAACGACACCCCCATCTGTCTATCCACTGGCCCCTGGATCTGCTGCCCAAACTAACTCCATGGTGACCCTGGGATGCCTGGTCAAGGGCTATTTCCCTGAGCCAGTGACAGTGACCTGGAACTCTGGATCCCTGTCCAGCGGTGTGCACACCTTCCCAGCTGTCCTGGAGTCTGACCTCTACACTCTGAGCAGCTCAGTGACTGTCCCCTCCAGCCCTCGGCCCAGCGAGACCGTCACCTGCAACGTTGCCCACCCGGCCAGCAGCACCAAGGTGGACAAGAAAATTGTGCCCAGGGATTGTGGTTGTAAGCCTTGCATATGTACAGTCCCAGAAGTATCATCTGTCTTCATCTTCCCCCCAAAGCCCAAGGATGTGCTCACCATTACTCTGACTCCTAAGGTCACGTGTGTTGTGGTAGACATCAGCAAGGATGATCCCGAGGTCCAGTTCAGCTGGTTTGTAGATGATGTGGAGGTGCACACAGCTCAGACGCAACCCCGGGAGGAGCAGTTCAACAGCACTTTCCGCTCAGTCAGTGAACTTCCCATCATGCACCAGGACTGGCTCAATGGCAAGGAGTTCAAATGCAGGGTCAACAGTGCAGCTTTCCCTGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGCAGACCGAAGGCTCCACAGGTGTACACCATTCCACCTCCCAAGGAGCAGATGGCCAAGGATAAAGTCAGTCTGACCTGCATGATAACAGACTTCTTCCCTGAAGACATTACTGTGGAGTGGCAGTGGAATGGGCAGCCAGCGGAGAACTACAAGAACACTCAGCCCATCATGAACACGAATGGCTCTTACTTCGTCTACAGCAAGCTCAATGTGCAGAAGAGCAACTGGGAGGCAGGAAATACTTTCACCTGCTCTGTGTTACATGAGGGCCTGCACAACCACCATACTGAGAAGAGCCTCTCCCACTCTCCTGGTAAA'
        IgG1Seq = IgG1Seq[:SeqLen]
        Seq = ('IgG1', IgG1Seq)
        ToClustalO.append(Seq)

        IgG2ASeq = 'CCAAAACAACAGCCCCATCGGTCTATCCACTGGCCCCTGTGTGTGGAGATACAACTGGCTCCTCGGTGACTCTAGGATGCCTGGTCAAGGGTTATTTCCCTGAGCCAGTGACCTTGACCTGGAACTCTGGATCCCTGTCCAGTGGTGTGCACACCTTCCCAGCTGTCCTGCAGTCTGACCTCTACACCCTCAGCAGCTCAGTGACTGTAACCTCGAGCACCTGGCCCAGCCAGTCCATCACCTGCAATGTGGCCCACCCGGCAAGCAGCACCAAGGTGGACAAGAAAATTGAGCCCAGAGGGCCCACAATCAAGCCCTGTCCTCCATGCAAATGCCCAGCACCTAACCTCTTGGGTGGACCATCCGTCTTCATCTTCCCTCCAAAGATCAAGGATGTACTCATGATCTCCCTGAGCCCCATAGTCACATGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGATGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTCTCCGGGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTCAACAACAAAGACCTCCCAGCGCCCATCGAGAGAACCATCTCAAAACCCAAAGGGTCAGTAAGAGCTCCACAGGTATATGTCTTGCCTCCACCAGAAGAAGAGATGACTAAGAAACAGGTCACTCTGACCTGCATGGTCACAGACTTCATGCCTGAAGACATTTACGTGGAGTGGACCAACAACGGGAAAACAGAGCTAAACTACAAGAACACTGAACCAGTCCTGGACTCTGATGGTTCTTACTTCATGTACAGCAAGCTGAGAGTGGAAAAGAAGAACTGGGTGGAAAGAAATAGCTACTCCTGTTCAGTGGTCCACGAGGGTCTGCACAATCACCACACGACTAAGAGCTTCTCCCGGACTCCGGGTAAA'
        IgG2ASeq = IgG2ASeq[:SeqLen]
        Seq = ('IgG2A', IgG2ASeq)
        ToClustalO.append(Seq)

        IgG2BSeq = 'CCAAAACAACACCCCCATCAGTCTATCCACTGGCCCCTGGGTGTGGAGATACAACTGGTTCCTCCGTGACCTCTGGGTGCCTGGTCAAGGGGTACTTCCCTGAGCCAGTGACTGTGACTTGGAACTCTGGATCCCTGTCCAGCAGTGTGCACACCTTCCCAGCTCTCCTGCAGTCTGGACTCTACACTATGAGCAGCTCAGTGACTGTCCCCTCCAGCACCTGGCCAAGTCAGACCGTCACCTGCAGCGTTGCTCACCCAGCCAGCAGCACCACGGTGGACAAAAAACTTGAGCCCAGCGGGCCCATTTCAACAATCAACCCCTGTCCTCCATGCAAGGAGTGTCACAAATGCCCAGCTCCTAACCTCGAGGGTGGACCATCCGTCTTCATCTTCCCTCCAAATATCAAGGATGTACTCATGATCTCCCTGACACCCAAGGTCACGTGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGACGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTATCCGGGTGGTCAGCACCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTGAACAACAAAGACCTCCCATCACCCATCGAGAGAACCATCTCAAAAATTAAAGGGCTAGTCAGAGCTCCACAAGTATACACTTTGCCGCCACCAGCAGAGCAGTTGTCCAGGAAAGATGTCAGTCTCACTTGCCTGGTCGTGGGCTTCAACCCTGGAGACATCAGTGTGGAGTGGACCAGCAATGGGCATACAGAGGAGAACTACAAGGACACCGCACCAGTTCTTGACTCTGACGGTTCTTACTTCATATATAGCAAGCTCAATATGAAAACAAGCAAGTGGGAGAAAACAGATTCCTTCTCATGCAACGTGAGACACGAGGGTCTGAAAAATTACTACCTGAAGAAGACCATCTCCCGGTCTCCGGGTAAA'
        IgG2BSeq = IgG2BSeq[:SeqLen]
        Seq = ('IgG2B', IgG2BSeq)
        ToClustalO.append(Seq)

        IgG2CSeq = 'CCAAAACAACAGCCCCATCGGTCTATCCACTGGCCCCTGTGTGTGGAGGTACAACTGGCTCCTCGGTGACTCTAGGATGCCTGGTCAAGGGTTATTTCCCTGAGCCAGTGACCTTGACCTGGAACTCTGGATCCCTGTCCAGTGGTGTGCACACCTTCCCAGCTCTCCTGCAGTCTGGCCTCTACACCCTCAGCAGCTCAGTGACTGTAACCTCGAACACCTGGCCCAGCCAGACCATCACCTGCAATGTGGCCCACCCGGCAAGCAGCACCAAAGTGGACAAGAAAATTGAGCCCAGAGTGCCCATAACACAGAACCCCTGTCCTCCACTCAAAGAGTGTCCCCCATGCGCAGCTCCAGACCTCTTGGGTGGACCATCCGTCTTCATCTTCCCTCCAAAGATCAAGGATGTACTCATGATCTCCCTGAGCCCCATGGTCACATGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGACGTCCAGATCAGCTGGTTTGTGAACAACGTGGAAGTACACACAGCTCAGACACAAACCCATAGAGAGGATTACAACAGTACTCTCCGGGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGTGGCAAGGAGTTCAAATGCAAGGTCAACAACAGAGCCCTCCCATCCCCCATCGAGAAAACCATCTCAAAACCCAGAGGGCCAGTAAGAGCTCCACAGGTATATGTCTTGCCTCCACCAGCAGAAGAGATGACTAAGAAAGAGTTCAGTCTGACCTGCATGATCACAGGCTTCTTACCTGCCGAAATTGCTGTGGACTGGACCAGCAATGGGCGTACAGAGCAAAACTACAAGAACACCGCAACAGTCCTGGACTCTGATGGTTCTTACTTCATGTACAGCAAGCTCAGAGTACAAAAGAGCACTTGGGAAAGAGGAAGTCTTTTCGCCTGCTCAGTGGTCCACGAGGTGCTGCACAATCACCTTACGACTAAGACCATCTCCCGGTCTCTGGGTAAA'
        IgG2CSeq = IgG2CSeq[:SeqLen]
        Seq = ('IgG2C', IgG2CSeq)
        ToClustalO.append(Seq)

        IgG3seq = 'CTACAACAACAGCCCCATCTGTCTATCCCTTGGTCCCTGGCTGCAGTGACACATCTGGATCCTCGGTGACACTGGGATGCCTTGTCAAAGGCTACTTCCCTGAGCCGGTAACTGTAAAATGGAACTATGGAGCCCTGTCCAGCGGTGTGCGCACAGTCTCATCTGTCCTGCAGTCTGGGTTCTATTCCCTCAGCAGCTTGGTGACTGTACCCTCCAGCACCTGGCCCAGCCAGACTGTCATCTGCAACGTAGCCCACCCAGCCAGCAAGACTGAGTTGATCAAGAGAATCGAGCCTAGAATACCCAAGCCCAGTACCCCCCCAGGTTCTTCATGCCCACCTGGTAACATCTTGGGTGGACCATCCGTCTTCATCTTCCCCCCAAAGCCCAAGGATGCACTCATGATCTCCCTAACCCCCAAGGTTACGTGTGTGGTGGTGGATGTGAGCGAGGATGACCCAGATGTCCATGTCAGCTGGTTTGTGGACAACAAAGAAGTACACACAGCCTGGACACAGCCCCGTGAAGCTCAGTACAACAGTACCTTCCGAGTGGTCAGTGCCCTCCCCATCCAGCACCAGGACTGGATGAGGGGCAAGGAGTTCAAATGCAAGGTCAACAACAAAGCCCTCCCAGCCCCCATCGAGAGAACCATCTCAAAACCCAAAGGAAGAGCCCAGACACCTCAAGTATACACCATACCCCCACCTCGTGAACAAATGTCCAAGAAGAAGGTTAGTCTGACCTGCCTGGTCACCAACTTCTTCTCTGAAGCCATCAGTGTGGAGTGGGAAAGGAACGGAGAACTGGAGCAGGATTACAAGAACACTCCACCCATCCTGGACTCAGATGGGACCTACTTCCTCTACAGCAAGCTCACTGTGGATACAGACAGTTGGTTGCAAGGAGAAATTTTTACCTGCTCCGTGGTGCATGAGGCTCTCCATAACCACCACACACAGAAGAACCTGTCTCGCTCCCCTGGTAAA'
        IgG3seq = IgG3seq[:SeqLen]
        Seq = ('IgG3', IgG3seq)
        ToClustalO.append(Seq)

        IgAseq = 'AGTCTGCGAGAAATCCCACCATCTACCCACTGACACTCCCACCAGTCCTGTGCAGTGATCCCGTGATAATCGGCTGCCTGATTCACGATTACTTCCCTTTCGGCACGATGAATGTGACCTGGGGAAAGAGTGGGAAGGATATAACCACCGTGAACTTTCCACCTGCCCTCGCCTCTGGGGGACGGTACACCATGAGCAGCCAGTTAACCCTGCCAGCTGTCGAGTGCCCAGAAGGAGAGTCCGTGAAATGTTCCGTGCAACATGACTCTAACCCCGTCCAAGAATTGGATGTGAATTGCTCTGGTCCTACTCCTCCTCCTCCTATTACTATTCCTTCCTGCCAGCCCAGCCTGTCACTGCAGCGGCCAGCTCTTGAGGACCTGCTCCTGGGTTCAGATGCCAGCATCACATGTACTCTGAATGGCCTGAGAAATCCTGAGGGAGCTGCTTTCACCTGGGAGCCCTCCACTGGGAAGGATGCAGTGCAGAAGAAAGCTGCGCAGAATTCCTGCGGCTGCTACAGTGTGTCCAGCGTCCTGCCTGGCTGTGCTGAGCGCTGGAACAGTGGCGCATCATTCAAGTGCACAGTTACCCATCCTGAGTCTGGCACCTTAACTGGCACAATTGCCAAAGTCACAGTGAACACCTTCCCACCCCAGGTCCACCTGCTACCGCCGCCGTCGGAGGAGCTGGCCCTGAATGAGCTCTTGTCCCTGACATGCCTGGTGCGAGCTTTCAACCCTAAAGAAGTGCTGGTGCGATGGCTGCATGGAAATGAGGAGCTGTCCCCAGAAAGCTACCTAGTGTTTGAGCCCCTAAAGGAGCCAGGCGAGGGAGCCACCACCTACCTGGTGACAAGCGTGTTGCGTGTATCAGCTGAAACCTGGAAACAGGGTGACCAGTACTCCTGCATGGTGGGCCACGAGGCCTTGCCCATGAACTTCACCCAGAAGACCATCGACCGTCTGTCGGGTAAACCCACCAATGTCAGCGTGTCTGTGATCATGTCAGAGGGAGATGGCATCTGCTAC'
        IgAseq = IgAseq[:SeqLen]
        Seq = ('IgA', IgAseq)
        ToClustalO.append(Seq)

        IgEseq = 'CCTCTATCAGGAACCCTCAGCTCTACCCCTTAAAGCCCTGTAAAGGCACTGCTTCCATGACCCTAGGCTGCCTAGTAAAGGACTACTTCCCTAATCCTGTGACTGTGACCTGGTATTCAGACTCCCTGAACATGAGCACTGTGAACTTCCCTGCCCTCGGTTCTGAACTCAAGGTCACCACCAGCCAAGTGACCAGCTGGGGCAAGTCAGCCAAGAACTTCACATGCCACGTGACACATCCTCCATCATTCAACGAAAGTAGGACTATCCTAGTTCGACCTGTCAACATCACTGAGCCCACCTTGGAGCTACTCCATTCATCCTGCGACCCCAATGCATTCCACTCCACCATCCAGCTGTACTGCTTCATTTATGGCCACATCCTAAATGATGTCTCTGTCAGCTGGCTAATGGACGATCGGGAGATAACTGATACACTTGCACAAACTGTTCTAATCAAGGAGGAAGGCAAACTAGCCTCTACCTGCAGTAAACTCAACATCACTGAGCAGCAATGGATGTCTGAAAGCACCTTCACCTGCAAGGTCACCTCCCAAGGCGTAGACTATTTGGCCCACACTCGGAGATGCCCAGATCATGAGCCACGGGGTGTGATTACCTACCTGATCCCACCCAGCCCCCTGGACCTGTATCAAAACGGTGCTCCCAAGCTTACCTGTCTGGTGGTGGACCTGGAAAGCGAGAAGAATGTCAATGTGACGTGGAACCAAGAGAAGAAGACTTCAGTCTCAGCATCCCAGTGGTACACTAAGCACCACAATAACGCCACAACTAGTATCACCTCCATCCTGCCTGTAGTTGCCAAGGACTGGATTGAAGGCTACGGCTATCAGTGCATAGTGGACCACCCTGATTTTCCCAAGCCCATTGTGCGTTCCATCACCAAGACCCCAGGCCAGCGCTCAGCCCCCGAGGTATATGTGTTCCCACCACCAGAGGAGGAGAGCGAGGACAAACGCACACTCACCTGTTTGATCCAGAACTTCTTCCCTGAGGATATCTCTGTGCAGTGGCTGGGGGATGGCAAACTGATCTCAAACAGCCAGCACAGTACCACAACACCCCTGAAATCCAATGGCTCCAATCAAGGCTTCTTCATCTTCAGTCGCCTAGAGGTCGCCAAGACACTCTGGACACAGAGAAAACAGTTCACCTGCCAAGTGATCCATGAGGCACTTCAGAAACCCAGGAAACTGGAGAAAACAATATCCACAAGCCTTGGTAACACCTCCCTCCGTCCCTCC'
        IgEseq = IgEseq[:SeqLen]
        Seq = ('IgE', IgEseq)
        ToClustalO.append(Seq)

        IgDseq = 'GTGATAAAAAGGAACCTGACATGTTCCTCCTCTCAGAGTGCAAAGCCCCAGAGGAAAATGAAAAGATAAACCTGGGCTGTTTAGTAATTGGAAGTCAGCCACTGAAAATCAGCTGGGAGCCAAAGAAGTCAAGTATAGTTGAACATGTCTTCCCCTCTGAAATGAGAAATGGCAATTATACAATGGTCCTCCAGGTCACTGTGCTGGCCTCAGAACTGAACCTCAACCACACTTGCACCATAAATAAACCCAAAAGGAAAGAAAAACCTTTCAAGTTTCCTGAGTCATGGGATTCCCAGTCCTCTAAGAGAGTCACTCCAACTCTCCAAGCAAAGAATCACTCCACAGAAGCCACCAAAGCTATTACCACCAAAAAGGACATAGAAGGGGCCATGGCACCCAGCAACCTCACTGTGAACATCCTGACCACATCCACCCATCCTGAGATGTCATCTTGGCTCCTGTGTGAAGTATCTGGCTTCTTCCCGGAAAATATCCACCTCATGTGGCTGGGTGTCCACAGTAAAATGAAGTCTACAAACTTTGTCACTGCAAACCCCACCGCCCAGCCTGGGGGCACATTCCAGACCTGGAGTGTCCTGAGACTACCAGTCGCTCTGAGCTCATCACTTGACACTTACACATGTGTGGTGGAACATGAGGCCTCAAAGACAAAGCTTAATGCCAGCAAGAGCCTAGCAATTAGTG'
        IgDseq = IgDseq[:SeqLen]
        Seq = ('IgD', IgDseq)
        ToClustalO.append(Seq)

        ABVectorSeq = 'CGTCGACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCTGTGACGGTCTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGC'
        ABVectorSeq = ABVectorSeq[:SeqLen]
        Seq = ('Ab-vector', ABVectorSeq)
        ToClustalO.append(Seq)

        hit = "Unknown"

        IgSeqs = [IgMSeq, IgG1Seq, IgG2ASeq, IgG2BSeq, IgG2CSeq, IgG3seq, IgAseq, IgDseq,IgEseq,ABVectorSeq]

        best = difflib.get_close_matches(IsoSeq, IgSeqs, 1, 0.6)

        if len(best) > 0:
            for item in ToClustalO:
                # test = item[1]
                if item[1] == best[0]:
                    hit = item[0]
                    break

    return hit



def IDPlacer(TestSeq, GermSeq):
    i = 0
    SeqLen = len(GermSeq)
    Seqlen2 = len(TestSeq)
    try:
        for i in range(i,SeqLen-1):
            if i < SeqLen-1:
                if TestSeq[i] == '-':  #start deletion
                    pos = i
                    IDSeq = ''
                    j = i
                    IDnuc = TestSeq[i]
                    while IDnuc == '-':
                        IDSeq += GermSeq[j]
                        j +=1
                        IDnuc = TestSeq[j]
                    if len(IDSeq) == 1:
                        if TestSeq[i+1] != GermSeq[i+1]:  #the next base downstream doesn't match
                            if TestSeq[i+1] == GermSeq[i]:
                                NewSeq = TestSeq[:i]
                                j =1

                                while TestSeq[i+j] == GermSeq[i+(j-1)]:
                                    NewSeq += TestSeq[i+j] #  #move one over until doesn't match
                                    j +=1
                                NewSeq += '-' + TestSeq[len(NewSeq)+1:]
                                TestSeq = NewSeq

                        elif TestSeq[i-1] != GermSeq[i-1]:  #the next base upstream doesn't match
                            test1 = TestSeq[i-1]
                            test2 = GermSeq[i]
                            if TestSeq[i-1] == GermSeq[i]:
                                NewSeq = TestSeq[i+1:]
                                j =1
                                moveBack = 0

                                while TestSeq[i-j] == GermSeq[i-(j-1)]:
                                    moveBack += 1
                                    NewSeq = TestSeq[i-j] + NewSeq #  #move one over until doesn't match
                                    j +=1

                                NewSeq = TestSeq[:i-moveBack] + '-' + NewSeq
                                TestSeq = NewSeq
                    else:
                        print('s')#We have deleted seq as IDSeq

    except:
        return TestSeq


    return TestSeq#, GermSeq



def ClustalO_new(SeqDict, wrapLength, ordered, working_prefix, bin):
    # input is a list of lists containing filename and sequence: ((fielname1, seq1),(fielname2, seq2))
    # import time
    # clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v
    # workingfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'ClustalOmega', 'my-in-seqs.fa')
    # workingfilename = '/Applications/ClustalOmega/my-in-seqs.fa'
    # workingdir, NameBase = os.path.split(DBname)
    #
    #
    #
    # NameBase = NameBase[:(len(NameBase) - 4)]

    import uuid

    NameBase  =str(uuid.uuid4())
    # NameBase = NameBase[:12]
    NameBase = NameBase.replace('-', '')

    NameBase  = NameBase.replace(' ', '')


    MyInFiles = NameBase + 'In.fa'
    MyOutFiles = NameBase + 'Out.fa'

    workingfilename = os.path.join(working_prefix, MyInFiles)
    savefilename = os.path.join(working_prefix,  MyOutFiles)

    workingdir, filename = os.path.split(workingfilename)

    os.chdir(workingdir)
    i = 1
    FASTAfile = ''
    for seq in SeqDict:
        for item in seq:
            if i % 2 != 0:
                FASTAfile += '>' + item + '\n'
                if i == 1:
                    SeqName = FASTAfile
            else:
                if item == '':
                    item = 'N'
                if SeqName == "Germline":
                    if len(item) < 50:
                        return False

                FASTAfile += item + '\n'
            i += 1

    with open(workingfilename, 'w') as currentFile:  # using with for this automatically closes the file even if you crash
        currentFile.write(FASTAfile)
    # todo with the ./ in front of clustalO all I need is the 'stand alone MAC binary' from http://www.clustal.org/omega/ no other installation
    # Mac users should rename the downloaded file to clustalo and
    #  place in the location of their choice. This file may need to
    # be made executable e.g.: chmod u+x clustalo
    # ClustalOCommandLine = "clustalo -i my-in-seqs.fa -o my-out-seqs.fa -v --full --force"  # --full uses better alignment but slower
    # with open('/Applications/ClustalOmega/my-out-seqs.fa', 'w') as currentFile:  #to set a test phrase so that it only progresses after clustal is doen
    with open(savefilename, 'w') as currentFile:
        currentFile.write('The clustal output is not done yet')

    if ordered == True:
        ClustalOCommandLine = bin + ' -i '+ MyInFiles + ' -o '+ MyOutFiles + ' -v --force --output-order=tree-order --outfmt=vie --resno --wrap=' \
                              + str(wrapLength)  #
    else:
        ClustalOCommandLine = bin + ' -i '+ MyInFiles + ' -o '+ MyOutFiles + ' -v --force --output-order=tree-order --outfmt=vie --resno --wrap=1000'
    ClustalOut = os.popen(ClustalOCommandLine)


    ClustalOutLines = []
    for line in ClustalOut:
        ClustalOutLines.append(line)

    send = False
    while send == False:
        # with open('/Applications/ClustalOmega/my-out-seqs.fa', 'r') as currentFile:  #using with for this automatically closes the file even if you crash
        with open(savefilename, 'r') as currentFile:
            testName = currentFile.readline()
            if testName != 'The clustal output is not done yet': send = True
            if send == False and len(ClustalOutLines) < 2:
                return False
                # print('first')

    os.remove(workingfilename)

    return savefilename