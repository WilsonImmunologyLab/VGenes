
import VGenesSQL
import VGenesSeq
from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
	setText
import re
import time
global FieldList
global RealNameList


FieldList = ['SeqName', 'SeqLen', 'GeneType', 'V1', 'V2', 'V3', 'D1', 'D2', 'D3', 'J1', 'J2', 'J3', 'StopCodon',
             'ReadingFrame', 'productive', 'Strand', 'VSeqend', 'VDJunction', 'Dregion', 'DJJunction', 'begJ',
             'VJunction', 'FR1From', 'FR1To', 'FR1length', 'FR1matches', 'FR1mis', 'FR1gaps', 'FR1PercentIdentity',
             'CDR1From', 'CDR1to', 'CDR1length', 'CDR1matches', 'CDR1mis', 'CDR1gaps', 'CDR1PercentIdentity', 'FR2From',
             'FR2To', 'FR2length', 'FR2matches', 'FR2mis', 'FR2gaps', 'FR2PercentIdentity', 'CDR2From', 'CDR2to',
             'CDR2length', 'CDR2matches', 'CDR2mis', 'CDR2gaps', 'CDR2PercentIdentity', 'FR3From', 'FR3To', 'FR3length',
             'FR3matches', 'FR3mis', 'FR3gaps', 'FR3PercentIdentity', 'TotMut', 'SeqAlignment', 'GVbeg', 'GVend',
             'GD1beg', 'GD1end', 'GD2beg', 'GD2end', 'GJbeg', 'GJend', 'Vbeg', 'Vend', 'D1beg', 'D1end', 'D2beg',
             'D2end', 'Jbeg', 'Jend', 'Project', 'Grouping', 'SubGroup', 'Species', 'Sequence', 'GermlineSequence',
             'CDR3DNA', 'CDR3AA', 'CDR3Length', 'CDR3beg', 'CDR3end', 'Specificity', 'Subspecificity', 'ClonalPool',
             'ClonalRank', 'VLocus', 'JLocus', 'DLocus', 'DateEntered', 'Comments', 'Quality', 'TotalMuts', 'Mutations',
             'IDEvent', 'CDR3MW', 'CDR3pI', 'Isotype', 'GCDR3beg', 'GCDR3end', 'Blank6', 'Blank7', 'Blank8', 'Blank9',
             'Blank10', 'Blank11', 'Blank12', 'Blank13', 'Blank14', 'Blank15', 'Blank16', 'Blank17', 'Blank18',
             'Blank19', 'Blank20', 'ID']

RealNameList = ["Name", "Length", "Type", "V gene", "V gene 2nd choice", "V gene 3rd choice", "D gene",
                " D gene 2nd choice", "D gene 3rd choice", "J gene", " J gene 2nd choice", "J gene 3rd choice",
                "Stop codons?", "Reading frame", "Productive?", "Strand", "End of V gene", "V to D Junction",
                "D region", "D to J junction", "Beginning of J", "V to J junction", "FWR1 first base", "FWR1 last base",
                "FWR1 length", "FWR1 matches", "FWR1 mismatches", "FWR1 gaps", "FWR1 percent identity",
                "CDR1 first base", "CDR1 last base", "CDR1 length", "CDR1 matches", "CDR1 mismatches", "CDR1 gaps",
                "CDR1 percent identity", "FWR2 first base", "FWR2 last base", "FWR2 length", "FWR2 matches",
                "FWR2 mismatches", "FWR2 gaps", "FWR2 percent identity", "CDR2 first base", "CDR2 last base",
                "CDR2 length", "CDR2 matches", "CDR2 mismatches", "CDR2 gaps", "CDR2 percent identity",
                "FWR3 first base", "FWR3 last base", "FWR3 length", "FWR3 matches", "FWR3 mismatches", "FWR3 gaps",
                "FWR3 percent identity", "IgBLAST mutation count", "IgBLAST Sequence Alignment", "Germline V begin",
                "Germline V end", "Germline D1 begin", "Germline D1 end", "Germline D2 begin", "Germline D2 end",
                "Germline J begin", "Germline J end", " V begin", " V end", " D1 begin", " D1 end", " D2 begin",
                " D2 end", " J begin", " J end", "Project", "Grouping", "Subgroup", "Species", "Sequence",
                " Germline sequence", "CDR3 DNA", "CDR3 peptide", "CDR3 length", "CDR3 first base", "CDR3 last base",
                "Specificity", "Subspecificity", "Clonal Pool", "Clonal Rank", "V locus", "J locus", "D locus",
                "Date and time entered", "Comments", "Quality", "Total Mutations", "Mutation list",
                "Insertions & deletions", "CDR3 molecular weight", "CDR3 isoelectric point", "Isotype",
                "Germlne CDR3 begin", "Germline CDR3 end", "Autoreactivity", "Blank7", "10xCluster", "Seuret_Cluster", "10xBarCode", "Population",
                "Blank12", "Blank13", "Blank14", "Blank15", "Blank16", "Blank17", "Blank18", "Blank19", "Blank20", "ID"]


def StandardReports(self, option, SequenceName, DBFilename):
    import os
    # first get list of seqs and info as tuple from DB:
    if option == 'FASTA Nucleotide file':
        fields = ['SeqName', 'Sequence']

        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1]
            Sequence = Sequence.replace('-', '')
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + Sequence + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)

    elif option == 'FASTA Nucleotide Rename file':
        ProjDict = {}
        i = 1
        fields = ['Project']  # create name dictionary
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        SQLStatement = 'SELECT DISTINCT' + SQLStatement[6:] + ' ORDER BY Project'
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        for item in DataIs:
            ProjName = item[0]
            ProjDict[ProjName] = 0


        fields = ['SeqName', 'Sequence', 'Project']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        FASTAFile = ''
        for item in DataIs:
            Seqname1 = item[0]
            # THE CLUSTAL OUTPUT IS NOT DONE
            Sequence = item[1]
            ProjName = item[2]
            ProjDict[ProjName] += 1

            Seqname = ProjName + '-' + str(ProjDict[ProjName])

            # validbases = 'agctAGCT'
            # Sequence =  Sequence1.join([char for char in Sequence1 if char in validbases])

            Sequence = Sequence.replace('-', '')
            # Sequence = Sequence.strip('-')
            checkit  = Sequence[0:11]


            if checkit != 'THE CLUSTAL':
                FASTAFile = FASTAFile + '>' + Seqname + '\n' + Sequence + '\n'
            else:
                print(Seqname1 + ' had THE CLUSTAL OUTPUT IS NOT DONE')

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)


    elif option == 'FASTA Amino Acid file':

        fields = ['SeqName', 'Sequence']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1]
            # Sequence  =
            AASeq, ErMessage = VGenesSeq.Translator(Sequence, 0)
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + AASeq + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)

    elif option == 'FASTA Germline Nucleotide file':

        fields = ['SeqName', 'GermlineSequence']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1]
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + Sequence + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)

    elif option == 'FASTA Germline Amino Acid file':

        fields = ['SeqName', 'GermlineSequence']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1]
            # Sequence  =
            AASeq, ErMessage = VGenesSeq.Translator(Sequence, 0)
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + AASeq + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)

    elif option == 'AbVec cloning PCR':
        CloningReport = ''
        CloningReportCSV = ''
        fields = ['SeqName', 'VLocus', 'JLocus', 'GeneType', 'Jend', 'Sequence']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        Heavy = []
        Kappa = []
        Lambda = []
        ErReport = ''
        # SeqNameCSV = ''
        for item in DataIs:
            SeqName = item[0]
            # if self.ui.ckReportCSV.isChecked():
            SeqNameCSV = SeqName

            Vlocus = item[1]
            Jlocus = item[2]
            Type = item[3]
            Jend = int(item[4])
            Sequence = item[5]
            Sequence.upper()

            Sequence = Sequence[0:Jend]

            Vlocus.upper()
            Jlocus.upper()
            Age1 = 0
            Sal1 = 0
            BsiW1 = 0
            Xho1 = 0
            if Type == 'Heavy':
                Age1 = Sequence.find('ACCGGT')
                Sal1 = Sequence.find('GTCGAC')
                if self.ui.ckReportText.isChecked():
                    if Age1 > 0:
                        SeqName += '\t**Age1 site detected at position ' + str(Age1) + '**'
                    if Sal1 > 0:
                        SeqName += '\t**Sal1 site detected at position ' + str(Sal1) + '**'

                if self.ui.ckReportCSV.isChecked():
                    if Age1 > 0:
                        SeqNameCSV += ',,,**Age1 site detected at position ' + str(Age1) + '**'
                    if Sal1 > 0:
                        SeqNameCSV += ',,,**Sal1 site detected at position ' + str(Sal1) + '**'

                if Vlocus[0:3] == 'VH1' or Vlocus[0:3] == 'VH5' or Vlocus[0:3] == 'VH7':
                    VPrimer = 'Age1-VH1/5'
                elif Vlocus == 'VH3-23':
                    VPrimer = 'Age1-VH3-23'
                elif Vlocus[0:6] == 'VH3-30' or Vlocus == 'VH3-33' or Vlocus == 'VH3-09':
                    VPrimer = 'Age1-VH3-30'
                elif Vlocus[0:3] == 'VH3':
                    VPrimer = 'Age1-VH3'
                elif Vlocus == 'VH4-34':
                    VPrimer = 'Age1-VH4-34'
                elif Vlocus[0:3] == 'VH4':
                    VPrimer = 'Age1-VH4'
                else:
                    VPrimer = 'None'

                if Jlocus[0:3] == 'JH1' or Jlocus[0:3] == 'JH2':
                    JPrimer = 'SalI-JH1/2'
                elif Jlocus == 'JH3':
                    JPrimer = 'SalI-JH3'
                elif Jlocus[0:3] == 'JH4' or Jlocus[0:3] == 'JH5':
                    JPrimer = 'JH4/5'
                elif Jlocus == 'JH6':
                    JPrimer = 'SalI-JH6'
                else:
                    JPrimer = 'None'

                Priming = (VPrimer, JPrimer)
                gene = (Priming, SeqName, SeqNameCSV)
                Heavy.append(gene)

            if Type == 'Kappa':
                Age1 = Sequence.find('ACCGGT')
                BsiW1 = Sequence.find('CGTACG')
                if self.ui.ckReportText.isChecked():
                    if Age1 > 0:
                        SeqName += '\t**Age1 site detected at position ' + str(Age1) + '**'
                    if BsiW1 > 0:
                        SeqName += '\t**BsiW1 site detected at position ' + str(BsiW1) + '**'
                if self.ui.ckReportCSV.isChecked():
                    if Age1 > 0:
                        SeqName += ',,,**Age1 site detected at position ' + str(Age1) + '**'
                    if BsiW1 > 0:
                        SeqName += ',,,**BsiW1 site detected at position ' + str(BsiW1) + '**'

                if Vlocus == 'VK1-09' or Vlocus == 'VK1-13' or Vlocus == 'VK1D-13':
                    VPrimer = 'Age1-Vk1-9'
                elif Vlocus == 'VK1-08' or Vlocus == 'VK1-43' or Vlocus == 'VK1D-43':
                    VPrimer = 'Age1-VK1D-43'
                elif Vlocus[0:3] == 'VK1':
                    VPrimer = 'Age1-VK1-5'
                elif Vlocus == 'VK2-28' or Vlocus == 'VK2D-28' or Vlocus == 'VK2-29' or Vlocus == 'VK2-30' or Vlocus == 'VK2D-30':
                    VPrimer = 'Age1-VK2-28'
                elif Vlocus[0:3] == 'Vk2':
                    VPrimer = 'Age1-VK2-24'
                elif Vlocus == 'VK3-11' or Vlocus == 'VK3D-11':
                    VPrimer = 'Age1-VK3-11'
                elif Vlocus == 'VK3-15' or Vlocus == 'VK3D-15':
                    VPrimer = 'Age1-VK3-15'
                elif Vlocus == 'VK3-20' or Vlocus == 'VK3D-20':
                    VPrimer = 'Age1-VK3-20'
                elif Vlocus[0:3] == 'Vk4':
                    VPrimer = 'Age1-VK4-1'

                if Jlocus[0:3] == 'JK1' or Jlocus[0:3] == 'JK2' or Jlocus[0:3] == 'JK4':
                    JPrimer = 'BsiW1-JK1/2/4'
                elif Jlocus == 'JK3':
                    JPrimer = 'BsiW1-JK3'
                elif Jlocus == 'JK5':
                    JPrimer = 'BsiW1-JK5'
                else:
                    JPrimer = 'None'

                Priming = (VPrimer, JPrimer)
                gene = (Priming, SeqName, SeqNameCSV)
                Kappa.append(gene)

            if Type == 'Lambda':
                Age1 = Sequence.find('ACCGGT')
                Xho1 = Sequence.find('CTCGAG')
                if self.ui.ckReportText.isChecked():
                    if Age1 > 0:
                        SeqName += '\t**Age1 site detected at position ' + str(Age1) + '**'

                    if Xho1 > 0:
                        SeqName += '\t**Xho1 site detected at position ' + str(Xho1) + '**'
                if self.ui.ckReportCSV.isChecked():
                    if Age1 > 0:
                        SeqNameCSV += ',,,**Age1 site detected at position ' + str(Age1) + '**'

                    if Xho1 > 0:
                        SeqNameCSV += ',,,**Xho1 site detected at position ' + str(Xho1) + '**'

                if Vlocus[0:3] == 'VL1':
                    VPrimer = 'Age1-VL1'
                elif Vlocus[0:3] == 'VL4' or Vlocus == 'VL5' or Vlocus == 'VL9':
                    VPrimer = 'Age1-VL4/5/9'
                elif Vlocus[0:3] == 'VL2':
                    VPrimer = 'Age1-VL2'
                elif Vlocus[0:3] == 'VL3':
                    VPrimer = 'Age1-VL3'
                elif Vlocus[0:3] == 'VL3':
                    VPrimer = 'Age1-VL3'
                elif Vlocus[0:3] == 'VL6':
                    VPrimer = 'Age1-VL6'
                elif Vlocus[0:3] == 'VL7' or Vlocus == 'VL8':
                    VPrimer = 'Age1-VL7/8'

                JPrimer = 'Xho1-JL'

                Priming = (VPrimer, JPrimer)
                gene = (Priming, SeqName, SeqNameCSV)
                Lambda.append(gene)

        if self.ui.ckReportText.isChecked():
            CloningReport = 'Cloning report for AbVec cloning PCR: ' + time.strftime('%c') + '\n'
        if self.ui.ckReportCSV.isChecked():
            CloningReportCSV = 'Cloning report for AbVec cloning PCR: ' + time.strftime('%c') + '\n'

        CurrentPrime = ('none', 'none')
        i = 0
        j = 0
        if len(Heavy) > 0:
            if self.ui.ckReportText.isChecked():
                CloningReport += '\n Heavy chain:'
            if self.ui.ckReportCSV.isChecked():
                CloningReportCSV += '\n Heavy chain: \n, Clone:, Sense:, Antisense:, Comments:'

            Heavy.sort()
            for gene in Heavy:
                SeqName = gene[1]
                SeqNameCSV = gene[2]
                Priming = gene[0]
                if Priming != CurrentPrime:
                    j += 1
                    i = 0
                    CurrentPrime = Priming
                    if self.ui.ckReportText.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',' + Priming[0] + ',' + Priming[1] + ':\n'
                    if self.ui.ckReportCSV.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',,' + Priming[0] + ',' + Priming[1] + ':\n'
                    if Priming[0] == "None" or Priming[1] == 'None':
                        NewPrimers = '*' + NewPrimers
                        ErReport = '*Either V or J primers could not be determined \n'
                    if self.ui.ckReportText.isChecked():
                        CloningReport += NewPrimers
                    if self.ui.ckReportCSV.isChecked():
                        CloningReportCSV += NewPrimers
                i += 1
                if self.ui.ckReportText.isChecked():
                    CloningReport += '\t' + str(i) + '. ' + SeqName + '\n'
                if self.ui.ckReportCSV.isChecked():
                    CloningReportCSV += str(i) + ',' + SeqNameCSV + '\n'

            CloningReport += ErReport
            CloningReportCSV += ErReport

        CurrentPrime = ('none', 'none')
        i = 0
        j = 0
        if len(Kappa) > 0:
            if self.ui.ckReportText.isChecked():
                CloningReport += '\n Kappa chain:'
            if self.ui.ckReportCSV.isChecked():
                CloningReportCSV += '\n Kappa chain: \n, Clone:, Sense:, Antisense:, Comments:'

            Kappa.sort()
            for gene in Kappa:
                SeqName = gene[1]
                SeqNameCSV = gene[2]
                Priming = gene[0]
                if Priming != CurrentPrime:
                    j += 1
                    i = 0
                    CurrentPrime = Priming
                    if self.ui.ckReportText.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',' + Priming[0] + ',' + Priming[1] + ':\n'
                    if self.ui.ckReportCSV.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',,' + Priming[0] + ',' + Priming[1] + ':\n'
                    if Priming[0] == "None" or Priming[1] == 'None':
                        NewPrimers = '*' + NewPrimers
                        ErReport = '*Either V or J primers could not be determined \n'
                    if self.ui.ckReportText.isChecked():
                        CloningReport += NewPrimers
                    if self.ui.ckReportCSV.isChecked():
                        CloningReportCSV += NewPrimers
                i += 1
                if self.ui.ckReportText.isChecked():
                    CloningReport += '\t' + str(i) + '. ' + SeqName + '\n'
                if self.ui.ckReportCSV.isChecked():
                    CloningReportCSV += str(i) + ',' + SeqNameCSV + '\n'

            CloningReport += ErReport
            CloningReportCSV += ErReport

        CurrentPrime = ('none', 'none')
        i = 0
        j = 0
        if len(Lambda) > 0:
            CloningReport += '\n Lambda chain: \n'
            CloningReportCSV += '\n Lambda chain: \n, Clone:, Sense:, Antisense:, Comments:\n'
            Lambda.sort()
            for gene in Lambda:
                SeqName = gene[1]
                SeqNameCSV = gene[2]
                Priming = gene[0]
                if Priming != CurrentPrime:
                    j += 1
                    i = 0
                    CurrentPrime = Priming
                    if self.ui.ckReportText.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',' + Priming[0] + ',' + Priming[1] + ':\n'
                    if self.ui.ckReportCSV.isChecked():
                        NewPrimers = '\n' + 'Master mix #' + str(j) + ',,' + Priming[0] + ',' + Priming[1] + ':\n'
                    if Priming[0] == "None" or Priming[1] == 'None':
                        NewPrimers = '*' + NewPrimers
                        ErReport = '*Either V or J primers could not be determined \n'
                    if self.ui.ckReportText.isChecked():
                        CloningReport += NewPrimers
                    if self.ui.ckReportCSV.isChecked():
                        CloningReportCSV += NewPrimers
                i += 1
                if self.ui.ckReportText.isChecked():
                    CloningReport += '\t' + str(i) + '. ' + SeqName + '\n'
                if self.ui.ckReportCSV.isChecked():
                    CloningReportCSV += str(i) + ',' + SeqNameCSV + '\n'

            CloningReport += ErReport
            CloningReportCSV += ErReport

        if self.ui.ckReportText.isChecked():
            Pathname = saveFile(self, 'text')
            if Pathname == None:
                return

            with open(Pathname, 'w') as currentfile:
                currentfile.write(CloningReport)

        if self.ui.ckReportDisplay.isChecked() and self.ui.ckReportText.isChecked():
            self.ShowVGenesText(Pathname)

        if self.ui.ckReportCSV.isChecked():
            Pathname = saveFile(self, 'csv')
            if Pathname == None:
                return

            with open(Pathname, 'w') as currentfile:
                currentfile.write(CloningReportCSV)

    elif option == 'HT-AbVec Cloning report':


        fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'Jend', 'Blank10']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        SQLStatement += ' ORDER BY Blank10'
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        NameList = []
        # for item in fields:
        #     NameList.append(str(self.TransLateFieldtoReal(item, False)))
        CSVOut = 'Sequence name, Sequence \n'


        preSeq = 'gcaactggtgtacattcc'  #'CTG CAA CCG GTG TAC ATT CA'


        for Seq in DataIs:
            answer = 'Yes'
            Barcode = Seq[8]
            Seqname1 = Seq[0]
            # fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'GermlineSequence', 'Blank10']
            SQLStatement = 'SELECT SeqName, GeneType, V1, J1, productive, Sequence, Vbeg, GermlineSequence, Blank10, IDEvent, Mutations FROM vgenesDB WHERE Blank10 = "'+ Barcode +'" ORDER BY GeneType'
            SeqClone = VGenesSQL.RunSQL(DBFilename, SQLStatement)

            if len(SeqClone) < 2:
                msg = 'For ' + Seqname1 + ' only a single variable gene was found with this barcode. You need both a heavy and light chain for each antibody'
                buttons = 'OK'

                answer = informationMessage(self, msg, buttons)
                answer = 'no'
            elif len(SeqClone) > 2:
                msg = 'For ' + Seqname1 + ' has more than one light chain associated with this sequence. Only one each H and L chanin should be designated with each barcode for this report to function.'
                buttons = 'OK'

                answer = informationMessage(self, msg, buttons)
                answer = 'no'


            if answer == 'Yes':


                for Record in SeqClone:
                    done = False
                    answer2 = 'Yes'
                    answer3 = 'Yes'
                    SeqName  = Record[0]
                    GeneType = Record[1]
                    Productive  = Record[4]
                    Sequence = Record[5]
                    Sequence = Sequence.upper()
                    Vbeg = int(Record[6])
                    GermSeq = Record[7]
                    #todo Figure GL for lambda and kappas and also fix Jend in database...
                    if GermSeq[len(GermSeq)-1] == 'G':
                        Jend = len(GermSeq)-1
                    else:
                        Jend = len(GermSeq)

                    IDevent = Record[9]
                    Mutations = Record[10]

                    # if Vbeg != 1:
                    #     msg = 'For ' + SeqName + ' the first base is not at position 1 possibly disrupting expression, continue with this sequence?'
                    #     buttons = 'YN'
                    #     answer2 = informationMessage(self, msg, buttons)

                    if Productive == 'No':
                        msg = SeqName + ' is a non-productive rearrangement and may not express, continue with this sequence?'
                        buttons = 'YN'
                        answer3 = informationMessage(self, msg, buttons)

                    if answer2 == 'Yes' and answer3 == 'Yes':
                        if IDevent == 'Insertion':
                            # msg = SeqName + 'Has an insertion and so needs to be prepared for synthesis by hand and will not be included.'
                            # buttons = 'OK'
                            # Sequence = ''
                            Mutset = Mutations.split(',')
                            for mutation in Mutset:
                                MutDetails = mutation.split('-')
                                for MutD in MutDetails:
                                    if MutD == 'Insertion':
                                        Insert = MutDetails[2]
                                        Jend += len(Insert)



                        Sequence = Sequence[:Jend]

                        if GeneType == 'Heavy':
                            postSeq = 'gcgtcgaccaagggcc' #'
                            HSeq = Sequence + postSeq

                        elif GeneType == 'Kappa':
                            postSeq = 'cgtacggtggcacagaaccggtgtccattcc'   #  CGT ACG gtg gca cag aAC CGG TGTcCATTCC 'gtacggtggctgcaccatctgtctt' # gtacggtggc'   #AA GAC AGA TGG TGC AGC CAC CGT ACG    CGTACGGTGGCTGCACCATCTGTCTT
                            LSeq = preSeq + Sequence + postSeq
                            done = True


                        elif GeneType == 'Lambda':
                            postSeq = 'ggtcagcccaaggccaaccccactgtcactctgttcccgccctcgaggtggcacagaaccggtgtccattcc'      #'ggtcagcccaaggctgccccctcggtcactctgttcccrccctcgagtgaggagcttcaagccaaca' #'ggtcagcccaaggctgccccctcggtcactctgttcccaccctcgagtgaggag'   #TG TTG GCT TGA AGC TCC TCA CTC GAG GGY GGG AAC AGA GTG  cactctgttcccrccctcgagtgaggagcttcaagccaaca
                            LSeq = preSeq + Sequence + postSeq
                            done = True

                        if done == True:
                            Sequence  = LSeq + HSeq
                            if IDevent == 'Deletion':
                                Sequence = Sequence.replace('-','')

                            SeqEntry = SeqName + ',' + Sequence + '\n'

                            CSVOut += SeqEntry





        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)
        #
    elif option == '10x Synthesis report':

        fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'Jend', 'Blank10']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        SQLStatement += ' ORDER BY Blank10'
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        NameList = []
        # for item in fields:
        #     NameList.append(str(self.TransLateFieldtoReal(item, False)))
        CSVOut = 'Sequence name, Sequence \n'

        preSeq = 'atcctttttctagtagcaactgcaaccggtgtacattca'  # 'CTG CAA CCG GTG TAC ATT CA'

        for Seq in DataIs:
            answer = 'Yes'
            Barcode = Seq[8]
            Seqname1 = Seq[0]
            # fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'GermlineSequence', 'Blank10']
            SQLStatement = 'SELECT SeqName, GeneType, V1, J1, productive, Sequence, Vbeg, GermlineSequence, Blank10, IDEvent, Mutations, Jend FROM vgenesDB WHERE Blank10 = "' + Barcode + '" ORDER BY GeneType'
            SeqClone = VGenesSQL.RunSQL(DBFilename, SQLStatement)
            if len(SeqClone) < 2:
                msg = 'For ' + Seqname1 + ' only a single variable gene was found with this barcode, continue with a single chain?'
                buttons = 'YN'
                answer = informationMessage(self, msg, buttons)

            # 'SELECT SeqName, GeneType, V1, J1, productive, Sequence, Vbeg, Jend, Blank10 FROM vgenesDB WHERE Blank10 = CACAAACCACGAAAGC-1 ORDER BY GeneType'
            # 'SELECT SeqName, GeneType, V1, J1, productive, Sequence, Vbeg, Jend, Blank10 FROM vgenesDB WHERE SeqName = "319_Cl100_H1" ORDER BY Blank10'

            if answer == 'Yes':

                for Record in SeqClone:
                    answer2 = 'Yes'
                    answer3 = 'Yes'
                    SeqName = Record[0]
                    GeneType = Record[1]
                    Productive = Record[4]
                    Sequence = Record[5]
                    Sequence = Sequence.upper()
                    Vbeg = int(Record[6])
                    GermSeq = Record[7]
                    Jend = Record[11]
                    # todo Figure GL for lambda and kappas and also fix Jend in database...
                    # if GermSeq[len(GermSeq) - 1] == 'G':
                    #     Jend = len(GermSeq) - 1
                    # else:
                    #     Jend = len(GermSeq)

                    IDevent = Record[9]
                    Mutations = Record[10]

                    if Vbeg != 1:
                        msg = 'For ' + SeqName + ' the first base is not at position 1 possibly disrupting expression, continue with this sequence?'
                        buttons = 'YN'
                        answer2 = informationMessage(self, msg, buttons)

                    if Productive == 'No':
                        msg = SeqName + ' is a non-productive rearrangement and may not express, continue with this sequence?'
                        buttons = 'YN'
                        answer3 = informationMessage(self, msg, buttons)

                    if answer2 == 'Yes' and answer3 == 'Yes':
                        if IDevent == 'Insertion':
                            # msg = SeqName + 'Has an insertion and so needs to be prepared for synthesis by hand and will not be included.'
                            # buttons = 'OK'
                            # Sequence = ''
                            Mutset = Mutations.split(',')
                            for mutation in Mutset:
                                MutDetails = mutation.split('-')
                                for MutD in MutDetails:
                                    if MutD == 'Insertion':
                                        Insert = MutDetails[2]
                                        Jend += len(Insert)

                        Sequence = Sequence[:Jend]

                        if GeneType == 'Heavy':
                            postSeq = 'gcgtcgaccaagggcccatcggtcttcc'  # G GAA GAC CGA TGG GCC CTT GGT CGA CGC     GCGTCGACCAAGGGCCCATCGGTCTTCC
                        elif GeneType == 'Kappa':
                            postSeq = 'cgtacggtggctgcaccatctgtctt'  # gtacggtggc'   #AA GAC AGA TGG TGC AGC CAC CGT ACG    CGTACGGTGGCTGCACCATCTGTCTT
                        elif GeneType == 'Lambda':
                            postSeq = 'ggtcagcccaaggctgccccctcggtcactctgttcccrccctcgagtgaggagcttcaagccaaca'  # 'ggtcagcccaaggctgccccctcggtcactctgttcccaccctcgagtgaggag'   #TG TTG GCT TGA AGC TCC TCA CTC GAG GGY GGG AAC AGA GTG  cactctgttcccrccctcgagtgaggagcttcaagccaaca

                        Sequence = preSeq + Sequence + postSeq
                        if IDevent == 'Deletion':
                            Sequence = Sequence.replace('-', '')

                        SeqEntry = SeqName + ',' + Sequence + '\n'

                        CSVOut += SeqEntry

        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)



    elif option == 'Sequence summary':
        fields = ['SeqName', 'Project', 'V1', 'D1', 'J1', 'VLocus', 'JLocus', 'productive', 'TotMut', 'CDR3DNA', 'CDR3AA',
                  'CDR3Length', 'CDR3pI', 'ClonalPool', 'Isotype', 'Sequence', 'Blank7']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        NameList = []
        for item in fields:
            NameList.append(str(self.TransLateFieldtoReal(item, False)))
        CSVOut = ''
        i = 0
        for item in NameList:

            CSVOut += str(item)
            CSVOut += ','


        CSVOut += '\n'

        for Record in DataIs:

            for item in Record:
                StringItem = str(item)
                CSVOut += str(StringItem)
                CSVOut += ','

            i+=1

            CSVOut += '\n'


        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)
        #
        # SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        # DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        # NameList = []
        # for item in fields:
        #     NameList.append(str(self.TransLateFieldtoReal(item, False)))
        #
        # if self.ui.ckReportText.isChecked():
        #     SeqSumm = ''
        #     for record in DataIs:
        #         i = 0
        #         for item in record:
        #             SeqSumm = SeqSumm + NameList[i] + ': \t' + str(item) + '\n'
        #             i += 1
        #         SeqSumm += '\n'
        #
        #     Pathname = saveFile(self, 'text')
        #
        #     with open(Pathname, 'w') as currentfile:
        #         currentfile.write(SeqSumm)
        #
        # if self.ui.ckReportDisplay.isChecked() and self.ui.ckReportText.isChecked():
        #     self.ShowVGenesText(Pathname)
        #
        # if self.ui.ckReportCSV.isChecked():
        #     SeqSumm = ''
        #     for item in NameList:
        #         SeqSumm = SeqSumm + str(item) + ','
        #     SeqSumm += '\n'
        #     for record in DataIs:
        #         for item in record:
        #             SeqSumm = SeqSumm + str(item) + ','
        #         SeqSumm += '\n'
        #
        #     Pathname = saveFile(self, 'csv')
        #
        #     with open(Pathname, 'w') as currentfile:
        #         currentfile.write(SeqSumm)






    elif option == 'Comma seperated values (.csv)':

        fields = 'All'
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        (dirname, filename) = os.path.split(DBFilename)
        filename = filename[:(len(filename)-4)]
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        CSVOut = ''
        Vfam = ''
        i = 0
        for item in RealNameList:
            if i != 58:
                if str(item) == 'Blank8':
                    item = 'Cluster'
                if str(item) == 'Blank9':
                    item = 'Seuret Cluster'
                if str(item) == 'Blank10':
                    item = 'Barcode'
                if str(item) == 'Blank11':
                    item = 'Population'

                CSVOut += str(item)
                CSVOut += ','
            i += 1
        CSVOut += 'V-family'
        CSVOut += ',Subject'

        CSVOut += '\n'
        for Record in DataIs:
            i = 0
            for item in Record:
                if i != 58:
                    if i == 90:
                        Vfam = item[:3]

                    StringItem = str(item)
                    if i == 97:
                        StringItem = StringItem.replace(',', '/').replace('\n', 'linefeed')
                    CSVOut += str(StringItem)
                    CSVOut += ','
                i += 1
            CSVOut += Vfam
            CSVOut += ','+filename
            CSVOut += '\n'

        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)

    elif option == 'CSV format Entire VDB':

        SQLSTATEMENT = 'SELECT Field,FieldNickName from fieldsname ORDER BY ID'
        DataIn = VGenesSQL.RunSQL(DBFilename, SQLSTATEMENT)
        fields = [i[0] for i in DataIn]
        field_names = [i[1] for i in DataIn]

        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        (dirname, filename) = os.path.split(DBFilename)
        filename = filename[:(len(filename)-4)]
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        CSVOut = ''

        CSVOut += ','.join(fields) + '\n'
        CSVOut += ','.join(field_names) + '\n'
        for record in DataIs:
            record_new = [str(x) for x in record]
            record_new[58] = re.sub(r'\n','#',record_new[58])
            record_new[97] = re.sub(',', '|', record_new[97])
            CSVOut += ','.join(record_new) + '\n'

        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)

    elif option == 'Clonal Analysis (.csv)':
        fields = 'All'
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        CSVOut = ''
        i = 0
        for item in RealNameList:
            if i != 58:
                CSVOut += str(item)
                CSVOut += ','
            i += 1
        CSVOut += '\n'
        for Record in DataIs:
            i = 0
            for item in Record:
                if i != 58:
                    StringItem = str(item)
                    if i == 97:
                        StringItem = StringItem.replace(',', '/').replace('\n', 'linefeed')
                    CSVOut += str(StringItem)
                    CSVOut += ','
                i += 1
            CSVOut += '\n'

        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)

    elif option == 'Custom report':

        print('custom report generator')

    print('done')