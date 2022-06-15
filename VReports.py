
import VGenesSQL
import VGenesSeq
from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
    setText
from VGenesMain import ProgressBar
from VGenesMain import GibsonDialog, PatentDialog, ExportOptionDialog
from PyQt5.QtWidgets import QMessageBox, QAbstractItemView, QTableWidgetItem, QTableWidget, QHeaderView, QTextEdit, QLineEdit, QCheckBox, QComboBox
from PyQt5.QtCore import QThread, pyqtSignal, Qt
from PyQt5 import QtGui

import re
import os
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


class CSVRep_thread(QThread):
    loadProgress = pyqtSignal(int, str)
    trigger = pyqtSignal(list)

    def __int__(self):
        super(CSVRep_thread, self).__init__()
        self.DBFilename = ''
        self.SequenceName = ''
        self.Pathname = ''
        self.vgene = ''

    def run(self):
        DBFilename = self.DBFilename
        SequenceName = self.SequenceName
        Pathname = self.Pathname
        Vgenes = self.vgene

        pct = 0
        label = "Fetching records from DB: " + DBFilename
        self.loadProgress.emit(pct, label)

        fields = 'All'
        SQLStatement = VGenesSQL.MakeSQLStatement(Vgenes, fields, SequenceName)
        (dirname, filename) = os.path.split(DBFilename)
        filename = filename[:(len(filename) - 4)]
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
        process = 1
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
            CSVOut += ',' + filename
            CSVOut += '\n'

            pct = int(process / len(DataIs) * 100)
            label = "Fetching records: " + str(process) + '/' + str(len(DataIs))
            self.loadProgress.emit(pct, label)
            process += 1

        pct = 100
        label = "Writing to file: " + Pathname
        self.loadProgress.emit(pct, label)
        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)

        Msg = 'Comma seperated values Report generated!'
        self.trigger.emit([0, Msg])

class CSVRepAll_thread(QThread):
    loadProgress = pyqtSignal(int, str)
    trigger = pyqtSignal(list)

    def __int__(self):
        super(CSVRepAll_thread, self).__init__()
        self.DBFilename = ''
        self.SequenceName = ''
        self.Pathname = ''
        self.vgene = ''
        self.list = []

    def run(self):
        DBFilename = self.DBFilename
        SequenceName = self.SequenceName
        Pathname = self.Pathname
        Vgenes = self.vgene
        list = self.list

        pct = 0
        label = "Fetching records from DB: " + DBFilename
        self.loadProgress.emit(pct, label)

        SQLSTATEMENT = 'SELECT Field,FieldNickName from fieldsname ORDER BY ID'
        DataIn = VGenesSQL.RunSQL(DBFilename, SQLSTATEMENT)
        fields = [i[0] for i in DataIn]
        field_names = [i[1] for i in DataIn]

        #SQLStatement = VGenesSQL.MakeSQLStatementNew(Vgenes, fields, SequenceName)
        if len(list) == 0:
            WHEREStatement = ' WHERE 1'
        else:
            WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(list) + '")'
        SQLStatement = 'SELECT * from vgenesDB' + WHEREStatement
        (dirname, filename) = os.path.split(DBFilename)
        filename = filename[:(len(filename) - 4)]
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        CSVOut = ''

        CSVOut += ','.join(fields) + '\n'
        CSVOut += ','.join(field_names) + '\n'
        process = 1
        for record in DataIs:
            record_new = [str(x) for x in record]
            record_new[58] = re.sub(r'\n', '#', record_new[58])
            record_new[97] = re.sub(',', '|', record_new[97])
            CSVOut += ','.join(record_new) + '\n'

            pct = int(process / len(DataIs) * 100)
            label = "Fetching records: " + str(process) + '/' + str(len(DataIs))
            self.loadProgress.emit(pct, label)
            process += 1

        pct = 100
        label = "Writing to file: " + Pathname
        self.loadProgress.emit(pct, label)

        with open(Pathname, 'w') as currentfile:
            currentfile.write(CSVOut)

        Msg = 'Comma seperated values Report generated!'
        self.trigger.emit([0, Msg])

class HCLC_thread(QThread):
    HCLC_progress = pyqtSignal(int, int, int)
    HCLC_finish = pyqtSignal(list)

    def __int__(self):
        super(HCLC_thread, self).__init__()
        self.DBFilename = ''
        self.Pathname = ''
        self.checkRecords = []

    def run(self):
        Msg = ''
        sign = ''
        if len(self.checkRecords) == 0:
            SQLStatement = 'SELECT DISTINCT(Blank10) FROM vgenesdb WHERE GeneType = "Heavy"'
            DataIn = VGenesSQL.RunSQL(self.DBFilename, SQLStatement)
            if len(DataIn) < 2:
                Msg = 'Your VGene DB do not have any barcode information!'
                sign = 1
            else:
                CSVOut = ''
                # make CSV header
                SQLSTATEMENT = 'SELECT Field,FieldNickName from fieldsname ORDER BY ID'
                DataInHeader = VGenesSQL.RunSQL(self.DBFilename, SQLSTATEMENT)
                fields = [i[0] for i in DataInHeader]
                field_names = ['HC_' + i[1] for i in DataInHeader] + ['LC_' + i[1] for i in DataInHeader]
                CSVOut += ','.join(field_names) + '\n'

                seq_num = 0
                progress = 0
                for record in DataIn:
                    barcode = record[0]
                    if barcode == 'Blank10' or barcode == '':
                        pass
                    else:
                        SQLStatement1 = 'SELECT * FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType = "Heavy"'
                        DataIn1 = VGenesSQL.RunSQL(self.DBFilename, SQLStatement1)

                        SQLStatement2 = 'SELECT ' + ",".join(fields) + ' FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType IN ("Kappa","Lambda")'
                        DataIn2 = VGenesSQL.RunSQL(self.DBFilename, SQLStatement2)

                        if len(DataIn1) == 1 and len(DataIn2) == 1:
                            data_hc = [str(x) for x in DataIn1[0]]
                            data_hc[58] = re.sub(r'\n', '#', data_hc[58])
                            data_hc[97] = re.sub(',', '|', data_hc[97])
                            data_lc = [str(x) for x in DataIn2[0]]
                            data_lc[58] = re.sub(r'\n', '#', data_lc[58])
                            data_lc[97] = re.sub(',', '|', data_lc[97])

                            CSVOut += ','.join(data_hc) + ',' + ','.join(data_lc) + '\n'
                            seq_num += 1

                    self.HCLC_progress.emit(progress, len(DataIn), int(progress/len(DataIn)*100))
                    progress += 1

                if seq_num > 0:
                    with open(self.Pathname, 'w') as currentfile:
                        currentfile.write(CSVOut)

                    Msg = 'Total ' + str(seq_num) + ' HC/LC pairs were found and exported!'
                    sign = 0
                else:
                    Msg = 'Did not find any HC/LC pair in your current DB!'
                    sign = 1
        else:
            list_str = '("' + '","'.join(self.checkRecords) + '")'
            SQLStatement = 'SELECT DISTINCT(Blank10) FROM vgenesdb WHERE SeqName IN ' + list_str
            DataIn = VGenesSQL.RunSQL(self.DBFilename, SQLStatement)
            if len(DataIn) < 2:
                Msg = 'Your VGene DB do not have any barcode information!'
                sign = 1
            else:
                CSVOut = ''
                # make CSV header
                SQLSTATEMENT = 'SELECT Field,FieldNickName from fieldsname ORDER BY ID'
                DataInHeader = VGenesSQL.RunSQL(self.DBFilename, SQLSTATEMENT)
                field_names = ['HC_' + i[1] for i in DataInHeader] + ['LC_' + i[1] for i in DataInHeader]
                CSVOut += ','.join(field_names) + '\n'

                seq_num = 0
                progress = 0
                for record in DataIn:
                    barcode = record[0]
                    if barcode == 'Blank10' or barcode == '':
                        pass
                    else:
                        SQLStatement1 = 'SELECT * FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType = "Heavy"'
                        DataIn1 = VGenesSQL.RunSQL(self.DBFilename, SQLStatement1)

                        SQLStatement2 = 'SELECT * FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType IN ("Kappa","Lambda")'
                        DataIn2 = VGenesSQL.RunSQL(self.DBFilename, SQLStatement2)

                        if len(DataIn1) == 1 and len(DataIn2) == 1:
                            data_hc = [str(x) for x in DataIn1[0]]
                            data_hc[58] = re.sub(r'\n', '#', data_hc[58])
                            data_hc[97] = re.sub(',', '|', data_hc[97])
                            data_lc = [str(x) for x in DataIn2[0]]
                            data_lc[58] = re.sub(r'\n', '#', data_lc[58])
                            data_lc[97] = re.sub(',', '|', data_lc[97])

                            CSVOut += ','.join(data_hc) + ',' + ','.join(data_lc) + '\n'
                            seq_num += 1

                    self.HCLC_progress.emit(progress, len(DataIn), int(progress / len(DataIn) * 100))
                    progress += 1

                if seq_num > 0:
                    with open(self.Pathname, 'w') as currentfile:
                        currentfile.write(CSVOut)

                    Msg = 'Total ' + str(seq_num) + ' HC/LC pairs were found and exported!'
                    sign = 0
                else:
                    Msg = 'Did not find any HC/LC pair in your current DB!'
                    sign = 1
        self.HCLC_finish.emit([sign, Msg])

def StandardReports(self, option, SequenceName, DBFilename):
    import os
    # first get list of seqs and info as tuple from DB:
    if option == 'FASTA Nucleotide file':
        '''
        fields = ['SeqName', 'Sequence']

        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        '''
        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,Sequence FROM vgenesdb' + WHEREStatement
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

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,Sequence FROM vgenesdb' + WHEREStatement
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1]
            Sequence = Sequence.replace('-', '')
            AASeq, ErMessage = VGenesSeq.Translator(Sequence, 0)
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + AASeq + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)
    elif option == 'FASTA Germline Nucleotide file':

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,Sequence FROM vgenesdb' + WHEREStatement
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1].upper()
            Sequence = Sequence.replace('-', '')
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + Sequence + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)
    elif option == 'FASTA Germline Amino Acid file':

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,Sequence FROM vgenesdb' + WHEREStatement
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        FASTAFile = ''
        for item in DataIs:
            Seqname = item[0]
            Sequence = item[1].upper()
            Sequence = Sequence.replace('-', '')
            AASeq, ErMessage = VGenesSeq.Translator(Sequence, 0)
            FASTAFile = FASTAFile + '>' + Seqname + '\n' + AASeq + '\n'

        Pathname = saveFile(self, 'FASTA')
        if Pathname == None:
            return

        with open(Pathname, 'w') as currentfile:
            currentfile.write(FASTAFile)

        self.ShowVGenesText(Pathname)
    elif option == 'FASTA raw Nucleotide sequences':
        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,Blank20 FROM vgenesdb' + WHEREStatement
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
    elif option == 'AbVec cloning PCR':
        CloningReport = ''
        CloningReportCSV = ''
        '''
        fields = ['SeqName', 'VLocus', 'JLocus', 'GeneType', 'Jend', 'Sequence']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        '''
        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,VLocus,JLocus,GeneType,Jend,Sequence FROM vgenesdb' + WHEREStatement
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
                else:
                    VPrimer = 'None'

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
                elif Vlocus[0:3] == 'VL4' or Vlocus[0:3] == 'VL5' or Vlocus[0:3] == 'VL9':
                    VPrimer = 'Age1-VL4/5/9'
                elif Vlocus[0:3] == 'VL2':
                    VPrimer = 'Age1-VL2'
                elif Vlocus[0:3] == 'VL3':
                    VPrimer = 'Age1-VL3'
                elif Vlocus[0:3] == 'VL3':
                    VPrimer = 'Age1-VL3'
                elif Vlocus[0:3] == 'VL6':
                    VPrimer = 'Age1-VL6'
                elif Vlocus[0:3] == 'VL7' or Vlocus[0:3] == 'VL8':
                    VPrimer = 'Age1-VL7/8'
                else:
                    VPrimer = 'None'
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

            Msg = 'AbVec Cloning Report generated!'
            self.ShowMessageBox([0, Msg])
    elif option == 'HT-AbVec Cloning report':
        '''
        fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', ',Jend', 'Blank10']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        SQLStatement += ' ORDER BY Blank10'
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        '''

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,GeneType,V1,J1,productive,Sequence,Vbeg,Jend,Blank10 FROM vgenesdb' + WHEREStatement
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
                    Sequence = Record[5].upper()
                    Vbeg = int(Record[6])
                    GermSeq = Record[7].upper()
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

        Msg = 'HT-AbVec Cloning Report generated!'
        self.ShowMessageBox([0, Msg])
    elif option == '10x Synthesis report':

        '''
        fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'Jend', 'Blank10']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        SQLStatement += ' ORDER BY Blank10'
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        '''

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,GeneType,V1,J1,productive,Sequence,Vbeg,Jend,Blank10 FROM vgenesdb' + WHEREStatement
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
                    Sequence = Record[5].upper()
                    Vbeg = int(Record[6])
                    GermSeq = Record[7].upper()
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

        Msg = '10x Synthesis Report generated!'
        self.ShowMessageBox([0, Msg])
    elif option == 'Sequence summary':
        '''
        fields = ['SeqName', 'Project', 'V1', 'D1', 'J1', 'VLocus', 'JLocus', 'productive', 'TotMut', 'CDR3DNA', 'CDR3AA',
                  'CDR3Length', 'CDR3pI', 'ClonalPool', 'Isotype', 'Sequence', 'Blank7']
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        '''

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        fields = ['SeqName', 'Project', 'V1', 'D1', 'J1', 'VLocus', 'JLocus', 'productive', 'TotMut', 'CDR3DNA',
                  'CDR3AA', 'CDR3Length', 'CDR3pI', 'ClonalPool', 'Isotype', 'Sequence', 'Blank7']
        SQLStatement = 'SELECT ' + ','.join(fields) + ' FROM vgenesdb' + WHEREStatement
        SQLStatement += ' ORDER BY Blank10'
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

        Msg = 'Sequence summary Report generated!'
        self.ShowMessageBox([0, Msg])
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
        '''
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

        Msg = 'Comma seperated values Report generated!'
        self.ShowMessageBox([0, Msg])
        '''
        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        # try multi-thread
        self.workThread = CSVRep_thread(self)
        self.workThread.DBFilename = DBFilename
        self.workThread.SequenceName = SequenceName
        self.workThread.Pathname = Pathname
        self.workThread.vgene = self
        self.workThread.start()
        self.workThread.trigger.connect(self.ShowMessageBox)
        self.workThread.loadProgress.connect(self.progressLabel)

        self.progress = ProgressBar(self)
        self.progress.show()
    elif option == 'CSV format full records':
        '''
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

        Msg = 'CSV format Entire VDB Report generated!'
        self.ShowMessageBox([0, Msg])
        '''

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'No':
                return

        Pathname = saveFile(self, 'csv')
        if Pathname == None:
            return

        # try multi-thread
        self.workThread = CSVRepAll_thread(self)
        self.workThread.DBFilename = DBFilename
        self.workThread.SequenceName = SequenceName
        self.workThread.Pathname = Pathname
        self.workThread.list = selected_list
        self.workThread.vgene = self
        self.workThread.start()
        self.workThread.trigger.connect(self.ShowMessageBox)
        self.workThread.loadProgress.connect(self.progressLabel)

        self.progress = ProgressBar(self)
        self.progress.show()
    elif option == 'Clonal Analysis (.csv)':
        '''
        fields = 'All'
        SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
        DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        '''

        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT * FROM vgenesdb' + WHEREStatement
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

        Msg = 'Clonal Analysis Report generated!'
        self.ShowMessageBox([0,Msg])
    elif option == 'Heavy/Light Chain pairs (.csv)':
        Pathname = saveFile(self.parent(), 'csv')
        if Pathname == None:
            return

        if self.ui.checkBoxAll.isChecked():
            listItems = []
        else:
            listItems = self.CheckedRecords
            if len(listItems) == 0:
                pass
            else:
                mode = 1
                msg = 'You selected part of records, will only identify HC/LC pairs for your selected records!'
                QMessageBox.information(self, 'Information', msg, QMessageBox.Ok, QMessageBox.Ok)

        self.HCLC_Thread = HCLC_thread(self)
        self.HCLC_Thread.DBFilename = DBFilename
        self.HCLC_Thread.Pathname = Pathname
        self.HCLC_Thread.checkRecords = listItems
        self.HCLC_Thread.HCLC_progress.connect(self.result_display)
        self.HCLC_Thread.HCLC_finish.connect(self.ShowMessageBox)
        self.HCLC_Thread.start()

        self.progress = ProgressBar(self)
        self.progress.show()

        '''
        SQLStatement = 'SELECT DISTINCT(Blank10) FROM vgenesdb WHERE GeneType = "Heavy"'
        DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        if len(DataIn) < 2:
            Msg = 'Your VGene DB do not have any barcode information!'
            QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok,
                                    QMessageBox.Ok)
            return

        CSVOut = ''
        # make CSV header
        SQLSTATEMENT = 'SELECT Field,FieldNickName from fieldsname ORDER BY ID'
        DataInHeader = VGenesSQL.RunSQL(DBFilename, SQLSTATEMENT)
        field_names = ['HC_' + i[1] for i in DataInHeader] + ['LC_' + i[1] for i in DataInHeader]
        CSVOut += ','.join(field_names) + '\n'

        seq_num = 0
        for record in DataIn:
            barcode = record[0]
            if barcode == 'Blank10' or barcode == '':
                pass
            else:
                SQLStatement1 = 'SELECT * FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType = "Heavy"'
                DataIn1 = VGenesSQL.RunSQL(DBFilename, SQLStatement1)

                SQLStatement2 = 'SELECT * FROM vgenesdb WHERE Blank10 = "' + barcode + '" AND GeneType IN ("Kappa","Lambda")'
                DataIn2 = VGenesSQL.RunSQL(DBFilename, SQLStatement2)

                if len(DataIn1) == 1 and len(DataIn2) == 1:
                    data_hc = [str(x) for x in DataIn1[0]]
                    data_hc[58] = re.sub(r'\n', '#', data_hc[58])
                    data_hc[97] = re.sub(',', '|', data_hc[97])
                    data_lc = [str(x) for x in DataIn2[0]]
                    data_lc[58] = re.sub(r'\n', '#', data_lc[58])
                    data_lc[97] = re.sub(',', '|', data_lc[97])

                    CSVOut += ','.join(data_hc) + ',' + ','.join(data_lc) + '\n'
                    seq_num += 1

        if seq_num > 0:
            Pathname = saveFile(self, 'csv')
            if Pathname == None:
                return
            with open(Pathname, 'w') as currentfile:
                currentfile.write(CSVOut)

            Msg = 'Total ' + str(seq_num) + ' HC/LC pairs were found and exported!'
            QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok,
                                QMessageBox.Ok)
        else:
            Msg = 'Did not find any HC/LC pair in your current DB!'
            QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok,
                                QMessageBox.Ok)
            return
        '''
    elif option == 'Custom report':
        pass
    elif option == 'Count AA mutations':
        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        SQLStatement = 'SELECT SeqName,GeneType,SeqAlignment FROM vgenesdb' + WHEREStatement

        DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
        if len(DataIn) > 0:

            Pathname = saveFile(self.parent(), 'csv')
            if Pathname == None:
                return
            with open(Pathname, 'w') as currentfile:
                out_str = 'SeqName,GeneType,Num of AA mutations of V\n'
                currentfile.write(out_str)
                for records in DataIn:
                    name = records[0]
                    type = records[1]
                    alignment = records[2]

                    ref_aa = ''
                    cur_aa = ''

                    lines = alignment.split('\n')

                    for i in range(len(lines)):
                        cur_line = lines[i]
                        if 'Query' in cur_line:
                            ref_aa += lines[i - 1]
                            try:
                                if lines[i + 2].find('          ',0,10) != -1:
                                    cur_aa += lines[i + 2]
                            except:
                                pass

                    cur_aa = cur_aa.rstrip()
                    ref_aa = ref_aa[0:len(cur_aa)]

                    num_aa_mutation_v = 0
                    for i in range(len(ref_aa)):
                        if ref_aa[i] != cur_aa[i]:
                            num_aa_mutation_v += 1
                    num_aa_mutation_j = 0
                    '''
                    if type != 'Heavy':
                        SQLStatement = 'SELECT Jbeg,Jend,Sequence,GermlineSequence FROM vgenesdb WHERE SeqName = "' + name +'"'
                        res = VGenesSQL.RunSQL(DBFilename, SQLStatement)
                        res = res[0]
                        j_start = int(res[0]) - 1
                        j_end = int(res[1])
                        ref_seq = res[3]
                        cur_seq = res[2]
                        ref_nt = ref_seq[j_start:j_end]
                        cur_nt = cur_seq[j_start:j_end]

                        # determine ORF
                        ruler_line = lines[2]
                        aa_line = lines[3]
                        nt_line = lines[4]

                        pos = ruler_line.find('<')
                        aa_line = aa_line[pos:]
                        nt_line = nt_line[pos:]

                        num_space = 0
                        while num_space < 10:
                            if aa_line[num_space] == ' ':
                                num_space += 1
                            else:
                                break
                        orf_offset = num_space % 3
                        orf_offset = orf_offset - 1
                        if orf_offset < 0:
                            orf_offset = orf_offset + 3

                        orf = j_start % 3
                        orf = orf_offset - orf
                        if orf < 0:
                            orf = orf + 3

                        #ref_1 = VGenesSeq.Translator(ref_nt, 0)[0]
                        #ref_2 = VGenesSeq.Translator(ref_nt, 1)[0]
                        #ref_3 = VGenesSeq.Translator(ref_nt, 2)[0]

                        #ref_11 = VGenesSeq.Translator(ref_seq, 0)[0]
                        #ref_21 = VGenesSeq.Translator(ref_seq, 1)[0]
                        #ref_31 = VGenesSeq.Translator(ref_seq, 2)[0]

                        ref_aa, msg1 = VGenesSeq.Translator(ref_nt, orf)
                        cur_aa, msg2 = VGenesSeq.Translator(cur_nt, orf)

                        for i in range(len(ref_aa)):
                            if ref_aa[i] != '~':
                                if ref_aa[i] != cur_aa[i]:
                                    num_aa_mutation_j += 1
                    '''
                    num_aa_mutation = num_aa_mutation_v + num_aa_mutation_j

                    #out_str = name + ',' + type + ',' + str(num_aa_mutation) + ',' + str(num_aa_mutation_v) + ',' + str(num_aa_mutation_j) + '\n'
                    out_str = name + ',' + type + ',' + str(num_aa_mutation) + '\n'
                    currentfile.write(out_str)
                    #res.append((name, type, num_aa_mutation, num_aa_mutation_v, num_aa_mutation_j))
            self.ShowVGenesText(Pathname)
    elif option == 'Sequence for GibsonClone':
        if len(self.AntibodyCandidates) == 0:
            selected_list = self.CheckedRecords
        else:
            selected_list = self.AntibodyCandidates

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        self.myGibsonDialog = GibsonDialog()

        SQLStatement = 'SELECT SeqName,GeneType,Sequence,Vbeg,Jend,Blank7,SeqAlignment FROM vgenesdb' + WHEREStatement
        DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

        horizontalHeader = ['Seq check', 'Name', 'GeneType', 'J end', 'V(D)J sequence', 'Translated AA']
        num_row = len(DataIn)
        num_col = len(horizontalHeader)
        self.myGibsonDialog.ui.tableWidget.setRowCount(num_row)
        self.myGibsonDialog.ui.tableWidget.setColumnCount(num_col)
        self.myGibsonDialog.ui.tableWidget.setHorizontalHeaderLabels(horizontalHeader)
        self.myGibsonDialog.ui.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.myGibsonDialog.ui.tableWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.myGibsonDialog.ui.tableWidget.setSelectionBehavior(QAbstractItemView.SelectRows)

        if len(DataIn) > 0:
            index = 0
            for records in DataIn:
                SeqName = records[0]
                GeneType = records[1]
                Sequence = records[2].upper()
                Vbeg = int(records[3])
                Jend = int(records[4])

                VDJseq = Sequence[Vbeg - 1:Jend]
                VDJseq = re.sub(r'\W+','',VDJseq)   # remove all alignment gaps

                JendSeq = VDJseq[-6:]
                checkRes = checkJend(GeneType, JendSeq)
                try:
                    ORF = int(records[5])
                except:
                    ORF = getORF(records[6])
                if ORF != 0:
                    checkRes = 'Check ORF'
                    VDJseq = VDJseq[ORF:]
                AAseq, msg = VGenesSeq.Translator(VDJseq, 0)
                if "*" in AAseq:
                    checkRes = 'ORF error'

                unit1 = QTableWidgetItem(checkRes)
                unit2 = QTableWidgetItem(SeqName)
                unit3 = QTableWidgetItem(GeneType)
                unit4 = QTableWidgetItem(JendSeq)
                #unit5 = QTableWidgetItem(VDJseq)
                unit1.setFlags(Qt.ItemIsEnabled)
                unit2.setFlags(Qt.ItemIsEnabled)
                unit3.setFlags(Qt.ItemIsEnabled)
                unit4.setFlags(Qt.ItemIsEnabled)
                self.myGibsonDialog.ui.tableWidget.setItem(index, 0, unit1)
                self.myGibsonDialog.ui.tableWidget.setItem(index, 1, unit2)
                self.myGibsonDialog.ui.tableWidget.setItem(index, 2, unit3)
                self.myGibsonDialog.ui.tableWidget.setItem(index, 3, unit4)
                #self.myGibsonDialog.ui.tableWidget.setItem(index, 4, unit5)

                cell_Text = QTextEdit()
                cell_Text.setPlainText(VDJseq)
                #cell_Text.resize(cell_Text.size().width(), 20)
                cell_Text.rowindex = SeqName
                cell_Text.textChanged.connect(self.myGibsonDialog.updateData)
                self.myGibsonDialog.ui.tableWidget.setCellWidget(index, 4, cell_Text)

                cell_TextAA = QTextEdit()
                cell_TextAA.setPlainText(AAseq)
                cell_TextAA.setReadOnly(True)
                #cell_TextAA.rowindex = index
                self.myGibsonDialog.ui.tableWidget.setCellWidget(index, 5, cell_TextAA)

                if checkRes == "Good":
                    self.myGibsonDialog.ui.tableWidget.item(index, 0).setBackground(Qt.green)
                elif checkRes == "Check ORF":
                    self.myGibsonDialog.ui.tableWidget.item(index, 0).setBackground(Qt.yellow)
                elif checkRes == "Jend Mut":
                    self.myGibsonDialog.ui.tableWidget.item(index, 0).setBackground(Qt.yellow)
                else:
                    self.myGibsonDialog.ui.tableWidget.item(index, 0).setBackground(Qt.red)
                index += 1

        # disable edit
        self.myGibsonDialog.ui.tableWidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        # resize table
        self.myGibsonDialog.ui.tableWidget.resizeColumnsToContents()
        self.myGibsonDialog.ui.tableWidget.resizeRowsToContents()
        for index in range(num_row):
            self.myGibsonDialog.ui.tableWidget.setRowHeight(index, 90)
        self.myGibsonDialog.ui.tableWidget.setColumnWidth(4, 600)
        # show sort indicator
        self.myGibsonDialog.ui.tableWidget.horizontalHeader().setSortIndicatorShown(True)
        # connect sort indicator to slot function
        self.myGibsonDialog.ui.tableWidget.horizontalHeader().sectionClicked.connect(self.myGibsonDialog.sort)
        # set signal
        #self.myGibsonDialog.ui.tableWidget.cellChanged.connect(self.myGibsonDialog.updateData)
        self.myGibsonDialog.ui.tableWidget.currentCellChanged.connect(self.myGibsonDialog.updateSelection)
        self.myGibsonDialog.GibsonUpdateSelectionSignal.connect(self.select_tree_by_name)
        self.myGibsonDialog.LogFileSignal.connect(self.displayLog)
        # show dialog
        self.myGibsonDialog.show()
    elif option == 'CSV format customized fields':
        if self.ui.tabWidget.currentIndex() == 11:
            if len(self.AntibodyCandidates) == 0:
                selected_list = self.CheckedRecords
            else:
                selected_list = self.AntibodyCandidates
        else:
            selected_list = self.CheckedRecords

        WHEREStatement = ' WHERE SeqName IN ("' + '","'.join(selected_list) + '")'
        if len(selected_list) == 0:
            question = 'You did not select any records, export all?'
            buttons = 'YN'
            answer = questionMessage(self, question, buttons)
            if answer == 'Yes':
                WHEREStatement = ' WHERE 1'
            else:
                return

        self.myExportOptionDialog = ExportOptionDialog()
        self.myExportOptionDialog.WHEREStatement = WHEREStatement
        self.myExportOptionDialog.DBFilename = DBFilename

        if DBFilename != '' and DBFilename != None and DBFilename != 'none':
            SQLStatement = 'SELECT display,Field,FieldNickName,FieldType,FieldComment,ID FROM fieldsname ORDER BY ID'
            DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
            header_list = ['Selected', 'Field', 'Field nickname', 'Field type', 'Field comment','ID']
            num_row = len(DataIn)
            num_col = len(header_list)
            self.myExportOptionDialog.ui.tableWidget.setRowCount(num_row)
            self.myExportOptionDialog.ui.tableWidget.setColumnCount(num_col)
            self.myExportOptionDialog.ui.tableWidget.setHorizontalHeaderLabels(header_list)

            #self.myExportOptionDialog.ui.tableWidget.setDragDropMode(QAbstractItemView.InternalMove)
            self.myExportOptionDialog.ui.tableWidget.setSelectionMode(QAbstractItemView.SingleSelection)
            self.myExportOptionDialog.ui.tableWidget.setSelectionBehavior(QAbstractItemView.SelectRows)

            for row_index in range(num_row):
                # col 0
                cell_checkBox = QCheckBox()
                cell_checkBox.setChecked(False)
                self.myExportOptionDialog.ui.tableWidget.setCellWidget(row_index, 0, cell_checkBox)

                # col 2:
                for col_index in range(1, num_col):
                    unit = QTableWidgetItem(str(DataIn[row_index][col_index]))
                    self.myExportOptionDialog.ui.tableWidget.setItem(row_index, col_index, unit)

            # disable edit
            self.myExportOptionDialog.ui.tableWidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
            # re-size column size
            self.myExportOptionDialog.ui.tableWidget.resizeColumnsToContents()

        # show dialog
        self.myExportOptionDialog.show()
    elif option == 'Antibody Patent report':
        if len(self.AntibodyCandidates) == 0:
            Msg = 'Nothing in Antibody candidate list!'
            self.ShowMessageBox([0, Msg])
        else:
            # check if all HC/LCs are correctly paired
            if self.ui.tableWidgetHC.rowCount() != self.ui.tableWidgetLC.rowCount():
                Msg = 'Number of HCs and LCs are not match!'
                QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
                return
            # fill HC data
            HC_dict = {}
            for index in range(self.ui.tableWidgetHC.rowCount()):
                barcode = self.ui.tableWidgetHC.item(index, 8).text()
                name = self.ui.tableWidgetHC.item(index, 0).text()
                if barcode in HC_dict.keys():
                    Msg = 'Multiple HCs shared same barcode!\n' + HC_dict[barcode] + ' and ' + name
                    QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
                    return
                HC_dict[barcode] = name
            # fill LC data
            LC_dict = {}
            for index in range(self.ui.tableWidgetLC.rowCount()):
                barcode = self.ui.tableWidgetLC.item(index, 7).text()
                name = self.ui.tableWidgetLC.item(index, 0).text()
                if barcode in LC_dict.keys():
                    Msg = 'Multiple LCs shared same barcode!\n' + LC_dict[barcode] + ' and ' + name
                    QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
                    return
                LC_dict[barcode] = name
            # make data
            Paired_data = []
            for barcode in HC_dict.keys():
                if barcode in LC_dict.keys():
                    Paired_data.append((barcode, HC_dict[barcode], LC_dict[barcode]))
                else:
                    Msg = 'A barcode can not find LC!\n' + barcode
                    QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
                    return
            # open dialog
            self.myPatentDialog = PatentDialog()
            self.myPatentDialog.DBFilename = DBFilename

            horizontalHeader = ['Antibody Name', 'Barcode', 'Heavy Chain', 'Light Chain']
            num_row = len(Paired_data)
            num_col = len(horizontalHeader)
            self.myPatentDialog.ui.tableWidget.setRowCount(num_row)
            self.myPatentDialog.ui.tableWidget.setColumnCount(num_col)
            self.myPatentDialog.ui.tableWidget.setHorizontalHeaderLabels(horizontalHeader)
            self.myPatentDialog.ui.tableWidget.horizontalHeader().setStretchLastSection(True)
            self.myPatentDialog.ui.tableWidget.setSelectionMode(QAbstractItemView.SingleSelection)
            self.myPatentDialog.ui.tableWidget.setSelectionBehavior(QAbstractItemView.SelectRows)
            index = 0
            for record in Paired_data:
                cell_Text = QLineEdit()
                cell_Text.setMinimumSize(250,20)
                cell_Text.setText(record[1])
                self.myPatentDialog.ui.tableWidget.setCellWidget(index, 0, cell_Text)

                unit1 = QTableWidgetItem(record[0])
                unit2 = QTableWidgetItem(record[1])
                unit3 = QTableWidgetItem(record[2])
                self.myPatentDialog.ui.tableWidget.setItem(index, 1, unit1)
                self.myPatentDialog.ui.tableWidget.setItem(index, 2, unit2)
                self.myPatentDialog.ui.tableWidget.setItem(index, 3, unit3)

                index += 1
            # disable edit
            self.myPatentDialog.ui.tableWidget.setEditTriggers(QAbstractItemView.NoEditTriggers)
            # resize table
            self.myPatentDialog.ui.tableWidget.resizeColumnsToContents()
            self.myPatentDialog.ui.tableWidget.resizeRowsToContents()
            self.myPatentDialog.show()
    print('done')

def checkJend(GeneType, JendSeq):
    if GeneType == 'Heavy':
        if JendSeq.upper() == 'TCCTCA':
            res = 'Good'
        else:
            if VGenesSeq.Translator(JendSeq.upper(), 0)[0] == 'SS':
                res = 'Jend Mut'
            else:
                res = 'Jend Error'
    elif GeneType == 'Kappa':
        if JendSeq.upper() in ['TCGAAC', 'ATTAAA', 'ATCAAA']:
            res = 'Good'
        else:
            if VGenesSeq.Translator(JendSeq.upper(), 0)[0] in ['SN', 'IK']:
                res = 'Jend Mut'
            else:
                res = 'Jend Error'
    elif GeneType == 'Lambda':
        if JendSeq.upper() == 'GTCCTA':
            res = 'Good'
        else:
            if VGenesSeq.Translator(JendSeq.upper(), 0)[0] == 'VL':
                res = 'Jend Mut'
            else:
                res = 'Jend Error'

    return res

def getORF(alignment):
    # get ORF info from alignment
    lines = alignment.split('\n')
    line_num = 0
    for line in lines:
        match = re.match(r'^(\s+)<-+FR1', line)
        if match:
            num1 = len(match.group(1))
            aa_line = lines[line_num + 1]
            match1 = re.match(r'^\s+', aa_line)
            num2 = len(match1.group())
            ORF = num2 - num1 - 1
            if ORF < 0:
                ORF = 0
            return ORF
        line_num += 1