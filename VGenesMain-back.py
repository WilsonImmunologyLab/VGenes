__author__ = 'wilsonp'
import os
import time

from PyQt5.QtCore import pyqtSlot, QTimer, QDateTime, Qt, QSortFilterProxyModel, QModelIndex
from PyQt5 import QtWidgets
import VReports
global OldName
from ui_VGenesMain import Ui_MainWindow
import IgBLASTer
from VgenesTextEdit import VGenesTextMain
import VGenesSQL
import VGenesSeq
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtPrintSupport import QPrintDialog, QPrinter
import VGenesCloneCaller
from ui_Import_Dialogue import Ui_DialogImport
global GLMsg
GLMsg = True
import csv
from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
	setText
from ui_VGenesStartUpDialogue import Ui_VGenesStartUpDialog

from ui_VGenesTextEdit import ui_TextEditor
from VGenesProgressBar import ui_ProgressBar
# from VGenesPYQTSqL import EditableSqlModel, initializeModel , createConnection

from PyQt5.QtCore import Qt, QObject, QEvent
from PyQt5.QtGui import QTextCursor, QFont, QPixmap, QTextCharFormat, QBrush, QColor, QTextCursor, QCursor
from PyQt5.QtWidgets import QApplication, QTableView
from PyQt5.QtSql import QSqlQuery, QSqlQueryModel
from operator import itemgetter
import itertools

from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtSql import QSqlDatabase, QSqlQuery

global LastPushed
LastPushed = ''

global data
global LastSelected
global DBFilename
global wasClicked
global LineLen
LastSelected = ()
global DontFindTwice
DontFindTwice = False
wasClicked = False
data = []
global JustMovedIt
JustMovedIt = True
global RefreshSQL
RefreshSQL = 'select * from vgenesdb ORDER BY Project, Grouping, SubGroup, SeqName'
# global data
global FieldList
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
global RealNameList
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
                "Germlne CDR3 begin", "Germline CDR3 end", "Blank6", "Blank7", "Blank8", "Blank9", "Blank10", "Blank11",
                "Blank12", "Blank13", "Blank14", "Blank15", "Blank16", "Blank17", "Blank18", "Blank19", "Blank20", "ID"]
global NameIndex
NameIndex = {}
global FieldsChanged
FieldsChanged = []
# global ChangesBuffer
# ChangesBuffer = []
global FieldChanged
FieldChanged = False
global FirstupdateF
FirstupdateF = True


class VGenesTextMain(QtWidgets.QMainWindow, ui_TextEditor):
	def __init__(self, parent=None):
		QtWidgets.QMainWindow.__init__(self, parent)
		# super(VGenesTextMain, self).__init__()
		self.setupUi()


# class VGenesProgressBar(QtWidgets.QMainWindow, ui_ProgressBar):
#     def __init__(self, parent=None):
#         QtWidgets.QMainWindow.__init__(self, parent)
#
#         self.setupUi()

class StartUpDialogue(QtWidgets.QDialog, Ui_VGenesStartUpDialog):
	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		self.setupUi(self)
		# self.TextEdit = VGenesTextMain()
		self.PopulateCombo()

		self.cboRecent.currentTextChanged.connect(self.cboRecentTextChanged)

	def PopulateCombo(self):
		# # todo need to make this filename fall in VGenes directory upon deployment
		try:
			filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'RecentPaths.vtx')

			with open(filename, 'r') as currentfile:
				self.cboRecent.clear()
				self.cboRecent.addItem('Select a recently opened VGenes database file')
				for line in currentfile:
					line = line.replace('\n', '')
					if os.path.isfile(line):
						self.cboRecent.addItem(line)
					else:
						Vgenes.UpdateRecentList(line, False)

		except:
			filename = '/Volumes/Promise Pegasus/Dropbox/VGenes/RecentPaths.vtx'

			with open(filename, 'r') as currentfile:
				self.cboRecent.clear()
				self.cboRecent.addItem('Select a recently opened VGenes database file')
				for line in currentfile:
					line = line.replace('\n', '')
					if os.path.isfile(line):
						self.cboRecent.addItem(line)
					else:
						Vgenes.UpdateRecentList(line, False)

	@pyqtSlot()
	def on_btnNew_clicked(self):

		global StartUpAnswer
		StartUpAnswer = 'New'
		Vgenes.StartUpClicked()

		self.close()

	@pyqtSlot()
	def on_cboOpen_clicked(self):

		global StartUpAnswer
		StartUpAnswer = 'Open'
		Vgenes.StartUpClicked()

		self.close()

	def cboRecentTextChanged(self):
		# self.cboRecent.currentTextChanged()
		# self.close()
		global StartUpAnswer

		StartUpAnswer = self.cboRecent.currentText()
		if os.path.isfile(self.cboRecent.currentText()):
			self.close()
			StartUpAnswer = 'Recent' + StartUpAnswer
			Vgenes.StartUpClicked()
			# return Answer

		else:

			Query = 'This file ' + self.cboRecent.currentText() + ' no longer exists at this location.'

			QtWidgets.QMessageBox.critical(None, "File error", Query, QtWidgets.QMessageBox.Cancel)
			StartUpAnswer += '\n'
			Vgenes.UpdateRecentList(StartUpAnswer, False)
			self.close()
			Vgenes.ApplicationStarted()


class ImportDialogue(QtWidgets.QDialog, Ui_DialogImport):
	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		self.setupUi(self)
		self.TextEdit = VGenesTextMain()

		global answer3
		answer3 = 'No'

	@pyqtSlot()
	def on_rdoChoose_clicked(self):

		if self.rdoChoose.isChecked():
			self.comboBoxProject.setEditable(True)
			self.comboBoxGroup.setEditable(True)
			self.comboBoxSubgroup.setEditable(True)
			self.comboBoxProject.setCurrentText('')
			self.comboBoxGroup.setCurrentText('')
			self.comboBoxSubgroup.setCurrentText('')

			fields = ['Project']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = VGenesSQL.MakeSQLStatement(self, fields, data[0])

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			if len(DataIn) > 0:
				for item in DataIn:
					self.comboBoxProject.addItem(item[0])
			DataIn.clear()

			fields = ['Grouping']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = VGenesSQL.MakeSQLStatement(self, fields, data[0])

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if len(DataIn) > 0:
				for item in DataIn:
					self.comboBoxGroup.addItem(item[0])
			DataIn.clear()

			fields = ['SubGroup']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = VGenesSQL.MakeSQLStatement(self, fields, data[0])

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			if len(DataIn) > 0:
				for item in DataIn:
					self.comboBoxSubgroup.addItem(item[0])


		else:
			self.comboBoxProject.setEditable(False)
			self.comboBoxGroup.setEditable(False)
			self.comboBoxSubgroup.setEditable(False)

	@pyqtSlot()
	def on_checkBoxFileStruc_clicked(self):

		if self.checkBoxFileStruc.isChecked():
			self.comboBoxProject.setEditable(False)
			self.comboBoxGroup.setEditable(False)
			self.comboBoxSubgroup.setEditable(False)
			self.comboBoxProject.setCurrentText('')
			self.comboBoxGroup.setCurrentText('')
			self.comboBoxSubgroup.setCurrentText('')

		else:
			self.comboBoxProject.setEditable(True)
			self.comboBoxGroup.setEditable(True)
			self.comboBoxSubgroup.setEditable(True)

	# Need code for function to populate identical to these on main form

	@pyqtSlot()
	def on_buttonBox_accepted(self):
		# Alldone = False
		Alldone = self.InitiateImport('none')

		if Alldone == True:
			self.close()

	@pyqtSlot()
	def on_buttonBox_rejected(self):

		self.close()

	def timingcheck(self):
		# TODO can delete this if gets work

		# import IgBLASTer

		# stmtS = "/Users/PCW-MacBookProRet/Dropbox/VGenes/Database/SFV-005H.nt"
		DirTo = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes')
		import timeit, os
		# todo change to app folder
		try:
			os.chdir(DirTo)
		except:
			os.chdir('/Volumes/Promise Pegasus/Dropbox/VGenes/VGenes')

		t = timeit.timeit('self.findTableViewRecord(self, FieldName)")', 'import IgBLASTer', number=10000)
		print(t)

	def InitiateImport(self, Filenamed):
		# need to transfer species grouping to IgBlaster
		answer = ''
		thetype = ''
		species = ''
		datalist = []
		global answer3
		# answerTo = answer3

		if self.radioButHuman.isChecked():
			species = 'Human'
		elif self.radioButMouse.isChecked():
			species = 'Mouse'

		if self.radioButFASTA.isChecked():
			thetype = 'FASTA'
		elif self.radioButIndSeq.isChecked():
			thetype = 'Sequence'

		if thetype == 'FASTA':
			if Filenamed == 'none':
				pathname = openFile(self, thetype)
			else:
				pathname = Filenamed[0]
			if pathname == None:
				return False

		elif thetype == 'Sequence':
			filenames = openFiles(self)
			pathname = self.ProcessSeqFiles(filenames)
			# filename = filenames[0]

		if self.rdoProductive.isChecked() == True:
			GetProductive = True
		else:
			GetProductive = False

		if pathname == None:
			return

		if self.checkBoxFileStruc.isChecked():

			(dirname, filename) = os.path.split(pathname)
			dirparts = dirname.split('/')
			NumParts = len(dirparts)

			subgroup = dirparts[NumParts - 1]

			for i in range((self.comboBoxSubgroup.count()) - 1):
				if self.comboBoxSubgroup.itemText(i) == subgroup:
					self.comboBoxSubgroup.setCurrentText(subgroup)
			if self.comboBoxSubgroup.currentText() != subgroup:
				self.comboBoxSubgroup.addItem(subgroup)
				self.comboBoxSubgroup.setCurrentText(subgroup)
			# need to check this code and add to other combos

			if NumParts > 2:
				grouping = dirparts[NumParts - 2]
				for i in range((self.comboBoxGroup.count()) - 1):
					if self.comboBoxGroup.itemText(i) == grouping:
						self.comboBoxGroup.setCurrentText(grouping)
				if self.comboBoxGroup.currentText() != grouping:
					self.comboBoxGroup.addItem(grouping)
					self.comboBoxGroup.setCurrentText(grouping)
			# need to check this code and add to other combos

			else:
				grouping = dirparts[NumParts - 1]
				for i in range((self.comboBoxGroup.count()) - 1):
					if self.comboBoxGroup.itemText(i) == grouping:
						self.comboBoxGroup.setCurrentText(grouping)
				if self.comboBoxGroup.currentText() != grouping:
					self.comboBoxGroup.addItem(grouping)
					self.comboBoxGroup.setCurrentText(grouping)

			if NumParts > 3:
				project = dirparts[NumParts - 3]
				for i in range((self.comboBoxProject.count()) - 1):
					if self.comboBoxProject.itemText(i) == project:
						self.comboBoxProject.setCurrentText(project)
				if self.comboBoxProject.currentText() != project:
					self.comboBoxProject.addItem(project)
					self.comboBoxProject.setCurrentText(project)

			elif NumParts > 2:
				project = dirparts[NumParts - 1]
				for i in range((self.comboBoxProject.count()) - 1):
					if self.comboBoxProject.itemText(i) == project:
						self.comboBoxProject.setCurrentText(project)
				if self.comboBoxProject.currentText() != project:
					self.comboBoxProject.addItem(project)
					self.comboBoxProject.setCurrentText(project)

			else:
				project = dirparts[NumParts - 1]
				self.comboBoxProject.setCurrentText(project)

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)

			# print(project + '\n' + grouping + '\n' + subgroup)

			# TODO make this query fire only if project or grouping is odd
			# TODO or alternative make open dialogue come up first...but that
			# type = 'YN'
			# Query = "Verify project and groupings based on file structure: continue?"
			# reply = questionMessage(self, Query,type)
			# if reply == 'Yes':
			#     # self.ShowVGenesText(ErLog)
			# ErReport, IgBLASTAnalysis = IgBLASTer.IgBLASTit(pathname, datalist)
			# BlastIsDone = False

			IgBLASTAnalysis = IgBLASTer.IgBLASTit(pathname, datalist)

			# while IgBLASTAnalysis is None:
			#     NotDone = True
			#


			# IgBlaster will load values into current database
			# else:
			#     return False
			#

		elif self.rdoChoose.isChecked():
			checklabel = {}
			# project = ''
			# grouping = ''
			# subgroup = ''
			if Filenamed == 'none':
				project = self.comboBoxProject.currentText()
				grouping = self.comboBoxGroup.currentText()
				subgroup = self.comboBoxSubgroup.currentText()
			else:
				project = Filenamed[1]
				grouping = Filenamed[2]
				subgroup = Filenamed[3]

			if project == '': project = 'none'
			if grouping == '': grouping = 'none'
			if subgroup == '': subgroup = 'none'

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)

			buildquery = 'No label was chosen for the following:'
			GiveWarning = False

			# BlastIsDone = False

			IgBLASTAnalysis = IgBLASTer.IgBLASTit(pathname, datalist)
			#
			# while IgBLASTAnalysis is None:
			#     NotDone = True
		elif self.rdoFunction.isChecked():
			project = 'ByFunction'
			grouping = ''
			subgroup = ''
			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)

			IgBLASTAnalysis = IgBLASTer.IgBLASTit(pathname, datalist)

		# print(pathname)
		Startprocessed = len(IgBLASTAnalysis)
		if Startprocessed == 0:
			self.close()

		Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)
		answer3 = answer

		# todo need code to verify a database is opne before you can import sequences.

		Vgenes.LoadDB(DBFilename)

		ErlogFile = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'IgBlast', 'database',
		                         'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
		i = 0
		if Startprocessed != Processed:
			newErLog = str(Processed) + ' sequences were analyzed by IgBLAST \n'
			with open(ErlogFile,
			          'r') as currentFile:  # using with for this automatically closes the file even if you crash
				for line in currentFile:
					if i > 0:
						newErLog += line
					i += 1

			with open(ErlogFile, 'w') as currentFile:
				currentFile.write(newErLog)

		self.close()
		self.ShowVGenesText(ErlogFile)
		# self.EditableSqlModel.refresh()
		# else:
		#     type = 'OK'
		#     Query = "Analysis failed: There were no valid Ig sequences input"
		#     reply = questionMessage(self, Query, type)

	@pyqtSlot()
	def on_btnImportOldVGenes_clicked(self):
		from operator import itemgetter
		msg = 'This function imports a comma separated values (CSV) file formatted as: Project, Group, Subgroup, Name, Sequence'
		buttons = 'OKC'
		global answer3
		answer3 = 'No'
		answer = informationMessage(self, msg, buttons)
		if answer == 'Cancel':
			print('no file')
			return
		Pathname = openFile(self, 'CSV')
		self.checkBoxFileStruc.setChecked(False)
		self.rdoChoose.setChecked(True)
		self.radioButFASTA.setChecked(True)
		SeqList = []
		fieldis = ''
		with open(Pathname, 'r') as currentfile:
			for line in currentfile:
				fields = line.split(',')
				if len(fields) == 5:
					for i in range(0, 4):
						fieldis = fields[i]
						fieldis = fieldis.upper()
						# if fieldis == '01B_STAN_17-006-1A03H_VH3':
						#     print('stop')
						fields[i] = fieldis.replace('"', '')

						if fields[i] == '' and i < 3:
							fields[i] = 'none'
					SeqList.append(fields)
		SeqList.sort(key=itemgetter(0, 1, 2, 3))
		LastP = ''
		LastG = ''
		LastSG = ''
		StartNewFASTA = True
		FASTAFile = ''
		NewFile = ''
		FirstOne = True
		DataPass = []
		FileNamed = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'IgBlast', 'database',
		                         'WorkingFile.nt')  # '/Applications/IgBlast/database/WorkingFile.nt'
		DataPass.append(FileNamed)
		DataPass.append('none')
		DataPass.append('none')
		DataPass.append('none')
		for seq in SeqList:
			project = seq[0]
			group = seq[1]
			subgroup = seq[2]
			SeqName = '>' + seq[3] + '\n'
			Sequence = seq[4] + '\n'

			if project != LastP:
				LastP = project
				StartNewFASTA = True
			# else:
			#     StartNewFASTA = False

			if group != LastG:
				LastG = group
				StartNewFASTA = True
			# else:
			#     StartNewFASTA = False

			if subgroup != LastSG:
				LastSG = subgroup
				StartNewFASTA = True
			# else:
			#     StartNewFASTA = False


			if StartNewFASTA == True:  # then write it and clear it clear
				StartNewFASTA = False
				if FirstOne == False:  # firstone is empty
					with open(FileNamed, 'w') as currentFile:
						currentFile.write(NewFile)

					Alldone = self.InitiateImport(DataPass)
					DataPass[1] = project
					if project == 'STAN-004':
						print('stop')
					DataPass[2] = group
					DataPass[3] = subgroup
				else:
					DataPass[1] = project
					DataPass[2] = group
					DataPass[3] = subgroup

				FirstOne = False
				NewFile = ''

			NewFile += SeqName
			NewFile += Sequence

		if NewFile != '':
			with open(FileNamed, 'w') as currentFile:
				currentFile.write(NewFile)
			Alldone = self.InitiateImport(DataPass)

		if Alldone == True:
			self.close()

	def ProcessSeqFiles(self, fileNames):
		if fileNames:
			FASTAfile = []
			FASTAparts = []
			for files in fileNames:
				WriteFASTA = True
				(dirname, filename) = os.path.split(files)  # parses filename from path
				(shortname, extension) = os.path.splitext(filename)  # parses filename into name and extension
				NameLine = '>' + shortname + '\n'

				# print(FASTAfile)
				os.chdir(dirname)

				with open(files,
				          'r') as currentFile:  # using with for this automatically closes the file even if you crash
					readLine = str(currentFile.read())

				CountNucs = readLine.count('a') + readLine.count('A') + readLine.count('g') + readLine.count(
					'G') + readLine.count('c') + readLine.count('C') + readLine.count('t') + readLine.count('T')
				PercentNuc = CountNucs / len(readLine)

				if readLine[0] == '>':
					query = files + ' appears to ba a FASTA file, process as such?'
					answer = questionMessage(self, query, 'YN')
					if answer == 'Yes':
						FASTAparts.append(readLine)
						WriteFASTA = False
					else:
						WriteFASTA = False

				# code to see if mostly a good sequence


				elif len(readLine) < 30:
					query = files + ' is a short sequence (<30 nucleotides), analyze anyways?'
					answer = questionMessage(self, query, 'YN')
					if answer == 'Yes':
						WriteFASTA = True
					else:
						WriteFASTA = False


				elif PercentNuc < 0.8:
					query = files + ' has over 20% of characters that are not nucleotides (A, G, C, or T), analyze anyways?'
					answer = questionMessage(self, query, 'YN')
					if answer == 'Yes':
						WriteFASTA = True
					else:
						WriteFASTA = False

				# print(readLine)
				if WriteFASTA == True:
					FASTAfile.append(NameLine)

					readLine = readLine.replace('\n', '').replace('\r', '')
					readLine += '\n'
					FASTAfile.append(readLine)

					# print(FASTAfile)
			FinalFASTA = ''.join(FASTAfile)

			now = 'FASTA from ' + time.strftime('%c') + '.nt'
			FASTAFileName = os.path.join(dirname, now)
			# need to test


			with open(FASTAFileName,
			          'w') as currentFile:  # using with for this automatically closes the file even if you crash
				currentFile.write(FinalFASTA)

			return FASTAFileName
			# ErLog = IgBLASTer.ProcessFASTA(filename)
			# if ErLog != '':
			#     type = 'YN'
			#     Query = "There were some bad sequences, would you like to see the Error log?"
			#
			#     reply = questionMessage(Query,type)
			#     if reply == 'Yes':
			#         self.ShowVGenesText(ErLog)

	def ShowVGenesText(self, filename):

		self.TextEdit.show()
		if filename != '':
			self.TextEdit.loadFile(filename)


class VGenesForm(QtWidgets.QMainWindow):
	def __init__(self):  # , parent=None):
		super(VGenesForm, self).__init__()  # parent)

		self.ui = Ui_MainWindow()

		self.ui.setupUi(self)

		self.ui.comboBoxSpecies.currentTextChanged.connect(self.on_comboBoxSpecies_editTextChanged)
		self.ui.cboTreeOp1.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboTreeOp2.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboTreeOp2.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboDecorate.currentTextChanged.connect(self.DecoratePeptide)
		self.ui.treeWidget.itemChanged.connect(self.handleChanged)
		self.ui.treeWidget.itemSelectionChanged.connect(self.TreeSelectChanged)
		self.ui.treeWidget.clicked['QModelIndex'].connect(self.treeWidgetClicked)
		self.ui.treeWidget.doubleClicked['QModelIndex'].connect(self.CheckMultiple)
		self.ui.cboReportOptions.currentTextChanged.connect(self.ReportOptions)

		self.ui.lcdNumber_max.display(self.ui.horizontalScrollBar.maximum())
		self.ui.dial.setMaximum(self.ui.horizontalScrollBar.maximum())

		self.TextEdit = VGenesTextMain()
		# self.VGProgress = VGenesProgressBar()
		self.ImportOptions = ImportDialogue()

	# def ShowProgressBar(self):
	#     self.VGProgress.handleButton()


	def ShowVGenesTextEdit(self, textToShow, style):

		if style == 'aligned':
			FontIs = self.TextEdit.textEdit.currentFont()
			font = QFont(FontIs)

			# FontSize = int(font.pointSize())
			font.setPointSize(12)
			font.setFamily('Courier New')

			self.TextEdit.textEdit.setFont(font)

		elif style == 'standard':
			FontIs = self.TextEdit.textEdit.currentFont()
			font = QFont(FontIs)

			# FontSize = int(font.pointSize())
			font.setPointSize(12)
			font.setFamily('Lucida Grande')

			self.TextEdit.textEdit.setFont(font)

		elif style == 'ProteinReport':
			FontIs = self.TextEdit.textEdit.currentFont()
			font = QFont(FontIs)

			# FontSize = int(font.pointSize())
			font.setPointSize(6)
			font.setFamily('Courier New')

			self.TextEdit.textEdit.setFont(font)

		self.TextEdit.show()

		self.TextEdit.textEdit.setText(textToShow)


	def initializeTreeView(self, SQLFields):

		self.ui.treeWidget.clear()

		self.TreeAddItems(self.ui.treeWidget.invisibleRootItem(), SQLFields)
		global wasClicked
		wasClicked = True

	def TreeAddItems(self, parent, SQLFields):
		dirnamed, filenamed = os.path.split(DBFilename)

		FieldNum = 0
		Field1 = ''
		Field2 = ''
		Field3 = ''
		i1 = 0
		i2 = 0
		i3 = 0
		i4 = 0
		for sfield in SQLFields:
			if sfield != 'None' and sfield is not None:
				if FieldNum == 0: Field1 = sfield
				if FieldNum == 1: Field2 = sfield
				if FieldNum == 2: Field3 = sfield
				FieldNum += 1
			else:
				if FieldNum == 0: Field1 = 'None'
				if FieldNum == 1: Field2 = 'None'
				if FieldNum == 2: Field3 = 'None'
				FieldNum += 1

		# todo make header actaul values not field names

		HeaderLabel = Field1 + ' / ' + Field2 + ' / ' + Field3 + ':'
		self.ui.treeWidget.setHeaderLabel(HeaderLabel)

		column = 0
		VTree_item = self.addParent(parent, column, filenamed, 'Top parent')

		# Field 1 selection
		if Field1 != 'None':  # some organizing based on item 1
			SQLStatement = 'SELECT DISTINCT ' + Field1 + ' FROM vgenesDB ORDER BY ' + Field1
			FieldValues = VGenesSQL.readData(DBFilename, SQLStatement)
			for Vcolumn in FieldValues:
				self.addChild(VTree_item, column, Vcolumn, 'Field-1')
				# i += 1
		else:  # All classification fields marked 'None'
			SQLStatement = 'SELECT SeqName FROM vgenesDB ORDER BY Project, Grouping, SubGroup, SeqName'
			FieldValues = VGenesSQL.readData(DBFilename, SQLStatement)
			for Vcolumn in FieldValues:
				self.addChild(VTree_item, column, Vcolumn, 'Field-1')
			self.addChild(VTree_item, column, ' ', ' ')
			return

		# Field 2
		I1 = 0
		if Field2 != 'None':  # There is selection criteria in field 2, iterates through for all Field1s
			for Vcolumn in FieldValues:
				SQLStatement2 = 'SELECT DISTINCT ' + Field2 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" ORDER BY Project, Grouping, SubGroup, SeqName'
				FieldValues2 = VGenesSQL.readData(DBFilename, SQLStatement2)
				for Vcolumn2 in FieldValues2:
					self.addChild(VTree_item.child(i1), column, Vcolumn2, 'Field-2')
				i1 += 1

		else:
			for Vcolumn in FieldValues:
				SQLStatement2 = 'SELECT SeqName FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" ORDER BY Project, Grouping, SubGroup, SeqName'
				FieldValues2 = VGenesSQL.readData(DBFilename, SQLStatement2)
				for Vcolumn2 in FieldValues2:
					self.addChild(VTree_item.child(i1), column, Vcolumn2, 'Field-2')
				i1 += 1
				# self.addChild( VTree_item, column, ' ', ' ')
				return

		# Field 3
		i1 = 0
		i2 = 0
		if Field3 != 'None':  # There is selection criteris in field 2, iterates through for all Field1s
			for Vcolumn in FieldValues:
				SQLStatement2 = 'SELECT DISTINCT ' + Field2 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" ORDER BY Project, Grouping, SubGroup, SeqName'
				FieldValues2 = VGenesSQL.readData(DBFilename, SQLStatement2)
				i2 = 0
				for Vcolumn2 in FieldValues2:
					SQLStatement3 = 'SELECT DISTINCT ' + Field3 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" AND ' + Field2 + ' = "' + Vcolumn2 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
					FieldValues3 = VGenesSQL.readData(DBFilename, SQLStatement3)
					for Vcolumn3 in FieldValues3:
						self.addChild(VTree_item.child(i1).child(i2), column, Vcolumn3, 'Field-3')
					i2 += 1
				i1 += 1

		else:
			for Vcolumn in FieldValues:
				SQLStatement2 = 'SELECT DISTINCT ' + Field2 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" ORDER BY Project, Grouping, SubGroup, SeqName'
				FieldValues2 = VGenesSQL.readData(DBFilename, SQLStatement2)
				i2 = 0
				for Vcolumn2 in FieldValues2:
					SQLStatement3 = 'SELECT SeqName FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" AND ' + Field2 + ' = "' + Vcolumn2 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
					FieldValues3 = VGenesSQL.readData(DBFilename, SQLStatement3)
					for Vcolumn3 in FieldValues3:
						self.addChild(VTree_item.child(i1).child(i2), column, Vcolumn3, 'Field-3')
					i2 += 1
				i1 += 1
				return

		# Last set only done if all 3 others selected
		i1 = 0
		i2 = 0
		i3 = 0
		for Vcolumn in FieldValues:
			SQLStatement2 = 'SELECT DISTINCT ' + Field2 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" ORDER BY Project, Grouping, SubGroup, SeqName'
			FieldValues2 = VGenesSQL.readData(DBFilename, SQLStatement2)
			i2 = 0
			for Vcolumn2 in FieldValues2:
				SQLStatement3 = 'SELECT DISTINCT ' + Field3 + ' FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" AND ' + Field2 + ' = "' + Vcolumn2 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
				FieldValues3 = VGenesSQL.readData(DBFilename, SQLStatement3)
				i3 = 0
				for Vcolumn3 in FieldValues3:
					SQLStatement4 = 'SELECT SeqName FROM vgenesDB WHERE ' + Field1 + ' = "' + Vcolumn + '" AND ' + Field2 + ' = "' + Vcolumn2 + '" AND ' + Field3 + ' = "' + Vcolumn3 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
					FieldValues4 = VGenesSQL.readData(DBFilename, SQLStatement4)
					for Vcolumn4 in FieldValues4:
						self.addChild(VTree_item.child(i1).child(i2).child(i3), column, Vcolumn4, 'Field-4')
						# self.addChild(VTree_item.child(i1).child(i2).child(i3), 1, 'test', 'Field-4')
					# self.addChild(VTree_item.child(i1).child(i2).child(i3), column, 'Blank', 'Field-4')
					i3 += 1
				i2 += 1
			i1 += 1

	def addParent(self, parent, column, title, data):
		item = QtWidgets.QTreeWidgetItem(parent, [title])
		item.setData(column, Qt.UserRole, data)
		item.setChildIndicatorPolicy(QtWidgets.QTreeWidgetItem.ShowIndicator)
		item.setExpanded(True)
		return item

	def addChild(self, parent, column, title, data):
		item = QtWidgets.QTreeWidgetItem(parent, [title])
		item.setData(column, Qt.UserRole, data)
		item.setCheckState(column, Qt.Unchecked)
		return item

	@pyqtSlot()
	def treeWidgetClicked(self):
		global wasClicked
		wasClicked = True

	def handleChanged(self, item, column):

		global wasClicked
		if wasClicked == False:
			return

		wasClicked = False

		Selected = self.ui.treeWidget.selectedItems()
		NumSelected = len(Selected)
		if len(Selected) > 1:
			self.CheckMultiple()
			if item.checkState(column) == Qt.Checked:
				item.setCheckState(0, Qt.Unchecked)
			elif item.checkState(column) == Qt.Unchecked:
				item.setCheckState(0, Qt.Checked)
			return

		if item.checkState(column) == Qt.Checked:
			self.CheckBelow(item, True)
			self.CheckUp(item, True)
			# childs = item.childCount()

		if item.checkState(column) == Qt.Unchecked:
			self.CheckBelow(item, False)
			self.CheckUp(item, False)
			print("unchecked"), item, item.text(column)

	def getTreePathUp(self, item):
		path = []
		while item is not None:
			item = item.parent()
			if item is not None:
				path.append(str(item.text(0)))

		return path  # '/'.join(reversed(path))

	def CheckBelow(self, item, ToCheck):  # checks all items below
		while item is not None:

			if item.childCount() != 0:
				for index in range(item.childCount()):  # will iterate through all children
					parent = item.child(index)
					namehere = parent.text(0)
					if ToCheck == True:
						parent.setCheckState(0, Qt.Checked)
						if parent.childCount() != 0:
							for index in range(parent.childCount()):
								child = parent.child(index)
								namehere = child.text(0)
								child.setCheckState(0, Qt.Checked)
								if child.childCount() != 0:
									for index in range(child.childCount()):
										childschild = child.child(index)
										namehere = childschild.text(0)
										childschild.setCheckState(0, Qt.Checked)
										if childschild.childCount() != 0:
											for index in range(childschild.childCount()):
												childschildschild = childschild.child(index)
												namehere = childschildschild.text(0)
												childschildschild.setCheckState(0, Qt.Checked)
												# child = parent.child(index)
					else:
						parent.setCheckState(0, Qt.Unchecked)
						if parent.childCount() != 0:
							for index in range(parent.childCount()):
								child = parent.child(index)
								namehere = child.text(0)
								child.setCheckState(0, Qt.Unchecked)
								if child.childCount() != 0:
									for index in range(child.childCount()):
										childschild = child.child(index)
										namehere = childschild.text(0)
										childschild.setCheckState(0, Qt.Unchecked)
										if childschild.childCount() != 0:
											for index in range(childschild.childCount()):
												childschildschild = childschild.child(index)
												namehere = childschildschild.text(0)
												childschildschild.setCheckState(0, Qt.Unchecked)

			item = item.child(0)
			# if item is not None:
			#     path.append(str(item.text(0)))

	def CheckUp(self, item, ToCheck):  # checks all items below

		root = self.ui.treeWidget.invisibleRootItem()
		for index in range(root.childCount()):
			fileIs = root.child(index)

		# itemsChecked = 0
		# itemsTotal = 0
		fileIsChecks = 0
		for index in range(fileIs.childCount()):

			project = fileIs.child(index)  # project level
			# while item is not None:  # iterate through project level items


			if project.childCount() != 0:  # if project level exists
				ProjectChecks = 0
				for Pindex in range(project.childCount()):  # will iterate through group level
					Group = project.child(Pindex)  # each group

					if Group.childCount() != 0:  # if there is a subgroup level
						GroupChecks = 0
						for Gindex in range(Group.childCount()):  # will iterate through subgroup level
							SubGroup = Group.child(Gindex)  # each subgroup

							if SubGroup.childCount() != 0:  # if there is a record level
								SubGroupChecks = 0
								ChildCount = SubGroup.childCount()
								for SGindex in range(SubGroup.childCount()):  # will iterate through record level
									Record = SubGroup.child(SGindex)  # each subgroup
									# itemsTotal +=1
									if Record.checkState(0) == Qt.Checked:
										# itemsChecked += 1
										SubGroupChecks += 1
								if SubGroupChecks == 0:
									SubGroup.setCheckState(0, Qt.Unchecked)
								elif SubGroupChecks == ChildCount:
									SubGroup.setCheckState(0, Qt.Checked)
								elif SubGroupChecks < ChildCount:
									SubGroup.setCheckState(0, Qt.PartiallyChecked)

							else:  # if group is lowest level
								# itemsTotal +=1
								if SubGroup.checkState(0) == Qt.Checked:
									# itemsChecked += 1
									GroupChecks += 1

								if GroupChecks == 0:
									Group.setCheckState(0, Qt.Unchecked)
								elif GroupChecks == Group.childCount():
									Group.setCheckState(0, Qt.Checked)
								elif GroupChecks < Group.childCount():
									Group.setCheckState(0, Qt.PartiallyChecked)

					else:  # if group is lowest level
						# itemsTotal +=1
						if Group.checkState(0) == Qt.Checked:
							# itemsChecked += 1
							ProjectChecks += 1
						if ProjectChecks == 0:
							project.setCheckState(0, Qt.Unchecked)
						elif ProjectChecks == project.childCount():
							project.setCheckState(0, Qt.Checked)
						elif ProjectChecks < project.childCount():
							project.setCheckState(0, Qt.PartiallyChecked)

			else:  # project is lowest (record) level
				# itemsTotal +=1
				if project.checkState(0) == Qt.Checked:
					# itemsChecked += 1
					fileIsChecks += 1

				if fileIsChecks == 0:
					fileIs.setCheckState(0, Qt.Unchecked)
				elif fileIsChecks == fileIs.childCount():
					fileIs.setCheckState(0, Qt.Checked)
				elif fileIsChecks < fileIs.childCount():
					fileIs.setCheckState(0, Qt.PartiallyChecked)

		itemsTotal = 0
		itemsChecked = 0



		while item is not None:
			itemsTotal = 0
			itemsChecked = 0
			partiallyChecked = 0
			item = item.parent()

			if item is not None:
				Checktxt = item.text(0)
				for index in range(item.childCount()):  # will iterate through all children
					itemsTotal += 1
					child = item.child(index)
					Checktxt = child.text(0)
					if child.checkState(0) == Qt.Checked:
						itemsChecked += 1
					if child.checkState(0) == Qt.PartiallyChecked:
						partiallyChecked += 1

				if item is not None:
					if itemsTotal == itemsChecked:
						item.setCheckState(0, Qt.Checked)
					elif itemsChecked < itemsTotal:
						if itemsChecked > 0 or partiallyChecked > 0:
							item.setCheckState(0, Qt.PartiallyChecked)
						else:
							item.setCheckState(0, Qt.Unchecked)
							# if item is not None:
							#     path.append(str(item.text(0)))

	def CheckMultiple(self):

		Selected = self.ui.treeWidget.selectedItems()
		# i = 0
		for item in Selected:
			# if i == 0:
			#     currentitemIs = item.text(0)
			# i += 1
			if item.checkState(0) == Qt.Checked:
				item.setCheckState(0, Qt.Unchecked)
			else:
				item.setCheckState(0, Qt.Checked)


				# self.findTreeItem(currentitemIs)

	def clearTreeChecks(self):

		global wasClicked
		wasClicked = False

		value = self.ui.treeWidget.selectedItems()
		currentitemIs = ''
		for item in value:
			currentitemIs = item.text(0)

		self.ui.treeWidget.selectAll()
		Selected = self.ui.treeWidget.selectedItems()
		for item in Selected:
			item.setCheckState(0, Qt.Unchecked)


		self.findTreeItem(currentitemIs)
		wasClicked = True

		return




	def getTreePathDown(self, item):
		path = []
		checkedkids = []
		while item is not None:

			if item.childCount() != 0:
				for index in range(item.childCount()):  # will iterate through all children
					parent = item.child(index)
					namehere = parent.text(0)

			item = item.child(0)
			if item is not None:
				path.append(str(item.text(0)))

		return path  # '/'.join(reversed(path))

	def getTreeChecked(self):
		root = self.ui.treeWidget.invisibleRootItem()
		for index in range(root.childCount()):
			fileIs = root.child(index)

		checkedkids = []
		checkedProjects = []
		checkedGroups = []
		checkedSubGroups = []

		for index in range(fileIs.childCount()):
			project = fileIs.child(index)  # project level
			# while item is not None:  # iterate through project level items


			if project.childCount() != 0:  # if project level exists
				if project.checkState(0) == Qt.Checked:
					checkedProjects.append(project.text(0))

				else:
					for Pindex in range(project.childCount()):  # will iterate through group level
						Group = project.child(Pindex)  # each group
						checkName = Group.text(0)
						# numkid =
						if Group.childCount() != 0:  # if there is a subgroup level
							if Group.checkState(0) == Qt.Checked:
								ProjName = project.text(0)
								GroupName = Group.text(0)
								SetGroup = (ProjName, GroupName)
								# checkedGroups.append(Group.text(0))
								checkedGroups.append(SetGroup)

							else:
								for Gindex in range(Group.childCount()):  # will iterate through subgroup level
									SubGroup = Group.child(Gindex)  # each subgroup
									checkName = SubGroup.text(0)
									if SubGroup.childCount() != 0:  # if there is a record level
										if SubGroup.checkState(0) == Qt.Checked:
											ProjName = project.text(0)
											GroupName = Group.text(0)
											SubGroupName = SubGroup.text(0)
											SetGroup = (ProjName, GroupName, SubGroupName)
											# checkedSubGroups.append(SubGroup.text(0))
											checkedSubGroups.append(SetGroup)
										else:
											for SGindex in range(
													SubGroup.childCount()):  # will iterate through record level
												Record = SubGroup.child(SGindex)  # each subgroup
												checkName = Record.text(0)
												if Record.checkState(0) == Qt.Checked:
													checkedkids.append(Record.text(0))


									else:  # if group is lowest level
										# for SGindex in range(SubGroup.childCount()):  # for each project if lowest level item
										#      Ritem = SubGroup.child(SGindex)
										if SubGroup.checkState(0) == Qt.Checked:
											checkedkids.append(SubGroup.text(0))

						else:  # if group is lowest level
							if Group.checkState(0) == Qt.Checked:
								checkedkids.append(Group.text(0))

			else:  # project is lowest (record) level
				if project.checkState(0) == Qt.Checked:
					checkedkids.append(project.text(0))

		return checkedProjects, checkedGroups, checkedSubGroups, checkedkids



	@pyqtSlot()
	def on_actionGL_triggered(self):
		global GLMsg
		if self.ui.actionGL.isChecked() == True:
			if GLMsg == True:
				question = 'Select a sequence to use for the predicted germline in a multiple alignment.'
				buttons = 'OK'
				answer = informationMessage(self, question, buttons)

	def ModuleFind(self):
		from modulefinder import ModuleFinder
		from VGenesDialogues import openFile
		import os

		Pathname = openFile(self, 'CSV')

		workingdir, filename = os.path.split(Pathname)

		os.chdir(workingdir)
		# import filename

		finder = ModuleFinder()
		finder.run_script(filename)

		Doc = 'Loaded modules:\n'
		for name, mod in finder.modules.items():
			# Doc('%s: ' % name, end=''\n)
			Doc += (','.join(list(mod.globalnames.keys())[:3]))
			Doc += '\n'
			# Doc
		Doc += ('-' * 50)

	# Doc +=('Modules not imported:')
	# Doc +=('\n'.join(finder.badmodules.keys()))


	@pyqtSlot()
	def on_actionPrint_triggered(self):

		FontIs = self.TextEdit.textEdit.currentFont()
		font = QFont(FontIs)
		if self.ui.tabWidget.currentIndex() == 0:
			fields = ['SeqName', 'V1', 'D1', 'J1', 'VLocus', 'JLocus', 'productive', 'TotMut', 'CDR3DNA', 'CDR3AA',
			          'CDR3Length', 'CDR3pI', 'ClonalPool', 'Isotype', 'Sequence']
			# SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
			DataIs = [data[0], data[3], data[6], data[9], data[90], data[91], data[14], data[57], data[81], data[82],
			          str(data[83]), str(data[100]), str(data[88]), data[101],
			          data[79]]  # VGenesSQL.RunSQL(DBFilename, SQLStatement)
			NameList = []
			for item in fields:
				NameList.append(str(self.TransLateFieldtoReal(item, False)))
			NameList[7] = 'Mutations'
			NameList[11] = 'CDR3 pI'

			Document = ''
			i = 0
			for item in DataIs:
				Document += (NameList[i] + ': \t' + str(item) + '\n')
				i += 1
			Document += ('Protein: ' + self.ui.txtAASeq.toPlainText() + '\n')
			Document += '\n'
			Document += self.windowTitle()
			font.setPointSize(10)
			font.setFamily('Lucida Grande')


		if self.ui.tabWidget.currentIndex() == 2:
			Document = data[0] + '\n'
			Document += ('DNA: ' + self.ui.txtDNASeq.toPlainText() + '\n')
			Document += ('Protein: ' + self.ui.txtAASeq.toPlainText() + '\n')
			Document += ('\n' + self.windowTitle())
			font.setPointSize(10)
			font.setFamily('Lucida Grande')

		elif self.ui.tabWidget.currentIndex() == 3:
			Document = data[0] + '\n'
			Document += self.ui.txtSeqAlignment.toPlainText()
			Document += ('\n' + self.windowTitle())
			font.setPointSize(7)
			font.setFamily('Courier New')

		self.TextEdit.textEdit.setFont(font)

		# self.TextEdit.show()

		self.TextEdit.textEdit.setText(Document)

		document = self.TextEdit.textEdit.document()
		printer = QPrinter()

		dlg = QPrintDialog(printer, self)
		if dlg.exec_() != QtWidgets.QDialog.Accepted:
			return

		if self.ui.tabWidget.currentIndex() == 3: printer.setOrientation(QPrinter.Landscape)
		document.print_(printer)

		self.statusBar().showMessage("Ready", 2000)

	@pyqtSlot()
	def on_btnDeleteRecord_clicked(self):
		self.on_actionDelete_record_triggered()

	@pyqtSlot()
	def on_actionDelete_record_triggered(self):
		fields = ['SeqName', 'ID']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		if SQLStatement == 'None':
			question = 'Check sequences to be deleted'
			buttons = 'OK'
			answer = informationMessage(self, question, buttons)
		else:
			question = 'Are you certain you want to delete the checked sequences?'
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)

		if answer == 'Yes':
			VGenesSQL.deleterecords(DBFilename, SQLStatement)
			if DBFilename != None:
				if os.path.isfile(DBFilename):
					self.LoadDB(DBFilename)

	@pyqtSlot()
	def on_actionFind_Clonal_triggered(self):
		# self.ShowProgressBar()
		from operator import itemgetter
		import itertools
		remove = False
		items = ('Clonal Pools', 'Annotate Duplicates', 'Remove Duplicates')
		title = 'Choose analysis:'
		item = setItem(self, items, title)
		if item == 'Cancel':
			return
		elif item[:3] == 'Ann':
			Duplicates = True
		elif item[:3] == 'Rem':
			Duplicates = True
			query = 'This function will delete all but one duplicated sequences, annotating ' \
			        'the remaining with depth (in the "Quality" field) and the names of ' \
			        'duplicated sequences (in the "Comments" field), \nProceed with delete (Yes), or just annotate (No)'
			answer = questionMessage(self, query, 'YNC')
			if answer == 'Yes':
				remove = True
			elif answer == 'Cancel':
				return
		else:
			Duplicates = False

		if self.ui.cboFindField.currentText() == 'Name': self.ui.cboFindField.setCurrentText('Project')
		# SeqName, Sequence, ClonalPool, GermlineSequence, Mutations
		answer = questionMessage(self,
		                         'Use a field to delineate multiple subjects (default = "Project")?\n\n "No" will process all selected as one subject.\n\n Press "Cancel" to choose field in the search panel before running analysis.',
		                         'YNC')
		if answer == 'Yes':
			field = self.ui.cboFindField.currentText()
			fieldsearch = self.TransLateFieldtoReal(field, True)
		elif answer == 'No':
			fieldsearch = 'None'
		elif answer == 'Cancel':
			return
		if fieldsearch == 'None':
			fields = ['SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence',
			          'ID']
		else:
			fields = ['SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence',
			          'ID', fieldsearch]

		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
		DataIs2 = []
		ErLog = ''
		Errs = 0
		for item in DataIs:
			SeqName = item[0]
			self.findTreeItem(SeqName)
			try:
				if int(data[83]) != 0:
					VGenesSQL.UpdateField(data[119], str(0), 'ClonalPool', DBFilename)
			except:
				print('stop')

			if int(item[6]) > 0:  # has CDR3 scored
				DataIs2.append(item)
			else:
				ErLog += SeqName + '\n'
				Errs += 1

		TotSeqs = len(DataIs2)
		# if fieldsearch == 'None':
		#     DataIs2.sort(key=itemgetter(3))
		# else:
		#     DataIs2.sort(key=itemgetter(17))
		ClonalPool = []
		ClonalPools = []

		# fields = ['SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence', 'ID', fieldsearch]

		if answer == 'Yes':
			DataIs2.sort(key=itemgetter(10))
			for k, v in itertools.groupby(DataIs2, key=itemgetter(10)):  # first split out seperate clonal pools
				# i = int(k)

				if len(k) != 0:
					for item in v:
						ClonalPool.append(item)
					CurrentPool = tuple(ClonalPool)
					ClonalPools.append(CurrentPool)
					ClonalPool.clear()
		else:
			ClonalPools.append(DataIs2)

		Currentrecord = self.ui.txtName.toPlainText()

		CPseqs = 0
		CPs = 0

		for pool in ClonalPools:
			Pool = list(pool)
			# result = VGenesSeq.Intraclonal(Pool)

			# if len(result) > 0:



			# CPList, DuplicateList = VGenesCloneCaller.CloneCaller(Pool)  # should return list of tuple with SeqNames in each Clonal Pool

			CPList = VGenesCloneCaller.CloneCaller(Pool, Duplicates)
			i = 1

			for record in CPList:
				CPs += 1
				j = 1
				DupList = 'Sequences identical: '
				for item in record:

					self.findTreeItem(item)
					if Duplicates == False:
						VGenesSQL.UpdateField(data[119], str(i), 'ClonalPool', DBFilename)
					else:
						if j == 1:
							SeqName = 'Duplicate of:  ' + item
							FirstOne = data[119]
						else:
							if remove == False:
								VGenesSQL.UpdateField(data[119], SeqName, 'Quality', DBFilename)
							else:
								VGenesSQL.UpdateField(data[119], 'Duplicate', 'Quality', DBFilename)
							DupList += (item + ', ')

						j += 1
					CPseqs += 1
				depth = 'Depth = ' + str(j - 1)
				if Duplicates == True:
					self.findTreeItem(FirstOne)
					if DupList[(len(DupList) - 2):] == ', ':
						DupList = DupList[:(len(DupList) - 2)]
					if data[94] != ' ' or data[94] != 'Comments':
						DupList = DupList + ', ' + data[94]

					if data[95] != ' ' or data[95] != 'Quality':
						depth = depth + '  ' + data[95]

					VGenesSQL.UpdateField(FirstOne, DupList, 'Comments', DBFilename)
					VGenesSQL.UpdateField(FirstOne, depth, 'Quality', DBFilename)
				i += 1

		model = self.ui.tableView.model()
		model.refresh()

		if answer == 'Yes':
			if DBFilename != None:
				if os.path.isfile(DBFilename):
					self.LoadDB(DBFilename)
		else:
			if self.ui.cboTreeOp1.currentText() == 'Clonal Pool' or self.ui.cboTreeOp2.currentText() == 'Clonal Pool' or self.ui.cboTreeOp3.currentText() == 'Clonal Pool':
				self.on_btnUpdateTree_clicked()

		self.findTreeItem(Currentrecord)
		ErLog2 = str(CPs) + ' clonal pools containing ' + str(CPseqs) + ' sequences were identified from ' + str(
			TotSeqs) + ' total sequences analyzed.\n'

		if remove == True:
			# self.clearTreeChecks()
			self.LoadDB(DBFilename)
			self.ui.txtFieldSearch.setPlainText('Duplicate')
			self.ui.cboFindField.setCurrentText('Quality')
			done = self.on_btnFieldSearch_clicked()
			done = self.on_actionDelete_record_triggered()

		if len(ErLog2) > 0:
			Erlog2 = ErLog2 + 'The following ' + str(
				Errs) + ' sequences could not be anaylzed for\nclonality because no CDR3s are indicated:\n' + ErLog
			ErlogFile = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'IgBlast', 'database',
			                         'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs

			with open(ErlogFile, 'w') as currentFile:
				currentFile.write(Erlog2)

			self.ShowVGenesText(ErlogFile)

			# VGenesSQL.UpdateField(DataRow,ItemValue,FieldName,DBFilename)
			# model.refresh()

	@pyqtSlot()
	def on_actionAnalyze_Mutations_triggered(self):
		import VMapHotspots
		# setItem(self, items, title):
		items = ('Intraclonal Diversity', 'HS-Summary', "Standard", "Hotspots", 'Mutation Frequencies', 'Heat Map',
		         'AnalysisDB Heat Map')
		title = 'Choose report type:'
		item = setItem(self, items, title)


		if item == "Cancel":
			return
		if item == 'Intraclonal Diversity':
			self.AnalyzeMutations()
		else:
			VMapHotspots.MapHotspots(self, item, DBFilename, data[0])


	@pyqtSlot()
	def AnalyzeMutations(self):

		if self.ui.cboFindField.currentText() == 'Name': self.ui.cboFindField.setCurrentText('Project')
		# SeqName, Sequence, ClonalPool, GermlineSequence, Mutations
		answer = questionMessage(self,
		                         'Use a field to delineate multiple subjects (default = "Project")?\n\n "No" will process all selected as one subject.\n\n Press "Cancel" to choose field in the search panel before running analysis.',
		                         'YNC')
		if answer == 'Yes':
			field = self.ui.cboFindField.currentText()
			fieldsearch = self.TransLateFieldtoReal(field, True)
		elif answer == 'No':
			fieldsearch = 'None'
		elif answer == 'Cancel':
			return

		if fieldsearch == 'None':
			fields = ['SeqName', 'Sequence', 'ClonalPool', 'GermlineSequence', 'Mutations', 'GVbeg', 'GVend', 'Species',
			          'V1']

		else:
			fields = ['SeqName', 'Sequence', 'ClonalPool', 'GermlineSequence', 'Mutations', 'GVbeg', 'GVend', 'Species',
			          'V1', fieldsearch]

		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first

		if answer == 'Yes':
			DataIn.sort(key=itemgetter(10))
		else:
			DataIn.sort(key=itemgetter(2))
		ClonalPool = []
		ClonalPools = []

		if answer == 'Yes':
			for k, v in itertools.groupby(DataIn, key=itemgetter(10)):  # first split out seperate clonal pools
				# i = int(k)

				if len(k) != 0:
					for item in v:
						ClonalPool.append(item)
					CurrentPool = tuple(ClonalPool)
					ClonalPools.append(CurrentPool)
					ClonalPool.clear()
		else:
			ClonalPools.append(DataIn)

		f = saveFile(self, 'CSV')
		with open(f, 'w') as currentfile:
			doc = 'Clonotype, Sequence 1, Sequence 2, Differences, R-Differences, S-Differences, Begin, End, Length,  Matches, Adjusted Matches, ' \
			      'Adjusted Differences, R-Matches, Adjusted R-Matches, Adjusted R-Diferences, ' \
			      'S-Matches, Adjusted S-Matches, Adjusted S-Differences\n'

			# CP, Name1, Name2, beg, end, lengthCompared, TotMatch, int(TotAdjMatch), TotDif, int(TotAdjDif), len(
			# 	RMatches),
			# int(RAdjMatch), len(RDifferences), int(RAdjDif), len(SMatches), int(SAdjMatch), len(SDifferences),
			# int(SAdjDif))


		for pool in ClonalPools:
			Pool = list(pool)
			result = VGenesSeq.Intraclonal(Pool, DBFilename)

			if len(result) > 0:

				with open(f, 'a') as currentfile:
					if answer == 'Yes':
						header = str(Pool[0][7]) + '\n'
						doc += header

					for item in result:
						for i in range(0, 18):
							doc += (str(item[i]) + ', ')
						doc += '\n'
					currentfile.write(doc)

		self.ShowVGenesText(f)

	@pyqtSlot()
	def on_actionCreateAnalysisDB_triggered(self):
		filename = openFile(self, 'FASTA')
		VGenesSQL.CreateAnalysisDB(filename, DBFilename)

	@pyqtSlot()
	def on_actionMultiple_Alignement_triggered(self):
		self.AlignSequences('none')

	@pyqtSlot()
	def on_actionImport_Vgenes_database_triggered(self):
		filename = openFile(self, 'db')
		VGenesSQL.ImportVDB(filename, DBFilename)

	@pyqtSlot()
	def on_actionIsotypes_triggered(self):
		DataIs = []

		DataIs.append(('IgMSeq',
		               'GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGGCCGTTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTGTCCTGGAAATACAAGAACAACTCTGACATCAGCAGTACCCGGGGCTTCCCATCAGTCCTGAGAGGGGGCAAGTACGCAGCCACCTCACAGGTGCTGCTGCCTTCCAAGGACGTCATGCAGGGCACAGACGAACACGTGGTGTGCAAAGTCCAGCACCCCAACGGCAACAAAGAAAAGAACGTGCCTCTTCCAG'))

		DataIs.append(('IgG1Seq',
		               'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACATCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAAAGTTGAGCCCAAATCTTGTGACAAAACTCACACATGCCCACCGTGCCCAGCACCTGAACTCCTGGGGGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACATGCGTGGTGGTGGACGTGAGCCACGAAGACCCTGAGGTCAAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTACCGGGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAATGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGATGAGCTGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTATCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'))

		DataIs.append(('IgG2Seq',
		               'CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCAGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGAGCGCAAATGTTGTGTCGAGTGCCCACCGTGCCCAGCACCACCTGTGGCAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGACACCCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCACGGGAGGAGCAGTTCAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTTGTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGGCAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACACCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'))

		DataIs.append(('IgG3Seq',
		               'CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACCCAGACCTACACCTGCAACGTGAATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGCTCAAAACCCCACTTGGTGACACAACTCACACATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCATGCCCACGGTGCCCAGAGCCCAAATCTTGTGACACACCTCCCCCGTGCCCAAGGTGCCCAGCACCTGAACTCCTGGGAGGACCGTCAGTCTTCCTCTTCCCCCCAAAACCCAAGGATACCCTTATGATTTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCACGAAGACCCCGAGGTCCAGTTCAAGTGGTACGTGGACGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTACAACAGCACGTTCCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGCCCTCCCAGCCCCCATCGAGAAAACCATCTCCAAAACCAAAGGACAGCCCCGAGAACCACAGGTGTACACCCTGCCCCCATCCCGGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAGCGGGCAGCCGGAGAACAACTACAACACCACGCCTCCCATGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAAGCTCACCGTGGACAAGAGCAGGTGGCAGCAGGGGAACATCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCGCTTCACGCAGAAGAGCCTCTCCCTGTCTCCGGGTAAATGA'))

		DataIs.append(('IgG4Seq',
		               'CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCCCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAGCTTGGGCACGAAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGAGAGTTGAGTCCAAATATGGTCCCCCATGCCCATCATGCCCAGCACCTGAGTTCCTGGGGGGACCATCAGTCTTCCTGTTCCCCCCAAAACCCAAGGACACTCTCATGATCTCCCGGACCCCTGAGGTCACGTGCGTGGTGGTGGACGTGAGCCAGGAAGACCCCGAGGTCCAGTTCAACTGGTACGTGGATGGCGTGGAGGTGCATAATGCCAAGACAAAGCCGCGGGAGGAGCAGTTCAACAGCACGTACCGTGTGGTCAGCGTCCTCACCGTCCTGCACCAGGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTCTCCAACAAAGGCCTCCCGTCCTCCATCGAGAAAACCATCTCCAAAGCCAAAGGGCAGCCCCGAGAGCCACAGGTGTACACCCTGCCCCCATCCCAGGAGGAGATGACCAAGAACCAGGTCAGCCTGACCTGCCTGGTCAAAGGCTTCTACCCCAGCGACATCGCCGTGGAGTGGGAGAGCAATGGGCAGCCGGAGAACAACTACAAGACCACGCCTCCCGTGCTGGACTCCGACGGCTCCTTCTTCCTCTACAGCAGGCTAACCGTGGACAAGAGCAGGTGGCAGGAGGGGAATGTCTTCTCATGCTCCGTGATGCATGAGGCTCTGCACAACCACTACACACAGAAGAGCCTCTCCCTGTCTCTGGGTAAATGA'))

		DataIs.append(('IgA1seq',
		               'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGGGCGTGACCGCCAGAAACTTCCCACCCAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCTAGCCGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCTCAACTCCACCTACCCCATCTCCCTCAACTCCACCTACCCCATCTCCCTCATGCTGCCACCCCCGACTGTCACTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCAGGTGTCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCGGGCTGTGCCGAGCCATGGAACCATGGGAAGACCTTCACTTGCACTGCTGCCTACCCCGAGTCCAAGACCCCGCTAACCGCCACCCTCTCAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGCGGCTTCAGCCCCAAGGACGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'))

		DataIs.append(('IgA2seq',
		               'CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAGCGGACAGAACGTGACCGCCAGAAACTTCCCACCTAGCCAGGATGCCTCCGGGGACCTGTACACCACGAGCAGCCAGCTGACCCTGCCGGCCACACAGTGCCCAGACGGCAAGTCCGTGACATGCCACGTGAAGCACTACACGAATCCCAGCCAGGATGTGACTGTGCCCTGCCCAGTTCCCCCACCTCCCCCATGCTGCCACCCCCGACTGTCGCTGCACCGACCGGCCCTCGAGGACCTGCTCTTAGGTTCAGAAGCGAACCTCACGTGCACACTGACCGGCCTGAGAGATGCCTCTGGTGCCACCTTCACCTGGACGCCCTCAAGTGGGAAGAGCGCTGTTCAAGGACCACCTGAGCGTGACCTCTGTGGCTGCTACAGCGTGTCCAGTGTCCTGCCTGGCTGTGCCCAGCCATGGAACCATGGGGAGACCTTCACCTGCACTGCTGCCCACCCCGAGTTGAAGACCCCACTAACCGCCAACATCACAAAATCCGGAAACACATTCCGGCCCGAGGTCCACCTGCTGCCGCCGCCGTCGGAGGAGCTGGCCCTGAACGAGCTGGTGACGCTGACGTGCCTGGCACGTGGCTTCAGCCCCAAGGATGTGCTGGTTCGCTGGCTGCAGGGGTCACAGGAGCTGCCCCGCGAGAAGTACCTGACTTGGGCATCCCGGCAGGAGCCCAGCCAGGGCACCACCACCTTCGCTGTGACCAGCATACTGCGCGTGGCAGCCGAGGACTGGAAGAAGGGGGACACCTTCTCCTGCATGGTGGGCCACGAGGCCCTGCCGCTGGCCTTCACACAGAAGACCATCGACCGCTTGGCGGGTAAACCCACCCATGTCAATGTGTCTGTTGTCATGGCGGAGGTGGACGGCACCTGCTACTGA'))

		DataIs.append(('IgEseq',
		               'CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCTGCAAAAACATTCCCTCCAATGCCACCTCCGTGACTCTGGGCTGCCTGGCCACGGGCTACTTCCCGGAGCCGGTGATGGTGACCTGCGACACAGGCTCCCTCAACGGGACAACTATGACCTTACCAGCCACCACCCTCACGCTCTCTGGTCACTATGCCACCATCAGCTTGCTGACCGTCTCGGGTGCGTGGGCCAAGCAGATGTTCACCTGCCGTGTGGCACACACTCCATCGTCCACAGACTGGGTCGACAACAAAACCTTCAGCGTCTGCTCCAGGGACTTCACCCCGCCCACCGTGAAGATCTTACAGTCGTCCTGCGACGGCGGCGGGCACTTCCCCCCGACCATCCAGCTCCTGTGCCTCGTCTCTGGGTACACCCCAGGGACTATCAACATCACCTGGCTGGAGGACGGGCAGGTCATGGACGTGGACTTGTCCACCGCCTCTACCACGCAGGAGGGTGAGCTGGCCTCCACACAAAGCGAGCTCACCCTCAGCCAGAAGCACTGGCTGTCAGACCGCACCTACACCTGCCAGGTCACCTATCAAGGTCACACCTTTGAGGACAGCACCAAGAAGTGTGCAGATTCCAACCCGAGAGGGGTGAGCGCCTACCTAAGCCGGCCCAGCCCGTTCGACCTGTTCATCCGCAAGTCGCCCACGATCACCTGTCTGGTGGTGGACCTGGCACCCAGCAAGGGGACCGTGAACCTGACCTGGTCCCGGGCCAGTGGGAAGCCTGTGAACCACTCCACCAGAAAGGAGGAGAAGCAGCGCAATGGCACGTTAACCGTCACGTCCACCCTGCCGGTGGGCACCCGAGACTGGATCGAGGGGGAGACCTACCAGTGCAGGGTGACCCACCCCCACCTGCCCAGGGCCCTCATGCGGTCCACGACCAAGACCAGCGGCCCGCGTGCTGCCCCGGAAGTCTATGCGTTTGCGACGCCGGAGTGGCCGGGGAGCCGGGACAAGCGCACCCTCGCCTGCCTGATCCAGAACTTCATGCCTGAGGACATCTCGGTGCAGTGGCTGCACAACGAGGTGCAGCTCCCGGACGCCCGGCACAGCACGACGCAGCCCCGCAAGACCAAGGGCTCCGGCTTCTTCGTCTTCAGCCGCCTGGAGGTGACCAGGGCCGAATGGGAGCAGAAAGATGAGTTCATCTGCCGTGCAGTCCATGAGGCAGCGAGCCCCTCACAGACCGTCCAGCGAGCGGTGTCTGTAAATCCCGGTAAATGA'))

		DataIs.append(('IgDseq',
		               'CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGGCATGCTTGATAACTGGGTACCACCCAACGTCCGTGACTGTCACCTGGTACATGGGGACACAGAGCCAGCCCCAGAGAACCTTCCCTGAGATACAAAGACGGGACAGCTACTACATGACAAGCAGCCAGCTCTCCACCCCCCTCCAGCAGTGGCGCCAAGGCGAGTACAAATGCGTGGTCCAGCACACCGCCAGCAAGAGTAAGAAGGAGATCTTCCGCTGGCCAGAGTCTCCAAAGGCACAGGCCTCCTCCGTGCCCACTGCACAACCCCAAGCAGAGGGCAGCCTCGCCAAGGCAACCACAGCCCCAGCCACCACCCGTAACACAGGAAGAGGAGGAGAAGAGAAGAAGAAGGAGAAGGAGAAAGAGGAACAAGAAGAGAGAGAGACAAAGACACCAGAGTGTCCGAGCCACACCCAGCCTCTTGGCGTCTACCTGCTAACCCCTGCAGTGCAGGACCTGTGGCTCCGGGACAAAGCCACCTTCACCTGCTTCGTGGTGGGCAGTGACCTGAAGGATGCTCACCTGACCTGGGAGGTGGCTGGGAAGGTCCCCACAGGGGGCGTGGAGGAAGGGCTGCTGGAGCGGCACAGCAACGGCTCCCAGAGCCAGCACAGCCGTCTGACCCTGCCCAGGTCCTTGTGGAACGCGGGGACCTCCGTCACCTGCACACTGAACCATCCCAGCCTCCCACCCCAGAGGTTGATGGCGCTGAGAGAACCCGCTGCGCAGGCACCCGTCAAGCTTTCTCTGAACCTGCTGGCCTCGTCTGACCCTCCCGAGGCGGCCTCGTGGCTCCTGTGTGAGGTGTCTGGCTTCTCGCCCCCCAACATCCTCCTGATGTGGCTGGAGGACCAGCGTGAGGTGAACACTTCTGGGTTTGCCCCCGCACGCCCCCCTCCACAGCCCAGGAGCACCACGTTCTGGGCCTGGAGTGTGCTGCGTGTCCCAGCCCCGCCCAGCCCTCAGCCAGCCACCTACACGTGTGTGGTCAGCCACGAGGACTCCCGGACTCTGCTCAACGCCAGCCGGAGCCTAGAAGTCAGCTACCTGGCCATGACCCCCCTGATCCCTCAGAGCAAGGATGAGAACAGCGATGACTACACGACCTTTGATGATGTGGGCAGCCTGTGGACCACCCTGTCCACGTTTGTGGCCCTCTTCATCCTCACCCTCCTCTACAGCGGCATTGTCACTTTCATCAAGGTGAAGTAG'))

		self.AlignSequences(DataIs)

	def AlignSequences(self, DataIn):

		QApplication.setOverrideCursor(Qt.WaitCursor)

		if DataIn == 'none':
			fields = ['SeqName', 'Sequence']
			SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
			DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
			global GLMsg
			if len(DataIs) == 1:
				GLMsg = False
				self.ui.actionGL.setChecked(True)
				GLMsg = True
				GermSeq = data[80]
				Germline = ('Germline', GermSeq)
				DataIs.append(Germline)
			else:
				if self.ui.actionGL.isChecked() == True:
					GLMsg = True
					GermSeq = data[80]
					Germline = ('Germline', GermSeq)
					DataIs.append(Germline)

		elif DataIn == 'edit':
			DataIn = 'none'
			DataIs = []


			SeqName = data[0]

			DNAseq = self.ui.txtDNASeq.toPlainText()
			Sequence = (SeqName, DNAseq)
			DataIs.append(Sequence)

			global GLMsg
			if len(DataIs) == 1:
				GLMsg = False
				self.ui.actionGL.setChecked(True)
				GLMsg = True
				GermSeq = data[80]
				Germline = ('Germline', GermSeq)
				DataIs.append(Germline)

		else:
			self.ui.actionGL.setChecked(False)
			DataIs = DataIn


		VGenesSeq.ClustalO(DataIs, 80, True)

		outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega',
		                           'my-out-seqs.fa')  # '/Applications/ClustalOmega/my-out-seqs.fa'
		lenName = 0
		longestName = 0
		alignmentText = ''
		germseq = ''
		germpeptide = ''

		each = ()
		all = []
		if self.ui.actionGL.isChecked() == False:
			# each.append('Consensus: ')
			longestName = 11
		else:
			# each.append('Germline: ')
			longestName = 10

		# each.append('') #DNA
		# each.append('') #AA
		peptide = ''
		SeqName = ''
		StartAll = False
		if os.path.isfile(outfilename):
			with open(outfilename, 'r') as currentfile:
				for line in currentfile:
					Readline = line.replace('\n', '').replace('\r', '').replace('-', '.')
					Readline = Readline.strip()
					if Readline[0] == '>':
						if StartAll == True:
							all.append(each)
						StartAll = True
						# each.clear()
						SeqName = Readline[1:] + ':'
						lenName = len(SeqName)
						if lenName > longestName:
							longestName = lenName + 2

							# if SeqName != 'Germline':
							#     each.append(SeqName)

					else:

						if self.ui.actionAA.isChecked() == True:
							AASeq, ErMessage = VGenesSeq.Translator(Readline, 0)
							peptide = ''
							for res in AASeq:
								peptide += (' ' + res + ' ')
						peptide = peptide[0:len(Readline)]
						if SeqName != 'Germline:':
							each = (SeqName, Readline, peptide)
							# each.append(Readline)
							# each.append(peptide)
						else:
							germseq = Readline
							germpeptide = peptide
							StartAll = False
			if StartAll == True:
				all.append(each)
		else:
			return
		# todo add header that says what germline based on
		if self.ui.actionGL.isChecked() == True:
			consensusDNA = germseq
			consensusAA = germpeptide

		else:
			firstOne = all[1]
			seqlen = len(firstOne[1])
			if self.ui.actionDNA.isChecked() == True:
				consensusDNA = ''
				tester = ''
				# testl = []
				for i in range(0, seqlen - 1):
					tester = ''
					Cnuc = ''
					for item in all:
						seq = item[1]
						tester += seq[i]

					frequencies = [(c, tester.count(c)) for c in set(tester)]
					Cnuc = max(frequencies, key=lambda x: x[1])[0]
					consensusDNA += Cnuc

			if self.ui.actionAA.isChecked() == True:
				consensusAA = ''
				tester = ''
				firstOne = all[1]
				seqlen = len(firstOne[1])
				# testl = []
				for i in range(0, seqlen - 1):
					tester = ''
					Caa = ''
					for item in all:
						seq = item[2]
						tester += seq[i]

					frequencies = [(c, tester.count(c)) for c in set(tester)]
					Caa = max(frequencies, key=lambda x: x[1])[0]
					consensusAA += Caa


					# need build numberring lines also
					# first record is germline or consensus whatever used and empty seq and AA
					# need to use ones produced above
					# also longestName is longest and need code to ensure all that long with ': '
					# build alignment with name and 50 per

		header = 'VGenes multiple alignment using Clustal Omega. \n'
		if self.ui.actionGL.isChecked() == False:
			ConName = 'Consensus: '

		else:
			ConName = 'Germline: '
			header += 'Alignment relative to the predicted germline gene for ' + data[0] + '.\n'

		while len(ConName) < longestName:
			ConName += ' '

		AASpaces = ''
		while len(AASpaces) < longestName:
			AASpaces += ' '

		if self.ui.actionDNA.isChecked() == False and self.ui.actionAA.isChecked() == False:
			self.ui.actionDNA.setChecked(True)

		alignmentText = header
		i = 0
		endSeg = 0
		done = False
		ConAdd = True

		# for j in range[0,longestName]:
		#     AASpaces += ' '
		if self.ui.actionDNA.isChecked() == True:
			maxLen = len(consensusDNA)
		else:
			NewConAA = consensusAA.replace(' ', '')

			# canAA = False
			maxLen = len(NewConAA)

		# canAA = True
		while endSeg <= maxLen - 1:
			if i + 60 < maxLen:
				# if i == 0:
				#     endSeg = 49
				# else:
				endSeg = i + 60
			else:
				endSeg = maxLen

			for seq in all:
				SeqName = seq[0]
				DNASeq = seq[1]
				AASeq = seq[2]
				NewAA = AASeq.replace(' ', '')
				while len(SeqName) < longestName:
					SeqName += ' '
				# todo can build num line even add CDR if align relative to germline instead just number as end
				toSpace = len(str(maxLen))
				endLabel = str(endSeg)
				while len(endLabel) < toSpace:
					endLabel += ' '
				endLabel = '  ' + endLabel

				if self.ui.actionDNA.isChecked() == True:

					ConSegDNA = consensusDNA[i:endSeg]
					DNASeqSeg = DNASeq[i:endSeg]
					DNAArt = ''
					for n in range(0, len(ConSegDNA)):
						if DNASeqSeg[n] == ConSegDNA[n]:
							if DataIn == 'none':
								DNAArt += '-'
							else:
								char = DNASeqSeg[n]
								char = char.upper()
								DNAArt += char
						else:
							if DataIn == 'none':
								DNAArt += DNASeqSeg[n]
							else:
								char = DNASeqSeg[n]
								char = char.lower()
								DNAArt += char

					ConSegDNA = ConName + ConSegDNA + endLabel
					DNASeqSeg = SeqName + DNAArt + endLabel
					if self.ui.actionAA.isChecked() == True:
						AArt = ''
						ConSegAA = consensusAA[i:endSeg]
						AASeqSeg = AASeq[i:endSeg]

						for n in range(0, len(ConSegAA)):
							if AASeqSeg[n] == ConSegAA[n]:
								AArt += ' '
							else:
								AArt += AASeqSeg[n]

						AASeqSeg = AASpaces + AArt  # + endLabel
						ConSegAA = AASpaces + ConSegAA
						if ConAdd == True:
							alignmentText += '\n' + ConSegAA + '\n'
							alignmentText += ConSegDNA + '\n'
							ConAdd = False
						alignmentText += AASeqSeg + '\n'
						alignmentText += DNASeqSeg + '\n'
					else:
						if ConAdd == True:
							alignmentText += '\n' + ConSegDNA + '\n'
							ConAdd = False
						alignmentText += DNASeqSeg + '\n'

				else:
					if self.ui.actionAA.isChecked() == True:
						AArt = ''
						ConSegAA = NewConAA[i:endSeg]
						AASeqSeg = NewAA[i:endSeg]

						for n in range(0, len(ConSegAA)):
							if AASeqSeg[n] == ConSegAA[n]:
								AArt += '-'
							else:
								AArt += AASeqSeg[n]

						AASeqSeg = SeqName + AArt + endLabel
						ConSegAA = ConName + ConSegAA
						if ConAdd == True:
							alignmentText += '\n' + ConSegAA + '\n'

							ConAdd = False
						alignmentText += AASeqSeg + '\n'

			i += 60
			ConAdd = True
			alignmentText += '\n'




		Style = 'aligned'

		self.ShowVGenesTextEdit(alignmentText, Style)


		QApplication.restoreOverrideCursor()

		# for item in Aligned:


		print('done')

	# do lengthy process


	def getTreeSelected(self):
		root = self.ui.treeWidget.invisibleRootItem()
		ListSelected = []
		for item in self.ui.treeWidget.selectedItems():
			ListSelected.append(str(item.text(0)))
			# (item.parent() or root).removeChild(item)

		return ListSelected

	@pyqtSlot()
	def TreeSelectChanged(self):



		value = self.ui.treeWidget.selectedItems()

		model = self.ui.tableView.model()
		while model.canFetchMore():
			model.fetchMore()

		name = ''
		# name2 = ''
		for item in value:
			name = item.text(0)
			# name2 = item.text(1)

		FieldCheck = self.FieldChangeCheck()
		# if FieldCheck == 'exit':
		#     return

		try:
			MatchingIndex = NameIndex[name]

			self.DialScroll(MatchingIndex)
		except:
			return
			# self.findTableViewRecord(name)
			# print(role)

	def MatchingValue(self, IndexIs):
		try:
			return list(NameIndex.keys())[list(NameIndex.values()).index(int(IndexIs))]
		except:
			return 'None'

	def TreeviewOptions(self):

		Option1 = self.ui.cboTreeOp1.currentText()
		self.ui.cboTreeOp1.setToolTip('Press update tree to implement changes')
		if Option1 == 'None':
			self.ui.cboTreeOp2.setCurrentText('None')
			# self.ui.cboTreeOp3.setCurrentText('None')
		Option2 = self.ui.cboTreeOp2.currentText()
		self.ui.cboTreeOp2.setToolTip('Press update tree to implement changes')
		Option3 = self.ui.cboTreeOp3.currentText()
		if Option2 == 'None':
			self.ui.cboTreeOp3.setCurrentText('None')
		self.ui.cboTreeOp3.setToolTip('Press update tree to implement changes')

		SqlStatement = 'SELECT DISTINCT '

	def ApplicationStarted(self):
		StartUpOptions = StartUpDialogue()
		StartUpOptions.exec_()

	@pyqtSlot()
	def StartUpClicked(self):
		self.show()

		global DBFilename
		if StartUpAnswer == 'New':
			self.on_action_New_triggered()

		elif StartUpAnswer == 'Open':
			self.on_action_Open_triggered()

		elif StartUpAnswer[0:6] == 'Recent':

			DBFilename = StartUpAnswer[6:]
			# Vgenes.on_action_Open_triggered()
			# self.LoadDB(DBFilename)

			# Vgenes.InputSeqQuery(self)
			if os.path.isfile(DBFilename):
				self.LoadDB(DBFilename)
			else:
				VGenesSQL.creatnewDB(DBFilename)
			self.UpdateRecentList(DBFilename, True)

		self.SaveBackup()
		# self.EditableSqlModel.refresh()

	# @pyqtSlot(int)
	# def on_spinBox_valueChanged(self, value): #how to handle spin box signals
	#     val = value # could also refer directly to control instead of value: self.ui.spinBox.value()
	#
	#     print("Changed to " + str(val))

	@pyqtSlot()
	def on_action_New_triggered(self):  # how to activate menu and toolbar actions!!!

		options = QtWidgets.QFileDialog.Options()
		global DBFilename
		# options |= QtWidgets.QFileDialog.DontUseNativeDialog
		DBFilename, _ = QtWidgets.QFileDialog.getSaveFileName(self,
		                                                      "New Database",
		                                                      "New database",
		                                                      "VGenes database Files (*.vdb);;All Files (*)",
		                                                      options=options)

		if DBFilename != None and DBFilename != '':
			(dirname, filename) = os.path.split(DBFilename)
			(shortname, extension) = os.path.splitext(filename)

			if extension != '.vdb':
				DBFilename = shortname + '.vdb'

			VGenesSQL.creatnewDB(DBFilename)

			self.UpdateRecentList(DBFilename, True)
			question = 'Would you like to enter sequences into your new database?'
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)
			if answer == 'Yes':
				self.ImportOptions.show()

			if os.path.isfile(DBFilename):
				self.LoadDB(DBFilename)

		else:
			self.hide()
			self.ApplicationStarted()

			# self.EditableSqlModel.refresh()

	def UpdateRecentList(self, DBFilename, AddOne):
		# todo need to make this filename fall in VGenes directory upon deployment
		# todo may need to switch this to configparser which is python modle to save ini files
		# todo change to app folder
		try:
			filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'RecentPaths.vtx')
			with open(filename, 'r') as currentfile:
				vv = currentfile

		except:
			filename = '/Volumes/Promise Pegasus/Dropbox/VGenes/RecentPaths.vtx'

		newfile = ''
		linemove = 0

		if AddOne == True:
			exists = False
			i = 0
			if os.path.isfile(DBFilename):
				with open(filename, 'r') as currentfile:
					for line in currentfile:

						line = line.replace('\n', '')
						if line == DBFilename:
							exists = True
							linemove = i

						i += 1
				i = 0

				if exists == False:
					DBFilename += '\n'
					newfile = DBFilename
					with open(filename, 'r') as currentfile:
						for line in currentfile:
							if i < 14:
								newfile += line
							i += 1

				elif exists == True:  # then move that file name to the top of the recents list, but no more then 15
					DBFilename += '\n'
					newfile = DBFilename

					with open(filename, 'r') as currentfile:
						for line in currentfile:
							if i != linemove and i < 14:
								newfile += line
							i += 1
				with open(filename, 'w') as currentfile:
					currentfile.write(newfile)
		elif AddOne == False:
			exists = False

			with open(filename, 'r') as currentfile:
				for line in currentfile:

					# line = line.replace('\n', '')
					if line == DBFilename:
						exists = True

			if exists == True:  # then move that file name to the top of the recents list, but no more then 15
				# DBFilename += '\n'
				newfile = ''

				with open(filename, 'r') as currentfile:
					for line in currentfile:
						if DBFilename != line:
							newfile += line

			with open(filename, 'w') as currentfile:
				currentfile.write(newfile)

	def DecoratePeptide(self):
		Decoration = self.ui.cboDecorate.currentText()

		AASeq = self.ui.txtAASeq.toPlainText()
		GDNAseq = data[80]

		if data[98] == 'Insertion' or data[98] == 'Both':
			mutate = data[97]
			mutations = mutate.split(',')
			for item in mutations:
				if item[:9] == 'Insertion':
					Ievent = item
					Iparts = Ievent.split('-')
					AddAt = int(Iparts[1])
					SeqToAdd = Iparts[2]
					GDNAseq = GDNAseq[:AddAt] + SeqToAdd + GDNAseq[AddAt:]

		GAASeq, ErMessage = VGenesSeq.Translator(GDNAseq, 0)

		if len(AASeq) > len(GAASeq):
			LenTo = len(GAASeq)
			AASeq = AASeq[:LenTo]
		else:
			LenTo = len(AASeq)

		for i in range(0, LenTo - 1):  # first replace bad codons with germline codons
			if AASeq[i] != GAASeq[i]:
				if AASeq[i] == '.' or AASeq[i] == '~' or AASeq[i] == '*':
					AASeq = AASeq[:i] + GAASeq[i] + AASeq[i + 1:]
			elif AASeq[i] == GAASeq[i]:
				if AASeq[i] == '.' or AASeq[i] == '~':
					AASeq = AASeq[:i] + AASeq[i + 1:] + '.'
					GAASeq = GAASeq[:i] + GAASeq[i + 1:] + '.'

		AASeq = AASeq.replace('~', '').replace('.', '')

		ColorMap = []
		WindowSize = self.ui.spinBox.value()
		cursor = self.ui.txtAASeq.textCursor()
		if Decoration == 'None':

			# Setup the desired format for matches
			format = QTextCharFormat()
			format.setForeground(QBrush(QColor("black")))
			format.setBackground(QBrush(QColor("white")))
			format.setFontUnderline(False)
			cursor.setPosition(0)
			cursor.setPosition(len(self.ui.txtAASeq.toPlainText()), QTextCursor.KeepAnchor)
			cursor.mergeCharFormat(format)

			cursor.setPosition(0)

			self.ui.lblScale.setEnabled(False)
			self.ui.lblScaleH.setEnabled(False)
			self.ui.lblScaleL.setEnabled(False)

		elif Decoration == 'Hydrophobicity':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize == 0: self.ui.spinBox.setValue(5)
			WindowSize = self.ui.spinBox.value()
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'Hydrophobicity', WindowSize, False)
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return
			Scale = (-4.5, 4.5)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)


		elif Decoration == 'Hydrophilicity':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize == 0: self.ui.spinBox.setValue(5)
			WindowSize = self.ui.spinBox.value()
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'Hydrophilicity', WindowSize, False)
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return
			Scale = (-3.4, 3.0)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)

		elif Decoration == 'Flexibility':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize == 0: self.ui.spinBox.setValue(9)
			WindowSize = 9
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'Flexibility', WindowSize, False)
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return
			Scale = (0.904, 1.102)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)

		elif Decoration == 'Surface':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize == 0: self.ui.spinBox.setValue(5)
			WindowSize = self.ui.spinBox.value()
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'Surface', WindowSize, False)
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return
			Scale = (0.394, 1.545)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)

		elif Decoration == 'Isoelectric Point (pI)':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize < 5: self.ui.spinBox.setValue(5)
			WindowSize = self.ui.spinBox.value()
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'MapAApI', WindowSize, False)
			if ColorMap == 0: return
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return
			Scale = (0, 14)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)

		elif Decoration == 'Instability':
			self.ui.lblScale.setEnabled(True)
			self.ui.lblScaleH.setEnabled(True)
			self.ui.lblScaleL.setEnabled(True)
			if WindowSize < 8: self.ui.spinBox.setValue(8)
			WindowSize = self.ui.spinBox.value()
			CurPos = (WindowSize // 2)
			ColorMap = VGenesSeq.OtherParam(AASeq, 'MapInstability', WindowSize, False)
			if len(ColorMap) != len(AASeq) - WindowSize + 1:
				sys.stderr.write('Sequence has errors and could not be decorated')
				return

			# for this need to scale relatively but so that anything>40 is in the red  as 40+ = unstable
			Highest = max(ColorMap)
			Lowest = min(ColorMap)
			maxi = ((40 - Lowest) / 8) * 11
			Scale = (Lowest, maxi)  # based on tests paramators
			self.DecorateText(ColorMap, Scale, CurPos, cursor)
			# Analysis = VGenesSeq.OtherParam(AASeq, Param, WindowSize)

	def DecorateText(self, ColorMap, Scale, CurPos, cursor):
		# o in colormap is black text on white background
		#  cursor is cursor from textbox being decorated, i.e.:
		# cursor = self.ui.txtAASeq.textCursor()   when from sequence panel
		#  need provide cursor strat as well...so starts color mid window:
		#          CurPos = (WindowSize // 2)      when from sequence panel
		#  CurPos and cursor will allow me to run through entire text of
		# any window with different paramaters and colormaps

		maxi = Scale[1]
		mini = Scale[0]
		increment = (maxi - mini) / 11

		# Setup the desired format for matches
		format = QTextCharFormat()

		for valueIs in ColorMap:
			if valueIs == 0:
				format.setBackground(QBrush(QColor("white")))
				format.setForeground(QBrush(QColor("black")))
			elif valueIs > mini + (increment * 10):
				format.setBackground(QBrush(QColor("darkred")))
				format.setForeground(QBrush(QColor("white")))
			elif valueIs > mini + (increment * 9) and valueIs <= (mini + (increment * 10)):
				format.setBackground(QBrush(QColor("darkMagenta")))
				format.setForeground(QBrush(QColor("white")))
			elif valueIs > mini + (increment * 8) and valueIs <= (mini + (increment * 9)):
				format.setBackground(QBrush(QColor("red")))
				format.setForeground(QBrush(QColor("black")))
			elif valueIs > mini + (increment * 7) and valueIs <= (mini + (increment * 8)):
				format.setBackground(QBrush(QColor("Magenta")))
				format.setForeground(QBrush(QColor("black")))
			elif valueIs > mini + (increment * 6) and valueIs <= (mini + (increment * 7)):
				format.setBackground(QBrush(QColor("yellow")))
				format.setForeground(QBrush(QColor("black")))
			elif valueIs > mini + (increment * 5) and valueIs <= (mini + (increment * 6)):
				format.setBackground(QBrush(QColor("white")))
				format.setForeground(QBrush(QColor("black")))
			elif valueIs > mini + (increment * 4) and valueIs <= (mini + (increment * 5)):
				format.setBackground(QBrush(QColor("green")))
				format.setForeground(QBrush(QColor("white")))
			elif valueIs > mini + (increment * 3) and valueIs <= (mini + (increment * 4)):
				format.setBackground(QBrush(QColor("darkGreen")))
				format.setForeground(QBrush(QColor("lightGray")))
			elif valueIs > mini + (increment * 2) and valueIs <= (mini + (increment * 3)):
				format.setBackground(QBrush(QColor("blue")))
				format.setForeground(QBrush(QColor("white")))
			elif valueIs > mini + (increment) and valueIs <= (mini + (increment * 2)):
				format.setBackground(QBrush(QColor("darkBlue")))
				format.setForeground(QBrush(QColor("white")))
			elif valueIs <= mini + increment:
				format.setBackground(QBrush(QColor("black")))
				format.setForeground(QBrush(QColor("white")))

			cursor.setPosition(CurPos)
			cursor.setPosition(CurPos + 1, QTextCursor.KeepAnchor)
			cursor.mergeCharFormat(format)

			CurPos += 1


			# Move to the next match

	@pyqtSlot(int)
	def on_spinBox_valueChanged(self):
		if self.ui.cboDecorate.currentText() != "None":
			self.DecoratePeptide()

	@pyqtSlot()
	def on_action_actionCopy_triggered(self):
		# shortcut=QKeySequence.Cut,
		# statusTip="Cut the current selection's contents to the clipboard",
		if self.ui.tabWidget.currentIndex() == 2:
			self.ui.txtDNASeq.copy()
		elif self.ui.tabWidget.currentIndex() == 1:
			self.ui.txtSeqAlignment.copy()

	@pyqtSlot()
	def on_action_Open_triggered(self):  # how to activate menu and toolbar actions!!!
		self.GOOpen(True)

	def GOOpen(self, GetName):
		global DBFilename
		if GetName == True:
			DBFilename = ''
			typeOpen = 'db'
			DBFilename = openFile(self, typeOpen)
			(dirname, filename) = os.path.split(DBFilename)

		if DBFilename != None:
			if os.path.isfile(DBFilename):
				self.LoadDB(DBFilename)
			else:
				VGenesSQL.creatnewDB(DBFilename)

			self.UpdateRecentList(DBFilename, True)
		else:
			self.hide()
			self.ApplicationStarted()

		if GetName == True:
			self.SaveBackup

			# self.LoadDB(DBFilename)
			# self.EditableSqlModel.refresh()

	@pyqtSlot()
	def on_action_Import_triggered(self):
		self.ImportOptions.show()

	@pyqtSlot()
	def on_actionGenerate_from_file_triggered(self):  # how to activate menu and toolbar actions!!!


		fileNames = openFiles()
		self.ProcessSeqFiles(fileNames)

	def ShowVGenesText(self, filename):

		self.TextEdit.show()
		if filename != '':
			self.TextEdit.loadFile(filename)

	def on_comboBoxSpecies_editTextChanged(self):  # had to connect above as couldn't figure signal
		# species = self.ui.comboBoxSpecies.currentText()
		print('clicked')

	@pyqtSlot(int)
	def on_dial_valueChanged(self, value):  # how to handle dial signals
		FieldCheck = self.FieldChangeCheck()
		if FieldCheck == 'exit':
			return

		val = value  # could also refer directly to control instead of value: self.ui.spinBox.value()
		self.ui.horizontalScrollBar.setValue(self.ui.dial.value())
		# self.ui.txtName.setText("Changed to " + str(val))

	@pyqtSlot(int)
	def on_horizontalScrollBar_valueChanged(self, value):  # how to handle slider and lcdNumber signals
		val = value  # could also refer directly to control instead of value: self.ui.spinBox.value()
		self.ui.dial.setValue(self.ui.horizontalScrollBar.value())
		self.ui.lcdNumber_current.display(val)
		self.DialScroll(value)

	@pyqtSlot(int)
	def DialScroll(self, value):
		currentRow = self.ui.tableView.currentIndex().row()

		if currentRow != value:
			currentColumn = self.ui.tableView.currentIndex().column()
			if currentColumn == -1: currentColumn = 0
			model = self.ui.tableView.model()
			records = model.rowCount()
			if value < records and value > -1:
				index = model.index(value, currentColumn)
				self.ui.tableView.setCurrentIndex(index)

				self.updateF(value)

	def findTableViewRecord(self, FieldName):

		# # todo way to find records in TableView but not connected
		model = self.ui.tableView.model()
		proxy = QSortFilterProxyModel
		proxy.setSourceModel(model)
		proxy.setFilterKeyColumn(0)
		proxy.setFilterFixedString(FieldName)
		MatchingIndex = proxy.mapToSource(proxy.index(0, 0))  # (QSortFilterProxyModel,0,0))
		# self.ui.tableView.scrollTo(MatchingIndex)
		self.ui.tableView.setCurrentIndex(MatchingIndex)

		# todo old code to update fields after tree clicked using Name Index now





		try:
			self.SeqButton(LastPushed)
		except:
			self.SeqButton('v')
		JustMoved = False

	@pyqtSlot()
	def on_actionMove_Up_triggered(self):  # how to activate menu and toolbar actions!!!
		self.MoveRecord('up')

	@pyqtSlot()
	def on_actionIncrease_font_size_triggered(self):

		if self.ui.tabWidget.currentIndex() == 3:
			val = int(self.ui.spnAlignFont.value())
			val += 1
			self.ui.spnAlignFont.setValue(val)
		elif self.ui.tabWidget.currentIndex() == 2:
			FontIs = self.ui.txtDNASeq.currentFont()
			font = QFont(FontIs)

			FontSize = int(font.pointSize())
			if FontSize < 36:
				FontSize += 1
			font.setPointSize(FontSize)
			font.setFamily('Lucida Grande')

			self.ui.txtDNASeq.setFont(font)
			self.ui.txtAASeq.setFont(font)

		elif self.ui.tabWidget.currentIndex() == 1:
			FontIs = self.ui.tableView.font()
			font = QFont(FontIs)

			FontSize = int(font.pointSize())
			if FontSize > 7:
				FontSize += 1
			font.setPointSize(FontSize)
			font.setFamily('Lucida Grande')

			self.ui.tableView.setFont(font)
			# self.ui.tableView.resizeColumnsToContents()

	@pyqtSlot()
	def on_actionDecrease_font_size_triggered(self):

		if self.ui.tabWidget.currentIndex() == 3:
			val = int(self.ui.spnAlignFont.value())
			val -= 1
			self.ui.spnAlignFont.setValue(val)
		elif self.ui.tabWidget.currentIndex() == 2:
			FontIs = self.ui.txtDNASeq.currentFont()
			font = QFont(FontIs)

			FontSize = int(font.pointSize())
			if FontSize > 7:
				FontSize -= 1
			font.setPointSize(FontSize)
			font.setFamily('Lucida Grande')

			self.ui.txtDNASeq.setFont(font)
			self.ui.txtAASeq.setFont(font)
		elif self.ui.tabWidget.currentIndex() == 1:
			FontIs = self.ui.tableView.font()
			font = QFont(FontIs)

			FontSize = int(font.pointSize())
			if FontSize > 7:
				FontSize -= 1
			font.setPointSize(FontSize)
			font.setFamily('Lucida Grande')

			self.ui.tableView.setFont(font)

			# self.ui.tableView.resizeColumnsToContents()

	def on_spnAlignFont_valueChanged(self, value):
		self.AlignFont()

	def AlignFont(self):
		FontSize = int(self.ui.spnAlignFont.text())
		font = QFont()
		font.setFamily("Courier New")
		font.setPointSize(FontSize)

		self.ui.txtSeqAlignment.setFont(font)

	@pyqtSlot()
	def on_actionMove_Down_triggered(self):  # how to activate menu and toolbar actions!!!
		self.MoveRecord('down')

	@pyqtSlot()
	def on_actionTop_triggered(self):  # how to activate menu and toolbar actions!!!
		self.MoveRecord('top')

	@pyqtSlot()
	def on_actionBottom_triggered(self):  # how to activate menu and toolbar actions!!!
		self.MoveRecord('bottom')

	def TransLateFieldtoReal(self, FieldName, ToField):
		i = 0
		if ToField == True:
			for field in RealNameList:
				if FieldName == field:
					return (FieldList[i])
				i += 1
		elif ToField == False:
			for field in FieldList:
				if FieldName == field:
					return (RealNameList[i])
				i += 1

	def FieldChangeCheck(self):

		if FieldChanged == True:
			# msg = 'Changes were made to this record, save them (if no the changes will be lost)?'
			# buttons = 'YN'
			# answer = questionMessage(self, msg, buttons)
			answer = 'Yes'
			if answer == 'Yes':
				self.on_action_Save_triggered()
				return "exit"
			else:
				FieldsChanged.clear

	def MoveRecord(self, direction):

		FieldCheck = self.FieldChangeCheck()
		if FieldCheck == 'exit':
			return

		currentRow = self.ui.tableView.currentIndex().row()
		currentColumn = self.ui.tableView.currentIndex().column()

		if currentRow == -1:
			currentRow = 0
		if currentColumn == -1:
			currentColumn = 0
		model = self.ui.tableView.model()
		records = model.rowCount()

		if direction == 'up':
			if currentRow > 0:
				currentRow -= 1

		elif direction == 'down':
			if currentRow < records - 1:
				currentRow += 1

			else:
				currentRow = records
		elif direction == 'top':
			currentRow = 0

		elif direction == 'bottom':
			currentRow = records - 1
		global JustMoved
		JustMoved = True
		index = model.index(currentRow, currentColumn)
		self.ui.tableView.setCurrentIndex(index)
		self.ui.radioButtonSeqView.setChecked(True)

		try:
			self.SeqButton(LastPushed)
		except:
			self.SeqButton('v')
		JustMoved = False


		# self.findTreeItem(data[0])

	def createView(self, model):
		# def createView(title, model):
		global views
		views = []
		view = self.ui.tableView
		# self.ui.tableView.keyPressEvent()
		views.append(view)
		view.setModel(model)

	def LoadDB(self, DBFilename):

		if not self.createConnection(DBFilename):
			sys.exit(1)

		editableModel = EditableSqlModel()
		self.initializeModel(editableModel)
		self.createView(editableModel)

		self.OnOpen()
		titletext = 'VGenes - ' + DBFilename
		self.setWindowTitle(titletext)

	def GenerateNameIndex(self):
		model = self.ui.tableView.model()
		index = model.index(0, 0)
		self.ui.tableView.setCurrentIndex(index)
		Maxi = model.rowCount()

		for i in range(0, Maxi):
			index = model.index(i, 0)
			NameIs = str(model.data(index))
			NameIndex[NameIs] = i

	def OnOpen(self):

		# self.ui.tableView.ColumnsToContents()
		self.ui.tableView.resizeColumnsToContents()

		self.updateF(-2)
		self.ui.lcdNumber_current.display(1)
		# self.ui.txtGotoRecord.setPlainText('1')
		model = self.ui.tableView.model()

		index = model.index(0, 0)
		self.ui.tableView.setCurrentIndex(index)
		self.GenerateNameIndex()

		self.ui.cboFindField.clear()
		self.ui.cboTreeOp1.clear()
		self.ui.cboTreeOp2.clear()
		self.ui.cboTreeOp3.clear()

		self.ui.cboFindField.addItem("Project")
		self.ui.cboFindField.addItem("Grouping")
		self.ui.cboFindField.addItem("Subgroup")

		for item in RealNameList:
			self.ui.cboFindField.addItem(item)
			self.ui.cboTreeOp1.addItem(item)
			self.ui.cboTreeOp2.addItem(item)
			self.ui.cboTreeOp3.addItem(item)
		self.ui.cboTreeOp1.addItem('None')
		self.ui.cboTreeOp2.addItem('None')
		self.ui.cboTreeOp3.addItem('None')
		self.ui.cboTreeOp1.setCurrentText('Project')
		self.ui.cboTreeOp2.setCurrentText('Grouping')
		self.ui.cboTreeOp3.setCurrentText('Subgroup')

		SQLFields = ('Project', 'Grouping', 'SubGroup')
		self.initializeTreeView(SQLFields)

		self.findTreeItem(data[0])

	def createConnection(self, DBFilename):

		db = QSqlDatabase.addDatabase('QSQLITE')
		# db.setDatabaseName(':memory:')
		db.setDatabaseName(DBFilename)
		if not db.open():
			QMessageBox.critical(None, "Cannot open database",
			                     "Unable to establish a database connection.\n",
			                     QMessageBox.Cancel)
			return False

		return True

	def initializeModel(self, model):

		model.setQuery('select * from vgenesdb ORDER BY Project, Grouping, SubGroup, SeqName')
		while model.canFetchMore():
			model.fetchMore()
		NumRows = model.rowCount()

		i = 0
		for item in FieldList:
			model.setHeaderData(i, Qt.Horizontal, item)

		model.setHeaderData(0, Qt.Horizontal, "SeqName")

	@pyqtSlot()
	def on_radioButtonGermView_clicked(self):
		self.ui.txtDNASeq.setText(data[80])
		try:
			self.SeqButton(LastPushed)
		except:
			self.SeqButton('v')

	def on_txtVgene_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVgene.toPlainText()
		LastSelected = ('V1', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtDgene_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDgene.toPlainText()
		LastSelected = ('D1', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtJgene_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJgene.toPlainText()
		LastSelected = ('J1', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	# todo need to put LastSelected for rest of fields

	def on_txtVLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVLocus.toPlainText()
		LastSelected = ('VLocus', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtDLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDLocus.toPlainText()
		LastSelected = ('DLocus', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtJLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJLocus.toPlainText()
		LastSelected = ('JLocus', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	@pyqtSlot()
	def on_txtProject_textChanged(self):

		valueTo = self.ui.txtProject.toPlainText()
		FieldIS = 'Project'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtProject_2_textChanged(self):

		valueTo = self.ui.txtProject_2.toPlainText()
		FieldIS = 'Project'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtGroup_textChanged(self):

		valueTo = self.ui.txtGroup.toPlainText()
		FieldIS = 'Grouping'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtGroup_2_textChanged(self):

		valueTo = self.ui.txtGroup_2.toPlainText()
		FieldIS = 'Grouping'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtSubGroup_textChanged(self):

		valueTo = self.ui.txtSubGroup.toPlainText()
		FieldIS = 'SubGroup'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtSubGroup_2_textChanged(self):

		valueTo = self.ui.txtSubGroup_2.toPlainText()
		FieldIS = 'SubGroup'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtName_textChanged(self):

		valueTo = self.ui.txtName.toPlainText()
		FieldIS = 'SeqName'
		# self.NameChange(valueTo)
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtName_2_textChanged(self):

		valueTo = self.ui.txtName_2.toPlainText()
		FieldIS = 'SeqName'

		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def NameChange(self, ToName):
		item = self.ui.treeWidget.currentItem()
		# global OldName
		OldName = str(item.text(0))
		TreeIndex = NameIndex[OldName]
		del NameIndex[OldName]
		NameIndex[ToName] = TreeIndex

		item.setText(0, ToName)
		# item.setData(0,)

	@pyqtSlot()
	def on_txtReadingFrame_textChanged(self):

		valueTo = self.ui.txtReadingFrame.toPlainText()
		FieldIS = 'ReadingFrame'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtID_textChanged(self):

		valueTo = self.ui.txtID.toPlainText()
		FieldIS = 'IDEvent'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtStop_textChanged(self):

		valueTo = self.ui.txtStop.toPlainText()
		FieldIS = 'StopCodon'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtProductive_textChanged(self):

		valueTo = self.ui.txtProductive.toPlainText()
		FieldIS = 'productive'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtVgene_textChanged(self):

		valueTo = self.ui.txtVgene.toPlainText()
		FieldIS = 'V1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtVgeneSeq_textChanged(self):

		valueTo = self.ui.txtVGeneSeq.toPlainText()
		FieldIS = 'V1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtDgene_textChanged(self):

		valueTo = self.ui.txtDgene.toPlainText()
		FieldIS = 'D1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtDgeneSeq_textChanged(self):

		valueTo = self.ui.txtDGeneSeq.toPlainText()
		FieldIS = 'D1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtJgene_textChanged(self):

		valueTo = self.ui.txtJgene.toPlainText()
		FieldIS = 'J1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtJgeneSeq_textChanged(self):

		valueTo = self.ui.txtJGeneSeq.toPlainText()
		FieldIS = 'J1'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtVLocus_textChanged(self):

		valueTo = self.ui.txtVLocus.toPlainText()
		FieldIS = 'VLocus'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtDLocus_textChanged(self):

		valueTo = self.ui.txtDLocus.toPlainText()
		FieldIS = 'DLocus'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtJLocus_textChanged(self):

		valueTo = self.ui.txtJLocus.toPlainText()
		FieldIS = 'JLocus'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtQuality_textChanged(self):

		valueTo = self.ui.txtQuality.toPlainText()
		FieldIS = 'Quality'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtQuality_2_textChanged(self):

		valueTo = self.ui.txtQuality_2.toPlainText()
		FieldIS = 'Quality'
		self.FieldChanger(valueTo, FieldIS)

	def FieldChanger(self, valueTo, FieldIS):
		global FieldChanged
		if JustMovedIt == False:
			FieldChanged = True
			self.UpdateFChanges(data[119], valueTo, FieldIS)
			if FieldIS == 'SeqName':
				self.NameChange(valueTo)

	def UpdateFChanges(self, IDIs, valueTo, FieldIs):
		IsThere = False
		i = 0
		ListItem = ()
		for item in FieldsChanged:
			if item[0] == IDIs and item[2] == FieldIs:
				if item[1] != valueTo:
					ListItem = (IDIs, valueTo, FieldIs)
					FieldsChanged[i] = ListItem
				IsThere = True
			i += 1

		if IsThere == False:
			ListItem = (IDIs, valueTo, FieldIs)
			FieldsChanged.append(ListItem)

	def on_action_Save_triggered(self):
		global FieldChanged
		F1 = self.ui.cboTreeOp1.currentText()
		FR1 = self.TransLateFieldtoReal(F1, True)
		F2 = self.ui.cboTreeOp2.currentText()
		FR2 = self.TransLateFieldtoReal(F2, True)
		F3 = self.ui.cboTreeOp3.currentText()
		FR3 = self.TransLateFieldtoReal(F3, True)
		NeedTree = False
		NeedOpen = False
		if len(FieldsChanged) > 0:
			model = self.ui.tableView.model()
			for item in FieldsChanged:
				ID = item[0]
				ItemValue = item[1]
				FieldName = item[2]
				if FieldName == FR1 or FieldName == FR2 or FieldName == FR3:
					NeedTree = True
				if FieldName == 'SeqName':
					NeedTree = False



					# NeedOpen = True
				VGenesSQL.UpdateField(ID, ItemValue, FieldName, DBFilename)
				model.refresh()

			FieldsChanged.clear()

			FieldChanged = False
		# self.GOOpen(False)
		if NeedTree == True:
			self.on_btnUpdateTree_clicked()
			# if NeedOpen == True:
			#     self.GOOpen(False)

	def on_txtProject_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtProject.toPlainText()
		LastSelected = ('Project', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtGroup_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtGroup.toPlainText()
		LastSelected = ('Grouping', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtSubGroup_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtSubGroup.toPlainText()
		LastSelected = ('SubGroup', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtVend_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVend.toPlainText()
		LastSelected = ('VSeqend', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtD_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtD.toPlainText()
		LastSelected = ('Dregion', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtVD_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVD.toPlainText()
		LastSelected = ('VDJunction', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtDJ_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDJ.toPlainText()
		LastSelected = ('DJJunction', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtJend_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJend.toPlainText()
		LastSelected = ('begJ', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtCDR3DNA_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3DNA.toPlainText()
		LastSelected = ('CDR3DNA', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtCDR3AA_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3AA.toPlainText()
		LastSelected = ('CDR3AA', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtCDR3Length_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3Length.toPlainText()
		LastSelected = ('CDR3Length', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtCDR3pI_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3pI.toPlainText()
		LastSelected = ('CDR3pI', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtCDR3MW_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3MW.toPlainText()
		LastSelected = ('CDR3MW', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtProductive_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtProductive.toPlainText()
		LastSelected = ('productive', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtReadingFrame_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtReadingFrame.toPlainText()
		LastSelected = ('ReadingFrame', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtStop_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtStop.toPlainText()
		LastSelected = ('StopCodon', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtQuality_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtQuality.toPlainText()
		LastSelected = ('Quality', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtQuality_2_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtQuality_2.toPlainText()
		LastSelected = ('Quality', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtID_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtID.toPlainText()
		LastSelected = ('IDEvent', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtClonalPool_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtClonalPool.toPlainText()
		LastSelected = ('ClonalPool', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	@pyqtSlot()
	def on_txtClonalPool_textChanged(self):
		valueTo = self.ui.txtClonalPool.toPlainText()
		FieldIS = 'ClonalPool'
		self.FieldChanger(valueTo, FieldIS)
		# if data[88] != 'None':
		# 	try:
		# 		if type(int(self.ui.txtClonalPool.toPlainText())) != int:
		# 			msg = 'Please enter integers only for clonal pool number.'
		# 			QtWidgets.QMessageBox.critical(None, "Integers only", msg, QtWidgets.QMessageBox.Cancel)
		# 			self.ui.txtClonalPool.setPlainText('0')
		# 	except:
		#
		# 		msg = 'Please enter integers only for clonal pool number.'
		# 		QtWidgets.QMessageBox.critical(None, "Integers only", msg, QtWidgets.QMessageBox.Cancel)
		# 		self.ui.txtClonalPool.setPlainText('0')

				# LastSelected = ('ClonalPool', valueis)

	def on_listViewSpecificity_selectionChanged(self):
		global LastSelected
		valueis = self.ui.listViewSpecificity.value()
		LastSelected = ('Specificity', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_listViewSpecificity_2_selectionChanged(self):
		global LastSelected
		valueis = self.ui.listViewSpecificity_2.value()
		LastSelected = ('Subspecificity', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	@pyqtSlot()
	def on_btnClearTreeChecks_clicked(self):
		self.clearTreeChecks()

	@pyqtSlot()
	def on_actionSuggestCanonical_triggered(self):





		self.ui.treeWidget.expandAll()

		fields = self.ui.cboTreeOp1.currentText()
		field1 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field1 == item:
				field1Value = data[i]
			i += 1

		fields = self.ui.cboTreeOp2.currentText()
		field2 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field2 == item:
				field2Value = data[i]
			i += 1

		fields = self.ui.cboTreeOp3.currentText()
		field3 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field3 == item:
				field3Value = data[i]
			i += 1

		if field1 == '': field1 = 'None'
		if field2 == '': field1 = 'None'
		if field3 == '': field1 = 'None'


		fields = ['SeqName', 'Sequence', 'VLocus', 'JLocus', 'Vbeg', 'Jend', 'SubGroup']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		# DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		# SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '"' # AND ' + Field3 + ' = "' + Vcolumn3 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
		foundRecs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		NumFound = len(foundRecs)
		# foundRecs.sort(key=itemgetter(0))
		i = 0
		Genes = []
		gene = []
		ToCLustalgene = []
		ToClustal = []
		LastSeq = 'FirstOne'  # ('FirstOne', 'x', 'x', 'x', 'x', 'x')
		ErrMes = ''
		self.clearTreeChecks()
		foundRecs.sort(key=itemgetter(0))
		NumRecords = len(foundRecs)
		NumDone = 0

		for item in foundRecs:
			NumDone += 1
			gene.clear()
			for i in range(0, 6):
				gene.append(item[i])
			SeqName = item[0]
			if NumDone == 1:
				currentitemIs = SeqName

			# Sequence = item[1]
			# Vgene = item[2]
			# Jgene = item[3]
			# Vbeg = item[4]
			# Jend = item[5]
			# if SeqName  == 'A116_1F02H-3':
			#     print('stop')
			if SeqName[len(SeqName) - 2] == '-':
				NameComp = SeqName[:len(SeqName) - 2]
			elif SeqName[len(SeqName) - 3] == '-':
				NameComp = SeqName[:len(SeqName) - 3]
			else:
				answer = informationMessage(self,
				                            'Make sure individual sequence sets are named the same\nwith only "-#" at the end as in: "045-2B06-1, 045-2B05-2"',
				                            'OK')
				return

			if LastSeq != 'FirstOne':
				# if LastSeq[len(LastSeq)-2] == '-':
				#     LastSeqC = LastSeq[:len(LastSeq)-2]
				# elif LastSeq[len(LastSeq)-3] == '-':
				#     LastSeqC = LastSeq[:len(LastSeq)-3]


				if NameComp == LastSeq:
					Genes.append(tuple(gene))

				else:
					#         code to compare sequences and then clear genes and start with current one
					if len(Genes) > 1:
						longest = 0
						Nucs = {'As': 0, 'Gs': 0, 'Cs': 0, 'Ts': 0, 'Ns': 0}
						LostCon = False
						Consensus = ''
						for i in range(0, len(Genes)):  # first get start and end
							SeqName = Genes[i][0]
							Sequence = Genes[i][1]
							Sequence = Sequence[(int(Genes[i][4])):(int(Genes[i][5]))]
							Sequence = Sequence.upper()
							# if len(Sequence)>longest: longest = len(Sequence)
							ToCLustalgene.append(SeqName)
							ToCLustalgene.append(Sequence)
							ToClustal.append(tuple(ToCLustalgene))
							ToCLustalgene.clear()
						ClustalOut = VGenesSeq.ClustalO(ToClustal, 1000, False)
						ToClustal.clear()
						outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega',
						                           'my-out-seqs.fa')
						Aligned = VGenesSeq.readClustalOutput(outfilename)
						Genes.clear()
						GeneDict = {}
						for i in range(0, len(Aligned)):
							SeqName = Aligned[i][0]
							Seq = Aligned[i][1]
							longest = len(Seq)
							GeneDict.update({SeqName: 0})
						# x.update({3:4})
						for j in range(0, longest):  # build consensus
							Nucs = {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'N': 0}
							for i in range(0, len(Aligned)):
								try:
									nuc = Aligned[i][1][j]
								except:
									print('stop')
								if nuc == 'A' or nuc == 'a':
									Nucs['A'] += 1
								elif nuc == 'G' or nuc == 'g':
									Nucs['G'] += 1
								elif nuc == 'C' or nuc == 'c':
									Nucs['C'] += 1
								elif nuc == 'T' or nuc == 't':
									Nucs['T'] += 1
								else:
									Nucs['N'] += 1
							Cnuc = max(Nucs, key=Nucs.get)
							Cfreq = Nucs[Cnuc]
							PerC = Cfreq / len(Aligned)
							if Cfreq / len(Aligned) >= 0.5:
								Consensus += Cnuc
							else:
								if LastSeq == 'LastOne':
									print('stop')

								ErrMes += LastSeq + ' has no consensus sequence.\n'
								LostCon = True
								Consensus += 'X'

						for j in range(0, longest):

							for i in range(0, len(Aligned)):
								SeqName = Aligned[i][0]
								try:
									nuc = Aligned[i][1][j]
								except:
									print('stop')
								ConNuc = Consensus[j]
								# SetVal = int(Genes[i][1])
								if nuc != ConNuc:
									GeneDict[SeqName] += 1
									# SetVal = int(Genes[i][1])
									# SetVal += 1
									# Genes[i][1] = SetVal

						NumRecs = len(Aligned)
						ConFound = False
						for i in range(0, NumRecs):
							SeqName = Aligned[i][0]
							# if SeqName  == 'A116_1E01H-1':
							#     print('it')
							Mutations = GeneDict[SeqName]
							if Mutations == 0:
								found = self.ui.treeWidget.findItems(SeqName, Qt.MatchRecursive, 0)
								ConFound = True
								for record in found:
									record.setCheckState(0, Qt.Checked)
								break
						if ConFound == False:
							if LastSeq == 'LastOne':
								print('stop')
							ErrMes += LastSeq + ' has no consensus sequence.\n'







					else:
						if LastSeq == 'LastOne':
							print('stop')
						ErrMes += LastSeq + ' has no consensus sequence.\n'
						LostCon = True

					Genes.clear()
					Genes.append(tuple(gene))

			else:
				Genes.append(tuple(gene))

			LastSeq = NameComp
			if NumDone == NumRecords - 1:
				LastSeq = "LastOne"


				# def on_btnFieldBulk_clicked(self):
				#     search = self.ui.txtFieldSearch.toPlainText()
				#     field = self.ui.cboFindField.currentText()

		QueryIS = 'Would you like to move the sequences selected for expression to a new Subgroup?'
		buttons = 'YN'
		answer = questionMessage(self, QueryIS, buttons)

		if answer == 'Yes':
			QueryIS = 'Enter text to be concatenated to the sub group name'
			DefaultText = '-Expressed'  # data[77] + '-Expressed'
			EditSubgroup = setText(self, QueryIS, DefaultText)

		elif answer == 'No':
			EditSubgroup = 'Cancelled Action'

		if EditSubgroup != 'Cancelled Action':
			fields = ['SeqName', 'SubGroup']
			SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
			foundRecs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for item in foundRecs:
				SeqName = item[0]
				Subgroup = item[1]
				NewSub = Subgroup + EditSubgroup

				SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

				WhereStart = SQLStatement.find('WHERE')
				WhereState = SQLStatement[WhereStart - 1:]  # + ' AND '
				SQLStatement = 'UPDATE vgenesDB SET Subgroup = "' + NewSub + '" WHERE SeqName = "' + SeqName + '"'

				foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

			model = self.ui.tableView.model()

			model.refresh()

			self.on_btnUpdateTree_clicked()
			# self.ui.txtFieldSearch.setPlainText(EditSubgroup)
			# self.ui.cboFindField.setCurrentText('Subgroup')
			# Doit = self.on_btnFieldBulk_clicked()
			# self.on_btnUpdateTree_clicked()

		if ErrMes != '':
			Style = 'standard'

			self.ShowVGenesTextEdit(ErrMes, Style)
		self.ui.treeWidget.collapseAll()
		if currentitemIs:
			self.findTreeItem(currentitemIs)

	@pyqtSlot()
	def on_btnCopyRecords_clicked(self):

		New = False
		self.MoveRecords(New)



	@pyqtSlot()
	def on_pushButtonSimilar_clicked(self):
		# print('searched')
		global wasClicked
		wasClicked = False

		value = self.ui.treeWidget.selectedItems()
		currentitemIs = ''

		for item in value:
			currentitemIs = item.text(0)

		self.clearTreeChecks()
		if LastSelected:

			fieldsearch = LastSelected[0]
		else:
			answer = informationMessage(self,
			                            'No field was seleceted.\nClick a field from the Record tab and records with similar values will be checked.',
			                            'OK')
			return

		search = LastSelected[1]

		fields = self.ui.cboTreeOp1.currentText()
		field1 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field1 == item:
				field1Value = data[i]
			i += 1

		fields = self.ui.cboTreeOp2.currentText()
		field2 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field2 == item:
				field2Value = data[i]
			i += 1

		fields = self.ui.cboTreeOp3.currentText()
		field3 = self.TransLateFieldtoReal(fields, True)
		i = 0
		for item in FieldList:
			if field3 == item:
				field3Value = data[i]
			i += 1

		if field1 == '': field1 = 'None'
		if field2 == '': field1 = 'None'
		if field3 == '': field1 = 'None'

		# global RefreshSQL
		if field1 == 'None' or field1 is None:
			SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '"'
		elif field2 == 'None' or field2 is None:
			SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '"'
		elif field3 == 'None' or field3 is None:
			SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '" AND ' + field2 + ' = "' + field2Value + '"'
		else:
			SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '" AND ' + field2 + ' = "' + field2Value + '" AND ' + field3 + ' = "' + field3Value + '"'

		# SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '"' # AND ' + Field3 + ' = "' + Vcolumn3 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
		foundRecs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		NumFound = len(foundRecs)
		i = 0
		for item in foundRecs:
			Seqname = item[0]
			found = self.ui.treeWidget.findItems(Seqname, Qt.MatchRecursive, 0)
			i += 1
			for record in found:
				if i == NumFound - 1:
					global wasClicked
					wasClicked = True
				# global wasClicked
				# wasClicked = True
				record.setCheckState(0, Qt.Checked)
		NewLbl = self.ui.label_Name.text()
		NewLbl += ', ' + str(NumFound) + ' selected'
		self.ui.label_Name.setText(NewLbl)

		self.findTreeItem(currentitemIs)

	@pyqtSlot()
	def on_btnFieldSearch_clicked(self):

		search = self.ui.txtFieldSearch.toPlainText()
		field = self.ui.cboFindField.currentText()
		fieldsearch = self.TransLateFieldtoReal(field, True)

		#     if search != '':
		#         SQLstatement = 'select * from vgenesdb WHERE ' + fieldsearch + ' LIKE ' + search + ' ORDER BY ' + fieldsearch + ', Project, Grouping, SubGroup, SeqName'
		# #         todo select one if exact or if no then all similar items in tree and table...best if tree checkable

		if self.ui.rdoLocal.isChecked():
			fields = self.ui.cboTreeOp1.currentText()
			field1 = self.TransLateFieldtoReal(fields, True)
			i = 0
			for item in FieldList:
				if field1 == item:
					field1Value = data[i]
				i += 1

			fields = self.ui.cboTreeOp2.currentText()
			field2 = self.TransLateFieldtoReal(fields, True)
			i = 0
			for item in FieldList:
				if field2 == item:
					field2Value = data[i]
				i += 1

			fields = self.ui.cboTreeOp3.currentText()
			field3 = self.TransLateFieldtoReal(fields, True)
			i = 0
			for item in FieldList:
				if field3 == item:
					field3Value = data[i]
				i += 1

			if field1 == '': field1 = 'None'
			if field2 == '': field1 = 'None'
			if field3 == '': field1 = 'None'

			# global RefreshSQL
			if field1 == 'None' or field1 is None:
				SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '"'
			elif field2 == 'None' or field2 is None:
				SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '"'
			elif field3 == 'None' or field3 is None:
				SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '" AND ' + field2 + ' = "' + field2Value + '"'
			else:
				SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '" AND ' + field2 + ' = "' + field2Value + '" AND ' + field3 + ' = "' + field3Value + '"'

		else:
			SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '"'

		# SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE ' + fieldsearch + ' = "' + search + '" AND ' + field1 + ' = "' + field1Value + '"' # AND ' + Field3 + ' = "' + Vcolumn3 + '" ORDER BY Project, Grouping, SubGroup, SeqName'
		foundRecs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		NumFound = len(foundRecs)
		answer = "No"
		global wasClicked
		wasClicked = False
		if NumFound > 1:
			question = 'More then one record was found with this search criteria. Check all (Y) or just navigate to the first (N)?'
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)
			if answer == 'Yes':
				self.clearTreeChecks()
		i = 0
		FindName = ''
		for item in foundRecs:
			i += 1
			Seqname = item[0]
			if i == 1:
				FindName = Seqname
			found = self.ui.treeWidget.findItems(Seqname, Qt.MatchRecursive, 0)
			if answer == 'Yes':
				for record in found:
					# global wasClicked
					# wasClicked = True
					if i == NumFound - 1:
						global wasClicked
						wasClicked = True
					record.setCheckState(0, Qt.Checked)

		self.findTreeItem(FindName)
		NewLbl = self.ui.label_Name.text()
		NewLbl += ', ' + str(NumFound) + ' selected'
		self.ui.label_Name.setText(NewLbl)

	@pyqtSlot()
	def on_btnFieldBulk_clicked(self):
		search = self.ui.txtFieldSearch.toPlainText()
		if search == '':
			answer = informationMessage(self, 'Enter a value to the edit/search field', 'OK')
			return

		field = self.ui.cboFindField.currentText()
		fieldsearch = self.TransLateFieldtoReal(field, True)
		if fieldsearch == "SeqName":
			answer = informationMessage(self, 'Sequence names cannot be changed in bulk', 'OK')
			return

		value = self.ui.treeWidget.selectedItems()
		currentitemIs = ''

		for item in value:
			currentitemIs = item.text(0)

		fields = self.ui.cboTreeOp1.currentText()
		field1 = self.TransLateFieldtoReal(fields, True)

		fields = self.ui.cboTreeOp2.currentText()
		field2 = self.TransLateFieldtoReal(fields, True)
		# field2Index = self.ui.cboTreeOp2.currentIndex()

		fields = self.ui.cboTreeOp3.currentText()
		field3 = self.TransLateFieldtoReal(fields, True)
		# field3Index = self.ui.cboTreeOp3.currentIndex()

		if field1 == '': field1 = 'None'
		if field2 == '': field1 = 'None'
		if field3 == '': field1 = 'None'
		TreeFields = []
		TreeFields.append(field1)
		TreeFields.append(field2)
		TreeFields.append(field3)

		fields = ['ID']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		WhereStart = SQLStatement.find('WHERE')
		WhereState = SQLStatement[WhereStart - 1:]  # + ' AND '
		SQLStatement = 'UPDATE vgenesDB SET ' + fieldsearch + ' = "' + search + '"' + WhereState  # -7.000000' WHERE locId = 173567"
		# ' WHERE SeqName = "A116_1B04H-2" OR SeqName = "A116_1B04H-3"'
		foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

		model = self.ui.tableView.model()

		model.refresh()

		# self.GOOpen(False)



		if fieldsearch == field1 or fieldsearch == field2 or fieldsearch == field3:
			self.on_btnUpdateTree_clicked()

		# self.LoadDB(DBFilename)



		self.findTreeItem(currentitemIs)



	@pyqtSlot()
	def ReportOptions(self):
		# self.ui.cboReportOptions.currentIndexChanged()
		option = self.ui.cboReportOptions.currentText()
		if option == 'FASTA Nucleotide file':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)

		if option == 'FASTA Amino Acid file':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)

		if option == 'FASTA Germline Nucleotide file':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)

		if option == 'FASTA Germline Amino Acid file':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)


		elif option == 'AbVec cloning PCR':
			self.ui.ckReportCSV.setDisabled(False)
			self.ui.ckReportCSV.setChecked(True)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(False)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(False)

		elif option == 'Sequence summary':
			self.ui.ckReportCSV.setDisabled(False)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)


		elif option == 'Comma seperated values (.csv)':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(True)
			self.ui.ckReportDisplay.setDisabled(True)
			self.ui.ckReportDisplay.setChecked(False)
			self.ui.ckReportText.setDisabled(True)
			self.ui.ckReportText.setChecked(False)


		elif option == 'Custom report':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(True)
			self.ui.ckReportDisplay.setChecked(False)
			self.ui.ckReportText.setDisabled(True)
			self.ui.ckReportText.setChecked(False)

	@pyqtSlot()
	def on_btnQuery_clicked(self):
		option = ''

		# if self.ui.rdoReport.isChecked():
		option = self.ui.cboReportOptions.currentText()

		VReports.StandardReports(self, option, data[0], DBFilename)



	@pyqtSlot()
	def SaveBackup(self):

		import shutil
		Backfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'BackUP.vdb')
		global DBFilename

		shutil.copy(DBFilename, Backfilename)

	@pyqtSlot()
	def on_actionRevert_to_previous_triggered(self):
		import shutil

		buttons = 'OKC'
		answer = informationMessage(self, 'Any work done since opening this instance of VGenes will be lost.', buttons)

		if answer == 'Cancel':
			return

		Backfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'BackUP.vdb')
		global DBFilename

		# DBFilename = filename
		self.GOOpen(False)
		shutil.move(Backfilename, DBFilename)

		self.GOOpen(False)

	@pyqtSlot()
	def on_btnExtractRecords_clicked(self):
		New = True
		self.MoveRecords(New)

	@pyqtSlot()
	def MoveRecords(self, New):

		IgBLASTAnalysis = []
		fields = ['*']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		if len(DataIs) > 0:
			for item in DataIs:
				IgBLASTAnalysis.append(list(item[:119]))

			if New == True:
				filename = saveFile(self, 'db')
			else:
				filename = openFile(self, 'db')

			answer3 = 'No'

			Startprocessed = len(IgBLASTAnalysis)
			if Startprocessed == 0: return
			if New == True:
				VGenesSQL.creatnewDB(filename)

			Processed, answer = VGenesSQL.enterData(self, filename, IgBLASTAnalysis, answer3)

			# todo need code to verify a database is open before you can import sequences.

			msg = 'Open ' + filename + '?'
			buttons = 'YN'
			answer = questionMessage(self, msg, buttons)
			if answer == 'Yes':
				global DBFilename
				DBFilename = filename
				self.GOOpen(False)

	@pyqtSlot()
	def on_chkSelectAllProt_clicked(self):
		if self.ui.chkSelectAllProt.isChecked():
			self.ui.chkFlexibility.setChecked(True)
			self.ui.chkHydrophilicity.setChecked(True)
			self.ui.chkHydrophobicity.setChecked(True)
			self.ui.chkInstability.setChecked(True)
			self.ui.chkpI.setChecked(True)
			self.ui.chkSurface.setChecked(True)
		else:
			self.ui.chkFlexibility.setChecked(False)
			self.ui.chkHydrophilicity.setChecked(False)
			self.ui.chkHydrophobicity.setChecked(False)
			self.ui.chkInstability.setChecked(False)
			self.ui.chkpI.setChecked(False)
			self.ui.chkSurface.setChecked(False)

	@pyqtSlot()
	def on_btnGenerateReport_clicked(self):
		# get info and seqs from checked
		# build text file and colormap to decorate or just CSV of colormap...
		# use 'repaired' aa sequence for color mapping...see decoratepeptide for example
		# as in:      1,0,.5,0.7,3.2,\n

		# first get list of seqs and info as tuple from DB:
		fields = ['SeqName', 'Sequence', 'GermlineSequence', 'CDR3Length', 'CDR1From', 'CDR1To', 'CDR2From', 'CDR2To',
		          'CDR3beg', 'CDR3end', 'Mutations', 'IDEvent', 'ID', 'Species', 'Jend']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		#  make array of fixed seqs to decorate
		# then run through DNASeq and make arrays for each decoration
		# build document with delineated regions and header line...need to run through
		# datais first to see where regions begin and end to get proper spacing...highest for each section (devided by 3)
		CDR1beg = 0
		CDR1end = 0
		CDR2beg = 0
		CDR2end = 0
		CDR3beg = 0
		CDR3end = 0
		CDR1len = 0
		CDR2len = 0
		CDR3len = 0
		FW1len = 0
		FW2len = 0
		FW3len = 0
		FW4len = 0

		NameLength = 0
		SeqLength = 0
		SeqArray = []
		AllSeqs = []

		# SeqArray has: SeqName, CDR1beg, CDR1end, CDR2beg, CDR2end, CDR3beg, CDR3end,
		for item in DataIs:
			SeqArray.clear()
			SeqName = item[0]
			SeqArray.append(SeqName)

			# make CDR1beg, CDR1end, just 3 Cs and NameLength

			DNASeq = item[1]
			GDNAseq = item[2]
			mutations = item[10]
			IDEvents = item[11]

			# unfixed version

			AASeq, ErMessage = VGenesSeq.Translator(DNASeq, 0)

			if IDEvents == 'Insertion' or IDEvents == 'Both':
				mutate = mutations
				mutations = mutate.split(',')
				for mut in mutations:
					if mut[:9] == 'Insertion':
						Ievent = mut
						Iparts = Ievent.split('-')
						AddAt = int(Iparts[1])
						SeqToAdd = Iparts[2]
						GDNAseq = GDNAseq[:AddAt] + SeqToAdd + GDNAseq[AddAt:]

			GAASeq, ErMessage = VGenesSeq.Translator(GDNAseq, 0)

			if int(item[4]) == 0 or int(item[5]) == 0 or int(item[6]) == 0 or int(item[7]) == 0 or int(item[8]) == 0:
				GCDRs = IgBLASTer.GetGLCDRs(GDNAseq, item[13])

			# GCDRs = [FW1beg, F1e,CDR1b, CDR1end, f2b,f2e, c2b.c2e,f3b,f3e]

			if int(item[4]) != 0:
				SeqArray.append((int(item[4]) - 1) / 3)  # 'c1b'
			else:
				SeqArray.append((int(GCDRs[2]) - 1) / 3)
			if int(item[5]) != 0:
				SeqArray.append((int(item[5])) / 3)  # c1e
			else:
				SeqArray.append(int(GCDRs[3]) / 3)

			if int(item[6]) != 0:
				SeqArray.append((int(item[6]) - 1) / 3)
			else:
				SeqArray.append((int(GCDRs[6]) - 1) / 3)

			if int(item[7]) != 0:
				SeqArray.append((int(item[7])) / 3)
			else:
				SeqArray.append(int(GCDRs[7]) / 3)

			if int(item[8]) != 0:
				SeqArray.append((int(item[8])) / 3)
			else:
				SeqArray.append(int(GCDRs[9]) / 3)

			if int(item[9]) != 0:
				SeqArray.append((int(item[9])) / 3)
			else:
				SeqArray.append(len(GAASeq))

			if int(item[9]) != 0:
				Jend = int(item[14]) / 3
			else:
				SeqArray.append(len(GAASeq))
			# SeqArray has: SeqName, CDR1beg, CDR1end, CDR2beg, CDR2end, CDR3beg, CDR3end,

			CDR1beg = int(SeqArray[1])
			CDR1end = int(SeqArray[2])
			CDR2beg = int(SeqArray[3])
			CDR2end = int(SeqArray[4])
			CDR3beg = int(SeqArray[5])
			CDR3end = int(SeqArray[6])

			if len(SeqName) > NameLength: NameLength = len(SeqName)

			if len(AASeq) > len(GAASeq):
				LenTo = len(GAASeq)
				AASeq = AASeq[:LenTo]
			else:
				LenTo = len(AASeq)

			SeqArray.append(AASeq)  # place original sequence without bad germ and seq regions for alignment

			for i in range(0, LenTo - 1):  # first replace bad codons with germline codons
				if AASeq[i] == GAASeq[i]:
					if AASeq[i] == '.' or AASeq[i] == '~':
						AASeq = AASeq[:i] + AASeq[i + 1:] + '.'
						GAASeq = GAASeq[:i] + GAASeq[i + 1:] + '.'

			for i in range(0, LenTo - 1):
				if AASeq[i] != GAASeq[i]:
					if AASeq[i] == '.' or AASeq[i] == '~' or AASeq[i] == '*':
						AASeq = AASeq[:i] + GAASeq[i] + AASeq[i + 1:]

			AASeq = AASeq.replace('~', '').replace('.', '')

			if len(AASeq) > SeqLength: SeqLength = len(AASeq)

			if CDR1beg > FW1len: FW1len = CDR1beg
			if (CDR1end - CDR1beg) > CDR1len: CDR1len = (CDR1end - CDR1beg)
			if (CDR2beg - CDR1end) > FW2len: FW2len = (CDR2beg - CDR1end)
			if (CDR2end - CDR2beg) > CDR2len: CDR2len = (CDR2end - CDR2beg)
			if (CDR3beg - CDR2end) > FW3len: FW3len = (CDR3beg - CDR2end)
			if (CDR3end - CDR3beg) > CDR3len: CDR3len = (CDR3end - CDR3beg)
			if (Jend - CDR3end) > FW4len: FW4len = (Jend - CDR3end)

			if self.ui.chkHydrophobicity.isChecked() == True:

				WindowSize = self.ui.spnHydrophobicity.value()
				if WindowSize < 2:
					WindowSize = 2
					self.ui.spnHydrophobicity.setValue(2)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnHydrophobicity.setValue(len(AASeq) - 1)
				PhobCurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'Hydrophobicity', WindowSize, True)

				PhobScale = (-4.5, 4.5)  # based on tests paramators

				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			if self.ui.chkHydrophilicity.isChecked() == True:
				WindowSize = self.ui.spnHydrophilicity.value()
				if WindowSize < 2:
					WindowSize = 2
					self.ui.spnHydrophilicity.setValue(2)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnHydrophilicity.setValue(len(AASeq) - 1)

				PhilCurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'Hydrophilicity', WindowSize, True)

				PhilScale = (-3.4, 3.0)  # based on tests paramators
				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			if self.ui.chkFlexibility.isChecked() == True:
				WindowSize = self.ui.spnFlexibility.value()
				if WindowSize < 9:
					WindowSize = 9
					self.ui.spnFlexibility.setValue(9)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnFlexibility.setValue(len(AASeq) - 1)
				FlexCurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'Flexibility', WindowSize, True)

				FlexScale = (0.904, 1.102)  # based on tests paramators
				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			if self.ui.chkSurface.isChecked() == True:
				WindowSize = self.ui.spnSurface.value()
				if WindowSize < 2:
					WindowSize = 2
					self.ui.spnSurface.setValue(2)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnSurface.setValue(len(AASeq) - 1)
				SurfCurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'Surface', WindowSize, True)

				SurfScale = (0.394, 1.545)  # based on tests paramators
				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			if self.ui.chkpI.isChecked() == True:
				WindowSize = self.ui.spnpI.value()
				if WindowSize < 2:
					WindowSize = 2
					self.ui.spnpI.setValue(2)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnpI.setValue(len(AASeq) - 1)
				pICurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'MapAApI', WindowSize, True)

				pIScale = (0, 14)  # based on tests paramators
				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			if self.ui.chkInstability.isChecked() == True:
				WindowSize = self.ui.spnInstability.value()
				if WindowSize < 8:
					WindowSize = 8
					self.ui.spnInstability.setValue(8)
				elif WindowSize > len(AASeq) - 1:
					WindowSize = len(AASeq) - 1
					self.ui.spnInstability.setValue(len(AASeq) - 1)
				InsCurPos = (WindowSize // 2)
				ColorMap = VGenesSeq.OtherParam(AASeq, 'MapInstability', WindowSize, True)

				# for this need to scale relatively but so that anything>40 is in the red  as 40+ = unstable
				if ColorMap != 0:
					Highest = max(ColorMap)
					Lowest = min(ColorMap)
					maxi = ((40 - Lowest) / 8) * 11
					InsScale = (Lowest, maxi)  # based on tests paramators
				SeqArray.append(ColorMap)
			else:
				SeqArray.append('None')

			AllSeqs.append(tuple(SeqArray))

		nameIS = '{message: <{width}}'.format(message='', width=NameLength)
		FW1 = '{message: <{width}}'.format(message='FWR1', width=FW1len)
		CW1 = '{message: <{width}}'.format(message='CDR1', width=CDR1len)
		FW2 = '{message: <{width}}'.format(message='FWR2', width=FW2len)
		CW2 = '{message: <{width}}'.format(message='CDR2', width=CDR2len)
		FW3 = '{message: <{width}}'.format(message='FWR3', width=FW3len)
		CW3 = '{message: <{width}}'.format(message='CDR3', width=CDR3len)
		FW4 = '{message: <{width}}'.format(message='FWR4', width=FW4len)
		HeaderLine = nameIS + ' | ' + FW1 + ' | ' + CW1 + ' | ' + FW2 + ' | ' + CW2 + ' | ' + FW3 + ' | ' + CW3 + ' | ' + FW4 + ' |\n'

		HLen = len(HeaderLine)

		# Colorhead = '{message:0<{width}}'.format(message='', width=HLen)

		NameLength = len(nameIS)

		HydroPhob = []
		HydroPhil = []
		flexi = []
		surface = []
		pI = []
		Instab = []
		CurrentSeq = ''
		if self.ui.chkHydrophobicity.isChecked() == True:
			for i in range(0, HLen - 1):
				HydroPhob.append(0)

		if self.ui.chkHydrophilicity.isChecked() == True:
			for i in range(0, HLen - 1):
				HydroPhil.append(0)

		if self.ui.chkFlexibility.isChecked() == True:
			for i in range(0, HLen - 1):
				flexi.append(0)

		if self.ui.chkSurface.isChecked() == True:
			for i in range(0, HLen - 1):
				surface.append(0)

		if self.ui.chkpI.isChecked() == True:
			for i in range(0, HLen - 1):
				pI.append(0)

		if self.ui.chkInstability.isChecked() == True:
			for i in range(0, HLen - 1):
				Instab.append(0)

		SeqSet = HeaderLine

		for Sequence in AllSeqs:
			PC1b = int(Sequence[1])
			PC1e = int(Sequence[2])
			PC2b = int(Sequence[3])
			PC2e = int(Sequence[4])
			PC3b = int(Sequence[5])
			PC3e = int(Sequence[6])

			AASeq = Sequence[7]

			SeqName = Sequence[0]
			SeqName.replace('|', '-')
			nameIS = '{message: <{width}}'.format(message=SeqName, width=NameLength)

			Segment = AASeq[0:PC1b]

			FW1 = '{message: <{width}}'.format(message=Segment, width=FW1len)
			FW1len2 = len(Segment)
			Segment = AASeq[PC1b:PC1e]
			CW1 = '{message: <{width}}'.format(message=Segment, width=CDR1len)
			CDR1len2 = len(Segment)
			Segment = AASeq[PC1e:PC2b]
			FW2 = '{message: <{width}}'.format(message=Segment, width=FW2len)
			FW2len2 = len(Segment)
			Segment = AASeq[PC2b:PC2e]
			CW2 = '{message: <{width}}'.format(message=Segment, width=CDR2len)
			CDR2len2 = len(Segment)
			Segment = AASeq[PC2e:PC3b]
			FW3 = '{message: <{width}}'.format(message=Segment, width=FW3len)
			FW3len2 = len(Segment)
			Segment = AASeq[PC3b:PC3e]
			CW3 = '{message: <{width}}'.format(message=Segment, width=CDR3len)
			CDR3len2 = len(Segment)
			Segment = AASeq[PC3e:]
			FW4 = '{message: <{width}}'.format(message=Segment, width=FW4len)
			FW4len2 = len(Segment)

			StartNew = len(SeqSet)
			CurrentSeq = nameIS + ' | ' + FW1 + ' | ' + CW1 + ' | ' + FW2 + ' | ' + CW2 + ' | ' + FW3 + ' | ' + CW3 + ' | ' + FW4 + ' |\n'
			SeqSet += CurrentSeq
			LenSeq = len(CurrentSeq)
			CurPos = ''

			for i in range(8, 14):

				if Sequence[i] != 'None':
					ColorPart = Sequence[i]

					# make list
					n = 0
					j = 0
					if i == 8:

						FW1len2 -= PhobCurPos

						while len(HydroPhob) != StartNew - 1:
							HydroPhob.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							HydroPhob.append(0)
							StartNew += 1

						for j in range(0, 2):  # for window size
							HydroPhob.append(0)
							StartNew += 1

						for j in range(0, FW1len2):
							if ColorPart != 0:
								HydroPhob.append(ColorPart[j])
							else:
								HydroPhob.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - PhobCurPos) + 3):  # for end of match and spacer
							HydroPhob.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								HydroPhob.append(ColorPart[j])
							else:
								HydroPhob.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							HydroPhob.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								HydroPhob.append(ColorPart[j])
							else:
								HydroPhob.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							HydroPhob.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								HydroPhob.append(ColorPart[j])
							else:
								HydroPhob.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							HydroPhob.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									HydroPhob.append(ColorPart[j])
								else:
									HydroPhob.append(0)
							except:
								print('tried2973')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							HydroPhob.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									HydroPhob.append(ColorPart[j])
								else:
									HydroPhob.append(0)
								n += 1
							except:
								print('tried2987')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							HydroPhob.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									HydroPhob.append(ColorPart[j])
								else:
									HydroPhob.append(0)
								n += 1
							except:
								print('done')


								# while len(HydroPhob) < LenSeq:
								#     HydroPhob.append(0)


					elif i == 9:

						FW1len2 -= 2

						while len(HydroPhil) != StartNew - 1:
							HydroPhil.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							HydroPhil.append(0)
							StartNew += 1

						for j in range(0, 2):  # for window size
							HydroPhil.append(0)
							StartNew += 1

						for j in range(0, FW1len2):
							if ColorPart != 0:
								HydroPhil.append(ColorPart[j])
							else:
								HydroPhil.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - PhilCurPos) + 3):  # for end of match and spacer
							HydroPhil.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								HydroPhil.append(ColorPart[j])
							else:
								HydroPhil.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							HydroPhil.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								HydroPhil.append(ColorPart[j])
							else:
								HydroPhil.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							HydroPhil.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								HydroPhil.append(ColorPart[j])
							else:
								HydroPhil.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							HydroPhil.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									HydroPhil.append(ColorPart[j])
								else:
									HydroPhil.append(0)
							except:
								print('tried3075')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							HydroPhil.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									HydroPhil.append(ColorPart[j])
								else:
									HydroPhil.append(0)
								n += 1
							except:
								print('tried3089')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							HydroPhil.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									HydroPhil.append(ColorPart[j])
								else:
									HydroPhil.append(0)
								n += 1
							except:
								print('done')

								# while len(HydroPhil) < len(SeqSet):
								#     HydroPhil.append(0)


					elif i == 10:

						FW1len2 -= 2

						while len(flexi) != StartNew - 1:
							flexi.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							flexi.append(0)
							StartNew += 1

						for j in range(0, 2):  # for window size
							flexi.append(0)
							StartNew += 1

						for j in range(0, FW1len2):
							if ColorPart != 0:
								flexi.append(ColorPart[j])
							else:
								flexi.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - FlexCurPos) + 3):  # for end of match and spacer
							flexi.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								flexi.append(ColorPart[j])
							else:
								flexi.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							flexi.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								flexi.append(ColorPart[j])
							else:
								flexi.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							flexi.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								flexi.append(ColorPart[j])
							else:
								flexi.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							flexi.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									flexi.append(ColorPart[j])
								else:
									flexi.append(0)
							except:
								print('tried3176')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							flexi.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									flexi.append(ColorPart[j])
								else:
									flexi.append(0)
								n += 1
							except:
								print('tried3190')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							flexi.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									flexi.append(ColorPart[j])
								else:
									flexi.append(0)
								n += 1
							except:
								print('done')
								# while len(flexi) < len(SeqSet):
								#     flexi.append(0)


					elif i == 11:

						FW1len2 -= 2

						while len(surface) != StartNew - 1:
							surface.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							surface.append(0)
							StartNew += 1

						for j in range(0, 2):  # for window size
							surface.append(0)
							StartNew += 1

						for j in range(0, FW1len2):
							if ColorPart != 0:
								surface.append(ColorPart[j])
							else:
								surface.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - SurfCurPos) + 3):  # for end of match and spacer
							surface.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								surface.append(ColorPart[j])
							else:
								surface.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							surface.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								surface.append(ColorPart[j])
							else:
								surface.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							surface.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								surface.append(ColorPart[j])
							else:
								surface.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							surface.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									surface.append(ColorPart[j])
								else:
									surface.append(0)
							except:
								print('tried3276')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							surface.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									surface.append(ColorPart[j])
								else:
									surface.append(0)
								n += 1
							except:
								print('tried3290')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							surface.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									surface.append(ColorPart[j])
								else:
									surface.append(0)
								n += 1
							except:
								print('done')
								# while len(surface) < len(SeqSet):
								#     surface.append(0)

					elif i == 12:

						FW1len2 -= 2

						while len(pI) != StartNew - 1:
							pI.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							pI.append(0)
							StartNew += 1

						for j in range(0, 2):  # for window size
							pI.append(0)
							StartNew += 1

						for j in range(0, FW1len2):
							if ColorPart != 0:
								pI.append(ColorPart[j])
							else:
								pI.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - pICurPos) + 3):  # for end of match and spacer
							pI.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								pI.append(ColorPart[j])
							else:
								pI.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							pI.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								pI.append(ColorPart[j])
							else:
								pI.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							pI.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								pI.append(ColorPart[j])
							else:
								pI.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							pI.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									pI.append(ColorPart[j])
								else:
									pI.append(0)
							except:
								print('tried3374')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							pI.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									pI.append(ColorPart[j])
								else:
									pI.append(0)
								n += 1
							except:
								print('tried3388')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							pI.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									pI.append(ColorPart[j])
								else:
									pI.append(0)
								n += 1
							except:
								print('done')
								# while len(pI) < len(SeqSet):
								#     pI.append(0)

					elif i == 13:

						FW1len2 -= 2

						while len(Instab) != StartNew - 1:
							Instab.append(0)

						while CurPos != '|':
							CurPos = SeqSet[StartNew]
							Instab.append(0)
							StartNew += 1

						for j in range(0, InsCurPos):  # for window size
							Instab.append(0)
							StartNew += 1

						for j in range(0, 2):
							if ColorPart != 0:
								Instab.append(ColorPart[j])
							else:
								Instab.append(0)
							n += 1

						for j in range(FW1len2, (len(FW1) - InsCurPos) + 3):  # for end of match and spacer
							Instab.append(0)

						for j in range(n, n + CDR1len2):
							if ColorPart != 0:
								Instab.append(ColorPart[j])
							else:
								Instab.append(0)
							n += 1
						for j in range(CDR1len2, len(CW1) + 3):  # for spacer
							Instab.append(0)

						for j in range(n, n + FW2len2):
							if ColorPart != 0:
								Instab.append(ColorPart[j])
							else:
								Instab.append(0)
							n += 1
						for j in range(FW2len2, len(FW2) + 3):  # for spacer
							Instab.append(0)

						for j in range(n, n + CDR2len2):
							if ColorPart != 0:
								Instab.append(ColorPart[j])
							else:
								Instab.append(0)
							n += 1
						for j in range(CDR2len2, len(CW2) + 3):  # for spacer
							Instab.append(0)

						for j in range(n, n + FW3len2):
							try:
								if ColorPart != 0:
									Instab.append(ColorPart[j])
								else:
									Instab.append(0)
							except:
								print('tried3472')
							n += 1
						for j in range(FW3len2, len(FW3) + 3):  # for spacer
							Instab.append(0)

						for j in range(n, n + CDR3len2):
							try:
								if ColorPart != 0:
									Instab.append(ColorPart[j])
								else:
									Instab.append(0)
								n += 1
							except:
								print('tried3486')
						for j in range(CDR3len2, len(CW3) + 3):  # for spacer
							Instab.append(0)

						for j in range(n, n + FW4len2 - 1):
							try:
								if ColorPart != 0:
									Instab.append(ColorPart[j])
								else:
									Instab.append(0)
								n += 1
							except:
								print('done')

								# while len(Instab) < len(SeqSet):
								#     Instab.append(0)






								# self.DecorateText(FinalMap, Scale, CurPos, cursor)

								# FW1 = '{message:0<{width}}'.format(message=Segment, width=FW1len) #the 0 is fill character
		FinalDoc = 'Protein properties report\n\n'

		FinalMap = []  # '000000000000000000000000000'
		FinalMap.clear
		for i in range(0, len(FinalDoc)):
			FinalMap.append(0)

		CurPos = 0

		if self.ui.chkHydrophobicity.isChecked() == True:
			FinalDoc += 'Hydrophobicity:\n' + SeqSet + '\n'
			# testLen = len(FinalDoc)
		if self.ui.chkHydrophilicity.isChecked() == True:
			FinalDoc += 'Hydrophilicity:\n' + SeqSet + '\n'
			# testLen = len(FinalDoc)
		if self.ui.chkFlexibility.isChecked() == True:
			FinalDoc += 'Flexibility:\n' + SeqSet + '\n'

		if self.ui.chkSurface.isChecked() == True:
			FinalDoc += 'Surface liklihood:\n' + SeqSet + '\n'

		if self.ui.chkpI.isChecked() == True:
			FinalDoc += 'Isoelectric point (pI):\n' + SeqSet + '\n'

		if self.ui.chkInstability.isChecked() == True:
			FinalDoc += 'Instability:\n' + SeqSet + '\n'

		FinalDoc += 'Scale: Low-> -5|-4|-3|-2|-1| 0 |+1|+2|+3|+4|+5  ->high'

		if self.ui.chkShowInEditor.isChecked() == True:
			Style = 'ProteinReport'
			self.ShowVGenesTextEdit(FinalDoc, Style)
			cursor = self.TextEdit.textEdit.textCursor()
		else:
			self.ui.txtProtein.setText(FinalDoc)
			cursor = self.ui.txtProtein.textCursor()

		if self.ui.chkHydrophobicity.isChecked() == True:

			for i in range(0, 16):
				FinalMap.append(0)

			for item in HydroPhob:
				FinalMap.append(item)

			CurPos += PhobCurPos
			Scale = PhobScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			CurPos += len(FinalMap)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1

			FinalMap.clear()

		if self.ui.chkHydrophilicity.isChecked() == True:

			for i in range(0, 16):
				FinalMap.append(0)
				# CurPos += 1

			for item in HydroPhil:
				FinalMap.append(item)

			testLen = len(FinalMap)
			CurPos += PhilCurPos
			Scale = PhilScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1
			FinalMap.clear()

		if self.ui.chkFlexibility.isChecked() == True:

			for i in range(0, 13):
				FinalMap.append(0)

			for item in flexi:
				FinalMap.append(item)

			CurPos += FlexCurPos
			Scale = FlexScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1
			FinalMap.clear()

		if self.ui.chkSurface.isChecked() == True:

			for i in range(0, 19):
				FinalMap.append(0)

			for item in surface:
				FinalMap.append(item)

			CurPos += SurfCurPos
			Scale = SurfScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1
			FinalMap.clear()

		if self.ui.chkpI.isChecked() == True:

			for i in range(0, 24):
				FinalMap.append(0)

			for item in pI:
				FinalMap.append(item)

			CurPos += pICurPos
			Scale = pIScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1
			FinalMap.clear()

		if self.ui.chkInstability.isChecked() == True:

			for i in range(0, 13):
				FinalMap.append(0)

			for item in Instab:
				FinalMap.append(item)

			CurPos += InsCurPos
			Scale = InsScale

			self.DecorateText(FinalMap, Scale, CurPos, cursor)
			testS = FinalDoc[CurPos]
			while testS != '|':
				testS = FinalDoc[CurPos]
				CurPos += 1
			CurPos += 1
			FinalMap.clear()

		Scale = (-5, 5)
		CurPos = len(FinalDoc) - 55
		# 'Scale: Low-> -5|-4|-3|-2|-1| 0 |+1|+2|+3|+4|+5  ->high'
		for i in range(0, 14):
			FinalMap.append(0)  # 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,10,0,0,0,0,0,0,0)
		for i in range(-5, 6):
			for j in range(0, 3):
				FinalMap.append(i)
		FinalMap.append(10)
		for i in range(0, 7):
			FinalMap.append(0)

		self.DecorateText(FinalMap, Scale, CurPos, cursor)



		# SeqArray has: SeqName, CDR1beg, CDR1end, CDR2beg, CDR2end, CDR3beg, CDR3end,
		# AAseq, Hydrophobicity, Hyrdorphilicity, flexibility, surface, pI, Instability,

	# now for decorated seqs build document in VGenetext editor (no wrap)
	# for fun add scale line made of spaces that are colormapped, can even give numbers instead spaces
	# and numerical color maps using VGenesSeq.OtherParam and then
	# color using VGenes text edit cursor and    self.DecorateText(ColorMap, Scale, CurPos, cursor)
	# where CurPos will allow iteration through seq with different map conditions





	@pyqtSlot()
	def on_btnUpdateTree_clicked(self):
		fields = self.ui.cboTreeOp1.currentText()
		field1Index = self.ui.cboTreeOp1.currentIndex()

		value = self.ui.treeWidget.selectedItems()
		currentitemIs = ''

		for item in value:
			currentitemIs = item.text(0)

		field1 = self.TransLateFieldtoReal(fields, True)

		fields = self.ui.cboTreeOp2.currentText()

		field2 = self.TransLateFieldtoReal(fields, True)
		field2Index = self.ui.cboTreeOp2.currentIndex()

		fields = self.ui.cboTreeOp3.currentText()

		# self.ui.lblTreeOtions3.setText(data[self.ui.cboTreeOp3.currentIndex()])
		field3 = self.TransLateFieldtoReal(fields, True)
		field3Index = self.ui.cboTreeOp3.currentIndex()
		if field1 == '': field1 = 'None'
		if field2 == '': field1 = 'None'
		if field3 == '': field1 = 'None'
		SQLFields = []
		SQLFields.append(field1)
		SQLFields.append(field2)
		SQLFields.append(field3)

		value = self.ui.treeWidget.selectedItems()
		currentitemIs = ''

		for item in value:
			currentitemIs = item.text(0)

		model = self.ui.tableView.model()

		# global RefreshSQL
		if field1 == 'None' or field1 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY SeqName'
		elif field2 == 'None' or field2 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', SeqName'
		elif field3 == 'None' or field3 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', SeqName'
		else:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', ' + field3 + ', SeqName'

		model.refresh()

		self.ui.tableView.sortByColumn(field1Index, Qt.AscendingOrder)
		self.ui.tableView.sortByColumn(0, Qt.AscendingOrder)
		self.ui.tableView.sortByColumn(field3Index, Qt.AscendingOrder)
		self.ui.tableView.sortByColumn(field2Index, Qt.AscendingOrder)

		self.initializeTreeView(SQLFields)
		global DontFindTwice
		DontFindTwice = True
		self.findTreeItem(currentitemIs)
		# self.ui.treeWidget.
		DontFindTwice = False

	def findTreeItem(self, ChildName):

		found = self.ui.treeWidget.findItems(ChildName, Qt.MatchRecursive, 0)

		for item in found:
			currentRecord = item.text(0)
			if len(found) > 1: break

			# print(currentRecord)
		if len(found) > 0: self.ui.treeWidget.setCurrentItem(item)
		try:
			return item
		except:
			print("line 7598 exception")

	def on_radioButtonSeqView_clicked(self):
		global JustMoved
		JustMoved = True
		self.ui.txtDNASeq.setText(data[79])

		try:
			self.SeqButton(LastPushed)
		except:
			self.SeqButton('v')

	def updateF(self, nID):
		ID =int(nID)
		if ID == -1 or ID == -2:
			global PreVID
			PreVID = 0

		global JustMovedIt
		global FirstupdateF
		JustMovedIt = True
		# JustMoved = True
		if PreVID != ID:
			# todo need to refill data whenever a record changes

			if ID != -1:
				data.clear()

				if FirstupdateF == False and ID > -1:
					# MatchingIndex = NameIndex[name]
					newID = list(NameIndex.keys())[list(NameIndex.values()).index(ID)]
					SQLStatement = 'SELECT * FROM vgenesDB WHERE SeqName = "' + str(
						newID) + '"'  # VGenesSQL.MakeSQLStatement(self, fields, data[0])
					DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
					for record in DataIs:
						for item in record:
							data.append(str(item))

				else:
					model = self.ui.tableView.model()
					if ID == -2: ID = 0

					for i in range(0, 120):
						index = model.index(ID, i)
						data.append(str(model.data(index)))

					# global FirstupdateF
					FirstupdateF = False

				PreVID = ID

				self.ui.txtName.setText(data[0])
				self.ui.txtName_2.setText(data[0])
				self.ui.label_Name.setText(data[0])

				self.ui.txtVgene.setText(data[3])
				self.ui.txtDgene.setText(data[6])
				self.ui.txtJgene.setText(data[9])
				self.ui.txtVGeneSeq.setText(data[3])
				self.ui.txtDGeneSeq.setText(data[6])
				self.ui.txtJGeneSeq.setText(data[9])
				self.ui.txtStop.setText(data[12])

				self.ui.txtVLocus.setText(data[90])
				self.ui.txtDLocus.setText(data[92])
				self.ui.txtJLocus.setText(data[91])

				self.ui.txtProject.setText(data[75])

				self.ui.txtProject_2.setText(data[75])
				self.ui.txtGroup.setText(data[76])
				self.ui.txtGroup_2.setText(data[76])

				self.ui.cboTreeOp1.setCurrentText(data[75])
				self.ui.cboTreeOp2.setCurrentText(data[76])
				self.ui.cboTreeOp3.setCurrentText(data[77])

				self.ui.txtSubGroup.setText(data[77])
				self.ui.txtSubGroup_2.setText(data[77])

				# self.ui.txtStop.setText(data[12])

				self.ui.txtReadingFrame.setText(data[13])
				self.ui.txtProductive.setText(data[14])
				self.ui.txtVend.setText(data[16])
				self.ui.txtVD.setText(data[17])
				self.ui.txtD.setText(data[18])

				self.ui.txtDJ.setText(data[19])
				self.ui.txtJend.setText(data[20])
				self.ui.txtSeqAlignment.setText(data[58])
				global JustMoved
				JustMoved = True

				if self.ui.btnEditSeq.isChecked():
					msg = 'Sequence edit mode was activated, do you want to save changes and re-analyze this sequence before proceeding?'
					buttons = 'YN'
					answer = questionMessage(self, msg, buttons)
					if answer == 'Yes':
						self.ui.btnEditSeq.setChecked(False)
						self.UpdateSeqAnalysis()

						self.ui.btnEditSeq.setText("Edit Mode")
						msg = 'Prese "Edit Mode" to edit sequence.'

						self.ui.lblSeq2.setText(msg)

					elif answer == 'No':
						self.ui.btnEditSeq.setChecked(False)
						self.ui.btnEditSeq.setText("Edit Mode")
						msg = 'Prese "Edit Mode" to edit sequence.'




				self.ui.txtDNASeq.setText(data[79])




				self.SeqButton(LastPushed)
				# self.ui.txtDNASeq.setReadOnly(True)




				JustMoved = False
				self.ui.txtCDR3DNA.setText(data[81])
				self.ui.txtCDR3AA.setText(data[82])
				self.ui.txtCDR3Length.setText(data[83])

				# if data[88] is int:
				self.ui.txtClonalPool.setPlainText(data[88])
				# if data[89] is int:
				self.ui.txtClonalRank.setPlainText(data[89])
				self.ui.txtDateTime.setText(data[93])
				self.ui.txtQuality.setText(data[95])
				self.ui.txtQuality_2.setText(data[95])
				self.ui.txtComments.setText(data[94])
				self.ui.txtID.setText(data[98])
				self.ui.txtCDR3MW.setText(data[99])
				self.ui.txtCDR3pI.setText(data[100])
				self.ui.txtIsotype.setText(data[101])
				self.ui.txtIsotypeSeq.setText(data[101])
				# self.IncrementDials()


				# if DontFindTwice == False:
				#     self.findTreeItem(data[0])
				currentRecord = self.ui.tableView.currentIndex().row()
				maxRecords = self.ui.tableView.model().rowCount()
				self.ui.horizontalScrollBar.setMaximum(maxRecords)
				self.ui.dial.setMaximum(maxRecords)

				self.ui.lcdNumber_max.display(maxRecords)
				self.ui.horizontalScrollBar.setValue(currentRecord)
				self.ui.dial.setValue(currentRecord)
				currentRecord += 1
				self.ui.lcdNumber_current.display(currentRecord)



				# self.ui.tableViewFeatures.setModel(model)

				# self.ui.tableViewFeatures.setModel(tabledata, header, self)
				# self.ui.tableView.setCurrentIndex(ID)

				JustMovedIt = False




			#  0 SeqName, 1 SeqLen, 2 GeneType, 3 V1, 4 V2, 5 V3, 6 D1, 7 D2, 8 D3, 9 J1, 10 J2, 11 J3,
			# 12 StopCodon, 13 ReadingFrame, 14 productive, 15 Strand, 16 VSeqend, 17 VDJunction,
			# 18 Dregion, 19 DJJunction, 20 begJ, 21 VJunction, 22 FR1From, 23 FR1To, 24 FR1length,
			# 25 FR1matches, 26 FR1mis, 27 FR1gaps, 28 FR1PercentIdentity, 29 CDR1From, 30 CDR1To, 
			# 31 CDR1length, 32 CDR1matches, 33 CDR1mis, 34 CDR1gaps, 35 CDR1PercentIdentity,
			# 36 FR2From, 37 FR2To, 38 FR2length, 39 FR2matches, 40 FR2mis, 41 FR2gaps, 42 FR2PercentIdentity,
			# 43 CDR2From, 44 CDR2To, 45 CDR2length, 46 CDR2matches, 47 CDR2mis, 48 CDR2gaps, 49 CDR2PercentIdentity,
			# 50 FR3From, 51 FR3To, 52 FR3length, 53 FR3matches, 54 FR3mis, 55 FR3gaps, 56 FR3PercentIdentity, 57 TotMut, 
			# 58 SeqAlignment, 59 GVbeg, 60 GVend, 61 GD1beg, 62 GD1end, 63 GD2beg, 64 GD2end, 65 GJbeg,
			#  66 GJend, 67 Vbeg, 68 Vend, 69 D1beg, 70 D1end, 71 D2beg, 72  D2end, 73 Jbeg, 74 Jend, 75 Project,
			#  76 Grouping, 77 SubGroup, 78 Species, 79 Sequence, 80 GermlineSequence, 81 CDR3DNA, 82 CDR3AA,
			#  83 CDR3Length, 84 CDR3beg, 85 CDRend, 86 Specificity, 87 Subspecificity, 88 ClonalPool, 89 ClonalRank,
			#  90 VLocus, 91 JLocus, 92 DLocus, 93 DateEntered, 94 Comments, 95 Quality, 96 TotalMuts, 97 Mutations, 98 IDEvent, CDR3MW, CDR3pI, Isotype, Blank4, Blank5, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, 99 ID

	# @pyqtSlot()
	# def IncrementDials(self):
	#
	#     if DontFindTwice == False:
	#         self.findTreeItem(data[0])
	#     currentRecord = self.ui.tableView.currentIndex().row()
	#     maxRecords = self.ui.tableView.model().rowCount()
	#     self.ui.horizontalScrollBar.setMaximum(maxRecords)
	#     self.ui.dial.setMaximum(maxRecords)
	#
	#     self.ui.lcdNumber_max.display(maxRecords)
	#     self.ui.horizontalScrollBar.setValue(currentRecord)
	#     self.ui.dial.setValue(currentRecord)
	#     currentRecord += 1
	#     self.ui.lcdNumber_current.display(currentRecord)



	@pyqtSlot()
	def on_actionCreateVDJdb_triggered(self):

		Pathname = openFile(self, 'FASTA')
		VGenesSQL.CreateVDJDB(Pathname)

	@pyqtSlot()
	def on_btnV_clicked(self):
		global LastPushed
		self.SeqButton('v')
		LastPushed = 'v'

	@pyqtSlot()
	def on_btnD_clicked(self):
		global LastPushed
		LastPushed = 'd'
		self.SeqButton('d')

	@pyqtSlot()
	def on_btnJ_clicked(self):
		global LastPushed
		LastPushed = 'j'
		self.SeqButton('j')

	@pyqtSlot()
	def on_btnFW_1_clicked(self):
		global LastPushed
		LastPushed = 'f1'
		self.SeqButton('f1')

	@pyqtSlot()
	def on_btnCW_1_clicked(self):
		global LastPushed
		LastPushed = 'c1'
		self.SeqButton('c1')

	@pyqtSlot()
	def on_btnFW_2_clicked(self):
		global LastPushed
		LastPushed = 'f2'
		self.SeqButton('f2')

	@pyqtSlot()
	def on_btnCW_2_clicked(self):
		global LastPushed
		LastPushed = 'c2'
		self.SeqButton('c2')

	@pyqtSlot()
	def on_btnFW_3_clicked(self):
		global LastPushed
		LastPushed = 'f3'
		self.SeqButton('f3')

	@pyqtSlot()
	def on_btnCW_3_clicked(self):
		global LastPushed
		LastPushed = 'c3'
		self.SeqButton('c3')

	@pyqtSlot()
	def on_btnFW_4_clicked(self):
		global LastPushed
		LastPushed = 'f4'
		self.SeqButton('f4')

	@pyqtSlot()
	def on_btnVDJ_clicked(self):
		global LastPushed
		LastPushed = 'vdj'
		self.SeqButton('vdj')

	@pyqtSlot()
	def on_btnC_clicked(self):
		global LastPushed
		LastPushed = 'c'
		self.SeqButton('c')

	@pyqtSlot()
	def SeqButton(self, button):
		global JustMoved
		cursor = self.ui.txtDNASeq.textCursor()
		AAcursor = self.ui.txtAASeq.textCursor()
		StartSel = 0
		EndSel = 0
		if button == 'v' or button == 'd' or button == 'j':
			self.ui.btnFW_1.setChecked(False)
			self.ui.btnCW_1.setChecked(False)
			self.ui.btnFW_2.setChecked(False)
			self.ui.btnCW_2.setChecked(False)
			self.ui.btnFW_3.setChecked(False)
			self.ui.btnCW_3.setChecked(False)
			self.ui.btnFW_4.setChecked(False)
			self.ui.btnVDJ.setChecked(False)
			self.ui.btnC.setChecked(False)

			if self.ui.btnV.isChecked() and self.ui.btnD.isChecked() and self.ui.btnJ.isChecked():
				self.ui.btnVDJ.setChecked(True)
				self.ui.btnFW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnCW_3.setChecked(True)
				self.ui.btnFW_4.setChecked(True)

				StartSel = int(data[67]) - 1
				EndSel = int(data[74])
			elif self.ui.btnV.isChecked() and self.ui.btnJ.isChecked():
				self.ui.btnD.setChecked(True)
				self.ui.btnVDJ.setChecked(True)
				self.ui.btnFW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnCW_3.setChecked(True)
				self.ui.btnFW_4.setChecked(True)

				StartSel = int(data[67]) - 1
				EndSel = int(data[74])



			elif self.ui.btnV.isChecked() and self.ui.btnD.isChecked():
				StartSel = int(data[67]) - 1
				EndSel = int(data[70])

			elif self.ui.btnD.isChecked() and self.ui.btnJ.isChecked():
				StartSel = int(data[69]) - 1
				EndSel = int(data[74])


			else:  # then not a combination
				if self.ui.btnV.isChecked():
					StartSel = int(data[67]) - 1
					EndSel = int(data[68])
				if self.ui.btnD.isChecked():
					StartSel = int(data[69]) - 1
					EndSel = int(data[70])
				if self.ui.btnJ.isChecked():
					StartSel = int(data[73]) - 1
					EndSel = int(data[74])
				else:
					cursor.setPosition(0)

					JustMoved = True
					self.ui.txtDNASeq.setTextCursor(cursor)
					self.ui.txtAASeq.setTextCursor(AAcursor)
					AAcursor.setPosition(0)
					AAcursor.setPosition(len(self.ui.txtAASeq.toPlainText()), QTextCursor.KeepAnchor)
					format = QTextCharFormat()
					format.setFontUnderline(False)
					AAcursor.mergeCharFormat(format)
					AAcursor.setPosition(0)
					JustMoved = False


		elif button == 'f1' or button == 'c1' or button == 'f2' or button == 'c2' or button == 'f3' or button == 'c3' or button == 'f4':
			self.ui.btnV.setChecked(False)
			self.ui.btnD.setChecked(False)
			self.ui.btnJ.setChecked(False)
			self.ui.btnC.setChecked(False)
			self.ui.btnVDJ.setChecked(False)

			if self.ui.btnFW_1.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[23]) - 1
				EndSel = int(data[74])
				self.ui.btnV.setChecked(True)
				self.ui.btnD.setChecked(True)
				self.ui.btnJ.setChecked(True)
				self.ui.btnVDJ.setChecked(True)

				self.ui.btnCW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_3.setChecked(True)

			# todo need fix this after cdr3 works
			elif self.ui.btnFW_1.isChecked() and self.ui.btnCW_3.isChecked():
				StartSel = int(data[22]) - 1
				EndSel = int(data[85])
				self.ui.btnV.setChecked(True)
				self.ui.btnD.setChecked(True)
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)

			# elif self.ui.btnFW_1.isChecked() and self.ui.btnCW_1.isChecked() and self.ui.btnFW_2.isChecked() and self.ui.btnCW_2.isChecked() and self.ui.btnFW_3.isChecked():
			#     StartSel = int(data[23])
			#     EndSel = int(data[51])
			#     self.ui.btnV.setChecked(True)
			elif self.ui.btnFW_1.isChecked() and self.ui.btnFW_3.isChecked():
				StartSel = int(data[22]) - 1
				EndSel = int(data[51])
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)


			elif self.ui.btnFW_1.isChecked() and self.ui.btnCW_2.isChecked():
				StartSel = int(data[22]) - 1
				EndSel = int(data[44])
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)


			elif self.ui.btnFW_1.isChecked() and self.ui.btnFW_2.isChecked():
				StartSel = int(data[22]) - 1
				EndSel = int(data[37])
				self.ui.btnCW_1.setChecked(True)

			elif self.ui.btnFW_1.isChecked() and self.ui.btnCW_1.isChecked():
				StartSel = int(data[22]) - 1
				EndSel = int(data[30])

			elif self.ui.btnCW_1.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[29]) - 1
				EndSel = int(data[74])
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_3.setChecked(True)

			# todo need fix cdr3 to get this work:
			elif self.ui.btnCW_1.isChecked() and self.ui.btnCW_3.isChecked():
				StartSel = int(data[29]) - 1
				EndSel = int(data[85])
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)

			elif self.ui.btnCW_1.isChecked() and self.ui.btnFW_3.isChecked():
				StartSel = int(data[29]) - 1
				EndSel = int(data[51])
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)

			elif self.ui.btnCW_1.isChecked() and self.ui.btnCW_2.isChecked():
				StartSel = int(data[29]) - 1
				EndSel = int(data[44])
				self.ui.btnFW_2.setChecked(True)

			elif self.ui.btnCW_1.isChecked() and self.ui.btnFW_2.isChecked():
				StartSel = int(data[29]) - 1
				EndSel = int(data[37])


			elif self.ui.btnFW_2.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[36]) - 1
				EndSel = int(data[74])

				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_3.setChecked(True)

			# todo need fix cdr3 to get this work:
			elif self.ui.btnFW_2.isChecked() and self.ui.btnCW_3.isChecked():
				StartSel = int(data[36]) - 1
				EndSel = int(data[85])
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)

			elif self.ui.btnFW_2.isChecked() and self.ui.btnFW_3.isChecked():
				StartSel = int(data[36]) - 1
				EndSel = int(data[51])
				self.ui.btnCW_2.setChecked(True)

			elif self.ui.btnFW_2.isChecked() and self.ui.btnCW_2.isChecked():
				StartSel = int(data[36]) - 1
				EndSel = int(data[44])

			elif self.ui.btnCW_2.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[43]) - 1
				EndSel = int(data[74])
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_3.setChecked(True)

			# todo need fix cdr3 to get this work:
			elif self.ui.btnCW_2.isChecked() and self.ui.btnCW_3.isChecked():
				StartSel = int(data[43]) - 1
				EndSel = int(data[85])
				self.ui.btnFW_3.setChecked(True)



			elif self.ui.btnCW_2.isChecked() and self.ui.btnFW_3.isChecked():
				StartSel = int(data[43]) - 1
				EndSel = int(data[51])

			elif self.ui.btnFW_3.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[50]) - 1
				EndSel = int(data[74])

				self.ui.btnCW_3.setChecked(True)

			# todo need fix cdr3 to get this work:
			elif self.ui.btnFW_3.isChecked() and self.ui.btnCW_3.isChecked():
				StartSel = int(data[50]) - 1
				EndSel = int(data[85])

			elif self.ui.btnCW_3.isChecked() and self.ui.btnFW_4.isChecked():
				StartSel = int(data[84]) - 1
				EndSel = int(data[74])

			else:
				if self.ui.btnFW_1.isChecked():
					StartSel = int(data[22]) - 1
					EndSel = int(data[23])
				if self.ui.btnCW_1.isChecked():
					StartSel = int(data[29]) - 1
					EndSel = int(data[30])
				if self.ui.btnFW_2.isChecked():
					StartSel = int(data[36]) - 1
					EndSel = int(data[37])
				if self.ui.btnCW_2.isChecked():
					StartSel = int(data[43]) - 1
					EndSel = int(data[44])
				if self.ui.btnFW_3.isChecked():
					StartSel = int(data[50]) - 1
					EndSel = int(data[51])
				if self.ui.btnCW_3.isChecked():
					StartSel = int(data[84]) - 1
					EndSel = int(data[85])
				if self.ui.btnFW_4.isChecked():
					StartSel = int(data[85]) - 1
					EndSel = int(data[74])


				else:
					cursor.setPosition(0)
					AAcursor.setPosition(0)

					JustMoved = True
					self.ui.txtDNASeq.setTextCursor(cursor)
					self.ui.txtAASeq.setTextCursor(AAcursor)
					AAcursor.setPosition(0)
					AAcursor.setPosition(len(self.ui.txtAASeq.toPlainText()), QTextCursor.KeepAnchor)
					format = QTextCharFormat()
					format.setFontUnderline(False)
					AAcursor.mergeCharFormat(format)
					AAcursor.setPosition(0)
					JustMoved = False


		elif button == 'c':
			self.ui.btnFW_1.setChecked(False)
			self.ui.btnCW_1.setChecked(False)
			self.ui.btnFW_2.setChecked(False)
			self.ui.btnCW_2.setChecked(False)
			self.ui.btnFW_3.setChecked(False)
			self.ui.btnCW_3.setChecked(False)
			self.ui.btnFW_4.setChecked(False)
			self.ui.btnV.setChecked(False)
			self.ui.btnD.setChecked(False)
			self.ui.btnJ.setChecked(False)
			self.ui.btnVDJ.setChecked(False)
			if int(data[1]) > int(data[74]):
				StartSel = int(data[74]) + 1
				EndSel = len(data[79])
		# todo need code foe canstant gene for this button to finish

		elif button == 'vdj':
			if self.ui.btnVDJ.isChecked():
				self.ui.btnFW_1.setChecked(True)
				self.ui.btnCW_1.setChecked(True)
				self.ui.btnFW_2.setChecked(True)
				self.ui.btnCW_2.setChecked(True)
				self.ui.btnFW_3.setChecked(True)
				self.ui.btnCW_3.setChecked(True)
				self.ui.btnFW_4.setChecked(True)
				self.ui.btnV.setChecked(True)
				self.ui.btnD.setChecked(True)
				self.ui.btnJ.setChecked(True)
				self.ui.btnC.setChecked(False)
				StartSel = int(data[67]) - 1
				EndSel = int(data[74])
			else:
				self.ui.btnFW_1.setChecked(False)
				self.ui.btnCW_1.setChecked(False)
				self.ui.btnFW_2.setChecked(False)
				self.ui.btnCW_2.setChecked(False)
				self.ui.btnFW_3.setChecked(False)
				self.ui.btnCW_3.setChecked(False)
				self.ui.btnFW_4.setChecked(False)
				self.ui.btnV.setChecked(False)
				self.ui.btnD.setChecked(False)
				self.ui.btnJ.setChecked(False)
				self.ui.btnC.setChecked(False)
				cursor.setPosition(0)
				AAcursor.setPosition(0)
				self.ui.txtDNASeq.setTextCursor(cursor)
				self.ui.txtAASeq.setTextCursor(AAcursor)
				AAcursor.setPosition(0)
				AAcursor.setPosition(len(self.ui.txtAASeq.toPlainText()), QTextCursor.KeepAnchor)
				format = QTextCharFormat()
				format.setFontUnderline(False)
				AAcursor.mergeCharFormat(format)
				AAcursor.setPosition(0)

		JustMoved = True
		cursor.setPosition(StartSel)
		cursor.setPosition(EndSel, QTextCursor.KeepAnchor)
		if StartSel != 0:
			AAStartSel = StartSel / 3
		else:
			AAStartSel = StartSel
		if EndSel != 0:
			AAEndSel = EndSel / 3
		else:
			AAEndSel = EndSel

		AAcursor.setPosition(AAStartSel)
		AAcursor.setPosition(AAEndSel, QTextCursor.KeepAnchor)

		self.ui.txtDNASeq.setTextCursor(cursor)
		if self.ui.cboDecorate.currentText() == 'None':
			self.ui.txtAASeq.setTextCursor(AAcursor)
		else:
			# if AAEndSel>0:
			self.ui.txtAASeq.setTextCursor(AAcursor)
			format = QTextCharFormat()
			format.setFontUnderline(True)
			format.setUnderlineColor(QColor("cyan"))
			format.setUnderlineStyle(QTextCharFormat.WaveUnderline)  # QTextCharFormat('WaveUnderline')
			AAcursor.mergeCharFormat(format)
			AAcursor.setPosition(0)
			self.ui.txtAASeq.setTextCursor(AAcursor)

		JustMoved = False

	@pyqtSlot()
	def get_text_selection(self):
		cursor = self.ui.txtDNASeq.textCursor()
		DNAseq = self.ui.txtDNASeq.toPlainText()
		lenSeq = len(DNAseq)
		return cursor.selectionStart(), cursor.selectionEnd(), lenSeq

	@pyqtSlot()
	def on_btnSearchSeq_clicked(self):
		search = self.ui.txtSearchSeq.toPlainText()
		self.FindSeq(search)

	@pyqtSlot()
	def FindSeq(self, SeqFind):
		FindSeq = SeqFind

		Found = self.ui.txtDNASeq.find(FindSeq)

		if Found == False:
			self.SeqButton('none')
			Found = self.ui.txtAASeq.find(FindSeq)

		if Found == False:
			msg = SeqFind + ' could not be found.'
			buttons = 'OK'
			answer = informationMessage(self, msg, buttons)

		# self.UpdatePositions()

		return Found

	@pyqtSlot()
	def on_txtDNASeq_textChanged(self):
		DNAseq = self.ui.txtDNASeq.toPlainText()
		DataIs = []
		global GLMsg
		global JustMoved
		if self.ui.radioButtonGermView.isChecked():
			return
		if self.ui.btnEditSeq.isChecked():


			if self.ui.radioButtonGermView.isChecked():
				self.on_radioButtonGermView_clicked()




			# self.ui.lblSeq2.setText(msg)
			frame = 0
			ErMes = 'Amino acid sequence: \n\n'



			AASeq, ErMessage = VGenesSeq.Translator(DNAseq, frame)
			if len(ErMessage) > 0:
				for mess in ErMessage:
					ErMes += mess + '\n'

			self.ui.lblAAseq.setText(ErMes)
			self.ui.txtAASeq.setPlainText(AASeq)

			if self.ui.cboDecorate.currentText() != "None":
				self.DecoratePeptide()

			self.UpdatePositions()

			self.clearTreeChecks()

			self.AlignSequences('edit')


			self.ui.txtDNASeq.setFocus()

		else:
			if JustMoved == False:
				msg = 'Enter "edit mode" to edit sequence and related fields (mutations, regions, etc.)?'

				buttons = 'OKC'
				answer = questionMessage(self, msg, buttons)
				if answer == 'OK':
					self.ui.btnEditSeq.setChecked(True)
					self.on_btnEditSeq_clicked()

					JustMoved = True

		# self.ui.lblSeq2.setText(msg)
		frame = 0
		ErMes = 'Amino acid sequence: \n\n'

		AASeq, ErMessage = VGenesSeq.Translator(DNAseq, frame)
		if len(ErMessage) > 0:
			for mess in ErMessage:
				ErMes += mess + '\n'

		self.ui.lblAAseq.setText(ErMes)
		self.ui.txtAASeq.setPlainText(AASeq)

		if self.ui.cboDecorate.currentText() != "None":
			self.DecoratePeptide()

		self.UpdatePositions()


	def get_Align_selection(self):
		cursor = self.ui.txtSeqAlignment.textCursor()
		DNAseq = self.ui.txtSeqAlignment.toPlainText()
		lenSeq = len(DNAseq)
		return cursor.selectionStart(), cursor.selectionEnd(), lenSeq

	def on_txtSeqAlignment_cursorPositionChanged(self):
		StartP, EndP, LenSeq = self.get_Align_selection()

		if StartP == EndP:
			lblText = 'Ig BLAST Alignment: position = ' + str(EndP)
		else:
			lblText = 'Ig BLAST Alignment: position = ' + str(StartP) + ' to ' + str(EndP) + ' equaling ' + str(
				EndP - StartP) + ' total nucleotides'

		self.ui.lblAlignment.setText(lblText)

	def on_txtDNASeq_cursorPositionChanged(self):
		# DNAseq = self.ui.txtDNASeq.toPlainText()
		try:
			if JustMoved == False:
				self.ui.btnFW_1.setChecked(False)
				self.ui.btnCW_1.setChecked(False)
				self.ui.btnFW_2.setChecked(False)
				self.ui.btnCW_2.setChecked(False)
				self.ui.btnFW_3.setChecked(False)
				self.ui.btnCW_3.setChecked(False)
				self.ui.btnFW_4.setChecked(False)
				self.ui.btnV.setChecked(False)
				self.ui.btnD.setChecked(False)
				self.ui.btnJ.setChecked(False)
				self.ui.btnC.setChecked(False)
		except:
			# todo need fix this waste
			print('tried4384')

		self.UpdatePositions()

	def UpdatePositions(self):
		StartP, EndP, LenSeq = self.get_text_selection()

		cursor = self.ui.txtDNASeq.textCursor()
		if StartP == EndP:
			lblText = 'DNA Sequence: position = ' + str(EndP) + ' of ' + str(LenSeq) + ' total nucleotides'
		else:
			lblText = 'DNA Sequence: ' + str(StartP + 1) + ' to ' + str(EndP + 1) + ' (' + str(
				EndP - StartP) + ' bases) ' ' selected of ' + str(LenSeq) + ' total nucleotides'

		self.ui.lblDNASeq.setText(lblText)

	# @pyqtSlot()
	# def on_btnEditSeq_clicked(self):

	@pyqtSlot()
	def on_btnEditSeq_clicked(self):
		if self.ui.btnEditSeq.isChecked():
			self.ui.btnEditSeq.setText("Save Changes")
			msg = 'Press button again to save changes'

			self.ui.lblSeq2.setText(msg)


			if self.ui.radioButtonGermView.isChecked():
				self.on_radioButtonGermView_clicked()


			# GLMsg = False
			# self.on_actionGL_triggered
			# GLMsg = True
			self.clearTreeChecks()
			self.AlignSequences('none')

			msg = 'Results of edits are aligned to the predicted germline in the text editor window. Changes will not be saved until Save Changes is pressed.'
			buttons = 'OK'
			answer = informationMessage(self, msg, buttons)

		else:

			msg = 'Saving changes to the sequence will cause it to be reanlyzed (mutations, regions) but will retain other information (clone, isotype). Continue?'
			buttons = 'OKC'
			answer = questionMessage(self, msg, buttons)
			if answer == 'OK':

				self.UpdateSeqAnalysis()


			self.ui.btnEditSeq.setText("Edit Mode")
			msg = 'Prese "Edit Mode" to edit sequence.'

			self.ui.lblSeq2.setText(msg)



			# todo code to update need update record function in VGenes SQL:

	def UpdateSeqAnalysis(self):
		datalist = []
		currentRow = self.ui.tableView.currentIndex().row()
		# todo change to app folder
		try:
			filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'UpdateRecord.nt')
			# filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'UpdateRecord.nt')
			with open(filename, 'r') as currentfile:
				vv = currentfile
		except:
			filename = '/Volumes/Promise Pegasus/Dropbox/VGenes/UpdateRecord.nt'

		sequence = self.ui.txtDNASeq.toPlainText()
		seqname = data[0]
		updateseq = '>' + seqname + '\n' + sequence + '\n'

		with open(filename, 'w') as currentfile:
			currentfile.write(updateseq)

		if filename == None:
			return
		project = self.ui.txtProject.toPlainText()
		group = self.ui.txtGroup.toPlainText()
		subgroup = self.ui.txtSubGroup.toPlainText()
		species = data[78]
		project = project.strip()
		group = group.strip()
		subgroup = subgroup.strip()
		species = species.strip()

		datalist.append(project)
		datalist.append(group)
		datalist.append(subgroup)
		datalist.append(species)
		datalist.append(False)
		# BlastIsDone = False

		IgBLASTAnalysis = IgBLASTer.IgBLASTit(filename, datalist)

		# while BlastIsDone == False:
		#     NotDone  = False
		# todo need to put any other preserved data into IgBLASTAnalysis
		i = 0
		model = self.ui.tableView.model()
		DataRow = int(data[119])

		for record in IgBLASTAnalysis:
			for item in record:
				FieldName = FieldList[i]
				ItemValue = str(item) #+ ' '
				if i != 119:
					VGenesSQL.UpdateField(DataRow, ItemValue, FieldName, DBFilename)


				i += 1


		model.refresh()

		global PreVID
		# NewID = int(data[119])+1
		# PreVID = NewID

		self.updateF(data[119])






		# Processed = VGenesSQL.updateData(self, DBFilename, IgBLASTAnalysis)


# todo need code to allow sequence to be edited then when save pushed record reanlyzed across the board

class EditableSqlModel(QSqlQueryModel):
	def flags(self, index):
		flags = super(EditableSqlModel, self).flags(index)

		if index.column() in (
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
		30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
		58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86,
		87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
		112, 113, 114, 115, 116, 117, 118):
			flags |= Qt.ItemIsEditable
		# print(str(index.row()))
		ID = Vgenes.ui.tableView.currentIndex().row()
		NameIs = ''
		if ID >= 0 and ID != PreVID:  # made this also not equal preVID before firing
			NameIs = Vgenes.MatchingValue(ID)
			Vgenes.findTreeItem(NameIs)
			ReportName = Vgenes.ui.txtName.toPlainText()  # tabbed this and next 2 lines in
			if ReportName != NameIs and ID >= 0:
				Vgenes.updateF(ID)

		return flags

	def setData(self, index, value, role):
		if index.column() not in (
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
		30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
		58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86,
		87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
		112, 113, 114, 115, 116, 117, 118):
			return False

		recordIs = Vgenes.ui.treeWidget.selectedItems()
		currentitemIs = ''

		for item in recordIs:
			currentitemIs = item.text(0)

		primaryKeyIndex = self.index(index.row(), 119)
		id = self.data(primaryKeyIndex)

		self.clear()

		i = 0
		for item in FieldList:
			if index.column() == i:
				field = item
				ok = self.setValue(id, value, field)
				self.refresh()
				if field == Vgenes.ui.cboTreeOp1.setCurrentText or field == Vgenes.ui.cboTreeOp2.setCurrentText or field == Vgenes.ui.cboTreeOp3.setCurrentText:
					Vgenes.on_btnUpdateTree_clicked()

				Vgenes.findTreeItem(currentitemIs)
				Vgenes.ui.tableView.setCurrentIndex(index)
				return
			i += 1

		self.refresh()

		RowID = Vgenes.ui.tableView.currentIndex()
		model = Vgenes.ui.tableView.model()
		# for i in range(0, 83):
		#     index = model.index(RowID, i)
		#     data.append (str(model.data(index)))


		return

	def refresh(self):

		self.setQuery(RefreshSQL)
		i = 0
		for item in FieldList:
			self.setHeaderData(i, Qt.Horizontal, item)
			i += 1

	def setValue(self, Id, Newvalue, field):
		query = QSqlQuery()
		queryText = 'update vgenesdb set ' + field + ' = ? where ID = ?'
		query.prepare(queryText)
		query.addBindValue(Newvalue)
		query.addBindValue(Id)
		return query.exec_()

		# @pyqtSlot("QModelIndex")
		# def ItemClicked(self,index):
		#     # QMessageBox.information(None,"Hello!","You Clicked: \n"+index.data().toString())
		#     # print(index.data().toString())
		#     ID = self.model.record(index).value('ID')
		#     SeqName = self.model.record(index).value('SeqName')
		#     print(ID)
		#     print(SeqName)


# todo need code to properly navigate tableView with arrow keys and tabs

if __name__ == '__main__':
	import sys

	app = QtWidgets.QApplication(sys.argv)
	Vgenes = VGenesForm()

	Vgenes.ApplicationStarted()
	sys.exit(app.exec_())
