__author__ = 'wilsonp'
import os

import sqlite3 as db
import time

from PyQt5.QtCore import pyqtSlot, QTimer, QDateTime, Qt, QSortFilterProxyModel, QModelIndex, pyqtSignal, QUrl
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

global UpdateSpecific
UpdateSpecific = True
import csv
from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
	setText, openfastq
from ui_VGenesStartUpDialogue import Ui_VGenesStartUpDialog

from ui_VGenesTextEdit import ui_TextEditor
from VGenesProgressBar import ui_ProgressBar
# from VGenesPYQTSqL import EditableSqlModel, initializeModel , createConnection

from PyQt5.QtCore import Qt, QObject, QEvent
from PyQt5.QtGui import QTextCursor, QFont, QPixmap, QTextCharFormat, QBrush, QColor, QTextCursor, QCursor
from PyQt5.QtWidgets import QApplication, QTableView
from PyQt5.QtSql import QSqlQuery, QSqlQueryModel
from operator import itemgetter
from PyQt5.QtWebEngine import *
from PyQt5.QtWebEngineWidgets import *
from pyecharts.charts import *
from pyecharts import options as opts
from pyecharts.globals import SymbolType
import itertools

from itertools import combinations
from collections import Counter
from subprocess import call, Popen, PIPE
from platform import system

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
global PreviewHTExp
PreviewHTExp = []
global PreviewHTcurrent
PreviewHTcurrent = 0
global PreviewCurrentType
PreviewCurrentType = 'H'
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
                "Germlne CDR3 begin", "Germline CDR3 end", "Autoreactivity", "Blank7", "10xCluster", "Seuret_Cluster", "10xBarCode", "Population",
                "Label", "Status", "Blank14", "Blank15", "Blank16", "Blank17", "Blank18", "Blank19", "Blank20", "ID"]
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

global GLMsg
GLMsg = True

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
			filename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'RecentPaths.vtx')

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
			filename = 'RecentPaths.vtx'

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
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			if len(DataIn) > 0:
				for item in DataIn:
					self.comboBoxProject.addItem(item[0])
			DataIn.clear()

			fields = ['Grouping']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if len(DataIn) > 0:
				for item in DataIn:
					self.comboBoxGroup.addItem(item[0])
			DataIn.clear()

			fields = ['SubGroup']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

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
		MaxNum = self.MaxImport.value()
		Alldone = self.InitiateImport('none', MaxNum)

		if Alldone == True:
			self.close()

	@pyqtSlot()
	def on_buttonBox_rejected(self):

		self.close()

	# def timingcheck(self):
	# 	# TODO can delete this if gets work
	#
	# 	# import IgBLASTer
	#
	# 	# stmtS = "/Users/PCW-MacBookProRet/Dropbox/VGenes/Database/SFV-005H.nt"
	# 	DirTo = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes')
	# 	import timeit, os
	# 	# todo change to app folder
	# 	try:
	# 		os.chdir(DirTo)
	# 	except:
	# 		os.chdir('/Volumes/Promise Pegasus/Dropbox/VGenes/VGenes')
	#
	# 	t = timeit.timeit('self.findTableViewRecord(self, FieldName)")', 'import IgBLASTer', number=10000)
	# 	print(t)

	def InitiateImport(self, Filenamed, MaxNum):
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
				pathname = openFiles(self, thetype)
			else:
				pathname = Filenamed[0]
			if pathname == None:
				return False


		elif thetype == 'Sequence':
			filenames = openFiles(self, 'seq')
			pathname1 = self.ProcessSeqFiles(filenames)
			pathname  = []
			pathname.append(pathname1)
			# filename = filenames[0]

		if self.rdoProductive.isChecked() == True:
			GetProductive = True
		else:
			GetProductive = False

		if pathname == None:
			return
		answer2 = ''

		ErlogFile = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'IgBlast', 'database',
		                         'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs


		ErlogFile2 = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'IgBlast', 'database',
		                         'ErLog2.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
		header = "Began input at " + time.strftime('%c')
		with open(ErlogFile2, 'w') as currentFile:
			currentFile.write(header)
		# firstOne = True

		if self.checkBoxFileStruc.isChecked():


			for item in pathname:

				(dirname, filename) = os.path.split(item)
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

					# if thetype != 'Sequence':
					if answer2 == '':
						if len(pathname) > 1:
							msg = 'More then 1 FASTA file was selected. Make each a seperate project based on the filenames?'
							buttons = 'YN'

							answer2 = informationMessage(self, msg, buttons)

							firstOne = False

				if answer2 == 'Yes':
					preproject = os.path.splitext(filename)
					project = preproject[0]


				datalist.clear()
				datalist.append(project)
				datalist.append(grouping)
				datalist.append(subgroup)
				datalist.append(species)
				datalist.append(GetProductive)
				datalist.append(MaxNum)




				IgBLASTAnalysis = IgBLASTer.IgBLASTit(item, datalist)

				Startprocessed = len(IgBLASTAnalysis)
				if Startprocessed == 0:
					self.close()

				Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)

				i = 0
				newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: ' + item + '\n'


				with open(ErlogFile,
				          'r') as currentFile:  # using with for this automatically closes the file even if you crash
					for line in currentFile:
						if i > 0:
							newErLog += line
						i += 1


				with open(ErlogFile2, 'a') as currentFile:
					currentFile.write(newErLog)






		elif self.rdoChoose.isChecked():
			# checklabel = {}

			for item in pathname:
				(dirname, filename) = os.path.split(item)


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

				# if thetype != 'Sequence':
				if answer2 == '':
					if len(pathname) > 1:
						msg = 'More then 1 FASTA file was selected. Make each a seperate project based on the filenames?'
						buttons = 'YN'

						answer2 = informationMessage(self, msg, buttons)

						firstOne = False

				if answer2 == 'Yes':
					preproject = os.path.splitext(filename)
					project = preproject[0]


				datalist.clear()

				datalist.append(project)
				datalist.append(grouping)
				datalist.append(subgroup)
				datalist.append(species)
				datalist.append(GetProductive)
				datalist.append(MaxNum)



				IgBLASTAnalysis = IgBLASTer.IgBLASTit(item, datalist)

				Startprocessed = len(IgBLASTAnalysis)
				if Startprocessed == 0:
					self.close()

				Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)

				i = 0
				newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: ' + item + '\n'


				with open(ErlogFile,
				          'r') as currentFile:  # using with for this automatically closes the file even if you crash
					for line in currentFile:
						if i > 0:
							newErLog += line
						i += 1


				with open(ErlogFile2, 'a') as currentFile:
					currentFile.write(newErLog)




		elif self.rdoFunction.isChecked():

			for item in pathname:
				(dirname, filename) = os.path.split(item)

				project = 'ByFunction'
				grouping = ''
				subgroup = ''

				# if thetype != 'Sequence':
				if answer2 == '':
					if len(pathname) > 1:
						msg = 'More then 1 FASTA file was selected. Make each a seperate project based on the filenames?'
						buttons = 'YN'

						answer2 = informationMessage(self, msg, buttons)

						firstOne = False

				if answer2 == 'Yes':
					preproject = os.path.splitext(filename)
					multiProject = preproject[0]
				else:
					multiProject = ''


				datalist.clear()



				datalist.append(project)
				datalist.append(grouping)
				datalist.append(subgroup)
				datalist.append(species)
				datalist.append(GetProductive)
				datalist.append(MaxNum)
				datalist.append(multiProject)

				# if thetype == 'Sequence':
				# 	item = pathname
				IgBLASTAnalysis = IgBLASTer.IgBLASTit(item, datalist)
				Startprocessed = 0
				try:

					Startprocessed = len(IgBLASTAnalysis)
				except:
					if Startprocessed == 0:
						self.close()

				Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)

				i = 0
				newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: ' + item + '\n'

				with open(ErlogFile,
				          'r') as currentFile:  # using with for this automatically closes the file even if you crash
					for line in currentFile:
						if i > 0:
							newErLog += line
						i += 1

				with open(ErlogFile2, 'a') as currentFile:
					currentFile.write(newErLog)


		Vgenes.LoadDB(DBFilename)

		self.ShowVGenesText(ErlogFile2)



	@pyqtSlot()
	def on_btnImportOldVGenes_clicked(self):
		from operator import itemgetter   #		SeqList.sort(key=itemgetter(0, 1, 2, 3))
		msg = 'This function imports a comma separated values (CSV) file formatted as: Project, Group, Subgroup, Name, Sequence'
		buttons = 'OKC'
		global answer3
		answer3 = 'No'
		answer = informationMessage(self, msg, buttons)
		if answer == 'Cancel':
			print('no file')
			return

		# self.rdoChoose.setChecked(True)
		# self.checkBoxFileStruc.setChecked(False)
		self.rdoFunction.setChecked(False)


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
		ItemP = []
		FileNamed = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'IgBlast', 'database',
		                         'WorkingFile.nt')  # '/Applications/IgBlast/database/WorkingFile.nt'
		ItemP.append(FileNamed)
		DataPass.append(ItemP)
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

					Alldone = self.InitiateImport(DataPass, 0)
					DataPass[1] = project
					# if project == 'STAN-004':
					# 	print('stop')
					DataPass[2] = group
					DataPass[3] = subgroup

					self.comboBoxProject.setCurrentText(project)
					self.comboBoxGroup.setCurrentText(group)
					self.comboBoxSubgroup.setCurrentText(subgroup)

				else:
					DataPass[1] = project
					DataPass[2] = group
					DataPass[3] = subgroup
					self.comboBoxProject.setCurrentText(project)
					self.comboBoxGroup.setCurrentText(group)
					self.comboBoxSubgroup.setCurrentText(subgroup)


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

class ResizeWidget(QWebEngineView):
	resizeSignal = pyqtSignal(int,int)
	def __init__(self,parent=None):
		super(ResizeWidget,self).__init__()
		self.id = 0
		self.h = 0
		self.w = 0

	def resizeEvent(self, evt):
		#w = evt.oldSize().width()
		#h = evt.oldSize().height()
		#print(f' size before :{w,h}')
		w = evt.size().width()
		h = evt.size().height()
		self.h = h
		self.w = w
		print(f' size now :{w, h, self.id}')
		self.resizeSignal.emit(w,h)

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
		self.ui.cboFindField.currentTextChanged.connect(self.on_cboFindField_currentTextChanged)
		self.ui.tabWidget.currentChanged['int'].connect(self.InitialGraphic)
		self.ui.toolButton.clicked.connect(self.GenerateFigure)
		self.ui.checkBoxFigLegend.clicked.connect(self.GenerateFigure)
		# self.ui.listViewSpecificity.highlighted['QString'].connect(self.SpecSet)
		# self.ui.listViewSpecificity.mouseDoubleClickEvent.connect(self.SpecSet)

		self.ui.lcdNumber_max.display(self.ui.horizontalScrollBar.maximum())
		self.ui.dial.setMaximum(self.ui.horizontalScrollBar.maximum())

		self.TextEdit = VGenesTextMain()
		# self.VGProgress = VGenesProgressBar()
		self.ImportOptions = ImportDialogue()

		self.ui.HTMLview = ResizeWidget(self)
		self.ui.gridLayoutStat.addWidget(self.ui.HTMLview, 2, 0, 10, 0)
		self.ui.HTMLview.resizeSignal.connect(self.resizeHTML)

	# def ShowProgressBar(self):
	#     self.VGProgress.handleButton()


	# @pyqtSlot()
	# def on_radioButton_21_clicked(self):
	# 	self.PopulateSpec()

	def resizeHTML(self, w, h):
		global DBFilename

		h = h - 20
		sender = self.sender()
		a = sender.id
		if sender.id == 1:
			pass

	def InitialGraphic(self):
		global DBFilename
		#Msg = str(self.ui.tabWidget.currentIndex())
		#QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
		if self.ui.tabWidget.currentIndex() == 6:
			fields_name = VGenesSQL.ColName(DBFilename)
			fields_name = [""] + fields_name
			self.ui.comboBoxPie.addItems(fields_name)
			self.ui.comboBoxCol1.addItems(fields_name)
			self.ui.comboBoxCol2.addItems(fields_name)
			self.ui.comboBoxBoxData.addItems(fields_name)
			self.ui.comboBoxBox1.addItems(fields_name)
			self.ui.comboBoxBox2.addItems(fields_name)
			self.ui.comboBoxRiver1.addItems(fields_name)
			self.ui.comboBoxRiver2.addItems(fields_name)
			self.ui.comboBoxWord.addItems(fields_name)

	def GenerateFigure(self):
		global DBFilename

		if self.ui.tabWidgetFig.currentIndex() == 0:
			# get data
			field = self.ui.comboBoxPie.currentText()
			if field == "":
				return
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			data = []
			for element in DataIn:
				data.append(element[0])
			result = Counter(data)
			labels = result.keys()
			values = result.values()

			my_pyecharts = (
				Pie(init_opts=opts.InitOpts(width="380px", height="380px"))
				.add('', [list(z) for z in zip(labels, values)], radius=["40%", "75%"])
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show = self.ui.checkBoxFigLegend.isChecked()
					),
				)
				.set_series_opts(label_opts=opts.LabelOpts(formatter=" {b}: {c} ({d}%)"))
			)
		elif self.ui.tabWidgetFig.currentIndex() == 1:
			# get data
			field1 = self.ui.comboBoxCol1.currentText()
			field2 = self.ui.comboBoxCol2.currentText()
			if field1 == "":
				return

			multi_factor = False
			if field2 == "":
				field = field1
			else:
				field = field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if multi_factor == True:
				if self.ui.checkBoxStack.isChecked():
					stack = "stack1"
				else:
					stack = None

				label_data = []
				for element in DataIn:
					label_data.append(element[0])

				result = Counter(label_data)
				labels = list(result.keys())
				values = list(result.values())

				data = {}
				for element in DataIn:
					if data.__contains__(element[1]):
						data[element[1]] = data[element[1]] + [element[0]]
					else:
						data[element[1]] = [element[0]]

				dic_keys = list(data.keys())

				my_bar = Bar(init_opts=opts.InitOpts(width="380px", height="380px"))\
					.add_xaxis(labels)\
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegend.isChecked()
						),
					)\
					.set_series_opts(label_opts=opts.LabelOpts(formatter=" {b}: {c}"))

				for group in dic_keys:
					cur_data = data[group]
					group_data = []
					for ele in labels:
						group_data.append(cur_data.count(ele))
					my_bar.add_yaxis(group, group_data,stack=stack)

				my_pyecharts = (
					my_bar
				)
			else:
				data = []
				for element in DataIn:
					data.append(element[0])

				result = Counter(data)
				labels = list(result.keys())
				values = list(result.values())

				my_pyecharts = (
					Bar(init_opts=opts.InitOpts(width="380px", height="380px"))
						.add_xaxis(labels)
						.add_yaxis(field1, values)
						.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegend.isChecked()
						),
					)
					.set_series_opts(label_opts=opts.LabelOpts(formatter=" {b}: {c}"))
				)
		elif self.ui.tabWidgetFig.currentIndex() == 2:
			# get data
			data_field = self.ui.comboBoxBoxData.currentText()
			field1 = self.ui.comboBoxBox1.currentText()
			field2 = self.ui.comboBoxBox2.currentText()
			if field1 == "":
				return
			multi_factor = False
			if field2 == "":
				field = data_field + "," + field1
			else:
				field = data_field + "," + field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			box_data = [i[0] for i in DataIn]
			box_data = list(map(float, box_data))

			if min(box_data) >= 0:
				null_data = [0,0,0,0,0]
			else:
				null_data = [min(box_data), min(box_data), min(box_data), min(box_data), min(box_data)]

			if multi_factor == True:
				label_data = []
				g2_label = []
				for element in DataIn:
					label_data.append(element[1])
					g2_label.append(element[2])

				result = Counter(label_data)
				labels = list(result.keys())

				g2_dict = {}
				i = 0
				for ele in g2_label:
					if g2_dict.__contains__(ele):
						g2_dict[ele] = g2_dict[ele] + [i]
					else:
						g2_dict[ele] = [i]
					i += 1
				g2_dict_keys = list(g2_dict.keys())

				my_bar = Boxplot(init_opts=opts.InitOpts(width="380px", height="380px"))\
					.add_xaxis(labels)\
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegend.isChecked()
						),
						#toolbox_opts = opts.ToolboxOpts()
					)

				# for each group in field 2
				for group in g2_dict_keys:
					cur_box_data = []
					cur_label_data = []
					for i in g2_dict[group]:
						cur_box_data.append(box_data[i])
						cur_label_data.append(label_data[i])

					sub_dict = {}
					i = 0
					for ele in cur_label_data:
						if sub_dict.__contains__(ele):
							sub_dict[ele] = sub_dict[ele] + [i]
						else:
							sub_dict[ele] = [i]
						i += 1

					data_v1 = []
					for ele in labels:
						if sub_dict.__contains__(ele):
							cur_data = []
							for i in sub_dict[ele]:
								cur_data.append(cur_box_data[i])
							data_v1.append(cur_data)
						else:
							data_v1.append(null_data)

					my_bar.add_yaxis(group, Boxplot.prepare_data(data_v1))

				my_pyecharts = (
					my_bar
				)
			else:
				data = []
				for element in DataIn:
					data.append(element[1])

				result = Counter(data)
				labels = list(result.keys())

				my_dict = {}
				i = 0
				for ele in data:
					if my_dict.__contains__(ele):
						my_dict[ele] = my_dict[ele] + [i]
					else:
						my_dict[ele] = [i]
					i += 1

				data_v1 = []
				for ele in labels:
					if my_dict.__contains__(ele):
						cur_data  = []
						for i in my_dict[ele]:
							cur_data.append(box_data[i])
						data_v1.append(cur_data)
					else:
						data_v1.append(null_data)

				my_pyecharts = (
					Boxplot(init_opts=opts.InitOpts(width="380px", height="380px"))
						.add_xaxis(labels)
						.add_yaxis(field1, Boxplot.prepare_data(data_v1))
						.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegend.isChecked()
						),
					)
				)
		elif self.ui.tabWidgetFig.currentIndex() == 3:
			# get data
			field = self.ui.comboBoxWord.currentText()
			if field == "":
				return
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			data = []
			for element in DataIn:
				data.append(element[0])
			result = Counter(data)
			keys = list(result)
			word_data = []
			for ele in keys:
				unit = (ele,result[ele])
				word_data.append(unit)
			my_pyecharts = (
				WordCloud(init_opts=opts.InitOpts(width="380px", height="380px"))
				.add("", word_data, word_size_range=[40, 200], shape=SymbolType.DIAMOND)
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show=self.ui.checkBoxFigLegend.isChecked()
					),
				)
			)

			# load figure
			html_path = '/Users/leil/Documents/Projects/VGenes/test.html'
			my_pyecharts.render(path=html_path)
			# adjust the window size seting
			file_handle = open(html_path, 'r')
			lines = file_handle.readlines()
			file_handle.close()
			# edit js line
			js_line = '<script type="text/javascript" src="' + \
			          os.path.join('echarts.min.js') + '"></script>\n' + \
			          '<script type="text/javascript" src="echarts-wordcloud.min.js"></script>'
			lines[5] = js_line
			# edit style line
			style_line = lines[10]
			style_pos = style_line.find('style')
			style_line = style_line[
			             0:style_pos] + 'style="position: fixed; top: 0px; left: 5%;width:90%; height:' + str(
				self.ui.HTMLview.h - 20) + 'px;"></div>'
			lines[9] = style_line
			content = '\n'.join(lines)
			file_handle = open(html_path, 'w')
			file_handle.write(content)
			file_handle.close()
			# show local HTML
			self.ui.HTMLview.load(QUrl('file://' + html_path))
			self.ui.HTMLview.show()
			return
		elif self.ui.tabWidgetFig.currentIndex() == 4:
			# get data
			field1 = self.ui.comboBoxRiver1.currentText()
			field2 = self.ui.comboBoxRiver2.currentText()
			if field1 == "" or field2 == "":
				return

			field = field1 + "," + field2
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			label_data = []
			time_data = []
			for element in DataIn:
				label_data.append(element[1])
				time_data.append(element[0])
			time_data = list(map(float, time_data))
			result = Counter(label_data)
			labels = list(result.keys())

			min_data = min(time_data)
			max_data = max(time_data)
			data_tick = 20
			step = (max_data - min_data)/data_tick
			
			data_river = []
			for i in range(0, 20):
				cur_tick = min_data + step * i
				cur_tick_end = min_data + step * (i + 1)
				data_in_range = []
				for element in DataIn:
					cur_vbal = float(element[0])
					if cur_vbal <= cur_tick_end and cur_vbal >= cur_tick:
						data_in_range.append(element[1])
				res = Counter(data_in_range)
				for ele in labels:
					if res.__contains__(ele):
						unit = [int(cur_tick), res[ele], ele]
					else:
						unit = [int(cur_tick), 0, ele]
					data_river.append(unit)

			my_pyecharts = (
				ThemeRiver()
				.add(
					labels,
					data_river,
					singleaxis_opts=opts.SingleAxisOpts(type_='value', min_='dataMin', max_='dataMax'),
				)
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show=self.ui.checkBoxFigLegend.isChecked()
					),
					#toolbox_opts=opts.ToolboxOpts()
				)
			)
		elif self.ui.tabWidgetFig.currentIndex() == 5:
			pass
		elif self.ui.tabWidgetFig.currentIndex() == 6:
			pass

		# load figure
		html_path = '/Users/leil/Documents/Projects/VGenes/test.html'
		my_pyecharts.render(path=html_path)
		# adjust the window size seting
		file_handle = open(html_path, 'r')
		lines = file_handle.readlines()
		file_handle.close()
		# edit js line
		js_line = '            <script type="text/javascript" src="' + \
		          os.path.join('echarts.min.js') + '"></script>'
		lines[5] = js_line
		# edit style line
		style_line = lines[9]
		style_pos = style_line.find('style')
		style_line = style_line[
		             0:style_pos] + 'style="position: fixed; top: 0px; left: 5%;width:90%; height:' + str(
			self.ui.HTMLview.h - 20) + 'px;"></div>'
		lines[9] = style_line
		content = '\n'.join(lines)
		file_handle = open(html_path, 'w')
		file_handle.write(content)
		file_handle.close()
		# show local HTML
		self.ui.HTMLview.load(QUrl('file://' + html_path))
		self.ui.HTMLview.show()


	def PopulateSpec(self):
		# if self.ui.radioButton_21.isChecked():
		self.ui.listViewSpecificity.setEditable(True)
		self.ui.listViewSpecificity_2.setEditable(True)



		SQLStatement = 'SELECT DISTINCT Specificity FROM vgenesDB'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		self.ui.listViewSpecificity.clear()


		entry  = ''
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		if len(DataIn) > 0:
			for item in DataIn:
				self.ui.listViewSpecificity.addItem(item[0])
		DataIn.clear()


		SQLStatement = 'SELECT DISTINCT Subspecificity FROM vgenesDB'

		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		self.ui.listViewSpecificity_2.clear()


		entry  = ''
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		if len(DataIn) > 0:
			for item in DataIn:
				self.ui.listViewSpecificity_2.addItem(item[0])
		DataIn.clear()
		global UpdateSpecific
		UpdateSpecific = False






	def MakeSQLStatement(self, fields):
		checkedProjects, checkedGroups, checkedSubGroups, checkedkids = self.getTreeChecked()

		SQLStatement = 'SELECT '

		if fields != 'All':
			fieldCount = len(fields)
			i = 1
			for field in fields:
				SQLStatement += field
				if i < fieldCount:
					SQLStatement += ', '
				else:
					SQLStatement += ' FROM vgenesDB'
				i += 1
		else:
			SQLStatement += '* FROM vgenesDB'

		firstmore = False

		if (len(checkedProjects) + len(checkedGroups) + len(checkedSubGroups) + len(
				checkedkids)) > 0:  # then somehting is seleected
			SQLStatement += ' WHERE '
			firstmore = True

		i = 1
		if len(checkedProjects) > 0:
			if firstmore == True:
				firstmore = False
			else:
				SQLStatement += ', '
				firstmore = False
			project = self.ui.cboTreeOp1.currentText()
			fieldname = self.TransLateFieldtoReal(project, True)
			SQLStatement = SQLStatement + fieldname + ' = '
			for item in checkedProjects:
				SQLStatement += ('"' + item)
				if i < len(checkedProjects):
					SQLStatement += '" OR '
				else:
					SQLStatement += '"'
				i += 1

			# if len(checkedGroups) > 0: SQLStatement += ', OR '

		i = 1
		if len(checkedGroups) > 0:
			if firstmore == True:
				firstmore = False
			else:
				if len(checkedProjects) > 0:
					SQLStatement += ' OR '
				else:
					SQLStatement += ', '

				firstmore = False
			group = self.ui.cboTreeOp2.currentText()
			fieldname = self.TransLateFieldtoReal(group, True)
			project = self.ui.cboTreeOp1.currentText()
			Projfieldname = self.TransLateFieldtoReal(project, True)

			# SQLStatement = SQLStatement + fieldname + ' = "'
			for item in checkedGroups:
				statement = fieldname + ' = "' + item[1] + '" AND ' + Projfieldname + ' = "' + item[0]
				SQLStatement += statement
				if i < len(checkedGroups):
					SQLStatement += '" OR '
				else:
					SQLStatement += '"'
				i += 1

		i = 1
		if len(checkedSubGroups) > 0:
			if firstmore == True:
				firstmore = False
			else:
				if len(checkedProjects) > 0 or len(checkedGroups) > 0:
					SQLStatement += ' OR '
				else:
					SQLStatement += ', '
				firstmore = False
			Subgroup = self.ui.cboTreeOp3.currentText()
			fieldname = self.TransLateFieldtoReal(Subgroup, True)
			group = self.ui.cboTreeOp2.currentText()
			Groupfieldname = self.TransLateFieldtoReal(group, True)
			project = self.ui.cboTreeOp1.currentText()
			Projfieldname = self.TransLateFieldtoReal(project, True)

			# SQLStatement = SQLStatement + fieldname + ' = "'
			for item in checkedSubGroups:
				statement = fieldname + ' = "' + item[2] + '" AND ' + Groupfieldname + ' = "' + item[
					1] + '" AND ' + Projfieldname + ' = "' + item[0]
				SQLStatement += statement
				if i < len(checkedSubGroups):
					SQLStatement += '" OR '
				else:
					SQLStatement += '"'
				i += 1

		i = 1
		if len(checkedkids) > 0:
			if firstmore == True:
				firstmore = False
			else:
				if len(checkedProjects) > 0 or len(checkedGroups) > 0 or len(checkedSubGroups) > 0:
					SQLStatement += ' OR '
				else:
					SQLStatement += ', '
				firstmore = False

			for item in checkedkids:
				SQLStatement += 'SeqName = "'
				SQLStatement += item
				if i < len(checkedkids):
					SQLStatement += '" OR '
				else:
					SQLStatement += '"'
				i += 1

		return SQLStatement


	def ShowVGenesTextEdit(self, textToShow, style):

		if style == 'aligned':
			FontIs = self.TextEdit.textEdit.currentFont()
			font = QFont(FontIs)

			# FontSize = int(font.pointSize())
			font.setPointSize(10)
			font.setFamily('Courier New')

			self.TextEdit.textEdit.setFont(font)

		elif style == 'standard':
			FontIs = self.TextEdit.textEdit.currentFont()
			font = QFont(FontIs)

			# FontSize = int(font.pointSize())
			font.setPointSize(10)
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
		# self.ui.treeWidget.selectedItems()
		currentitemIs = ''
		for item in value:
			currentitemIs = item.text(0)

		root = self.ui.treeWidget.invisibleRootItem()
		for index in range(root.childCount()):
			fileIs = root.child(index)
			fileIs.setCheckState(0,Qt.Unchecked)



		for index in range(fileIs.childCount()):
			project = fileIs.child(index)  # project level
			project.setCheckState(0,Qt.Unchecked)

			if project.childCount() != 0:  # if project level exists

				for Pindex in range(project.childCount()):  # will iterate through group level
					Group = project.child(Pindex)  # each group
					Group.setCheckState(0,Qt.Unchecked)

					if Group.childCount() != 0:  # if there is a subgroup level
						for Gindex in range(Group.childCount()):  # will iterate through subgroup level
							SubGroup = Group.child(Gindex)  # each subgroup
							SubGroup.setCheckState(0,Qt.Unchecked)

							if SubGroup.childCount() != 0:  # if there is a record level

								for SGindex in range(SubGroup.childCount()):  # will iterate through record level

									Record = SubGroup.child(SGindex)  # each subgroup
									Record.setCheckState(0,Qt.Unchecked)



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

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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


	def RemoveDuplicates(self, DataList):
		from operator import itemgetter
		import itertools

		DataList.sort(key=itemgetter(1, 2, 3))

	@pyqtSlot()
	def SharedClones(self):
		ProjDict  = {}
		i = 1
		fields = ['Project'] #create name dictionary
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		SQLStatement = 'SELECT DISTINCT' + SQLStatement[6:] + ' ORDER BY Project'
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		for item in DataIs:
			NameP = 'P'+ str(i)
			ProjName = item[0]
			DayName =  ProjName[6:8]
			if DayName == 'DS':
				NameP = 'Pre'
			elif DayName == 'D0':
				NameP = 'D0'
			elif DayName == 'D1':
				NameP = 'D1'
			elif DayName == 'D3':
				NameP = 'D3'
			elif DayName == 'D7':
				NameP = 'D7'

			ProjDict[ProjName] = NameP
			i += 1

		fields = ['ClonalPool']
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		SQLStatement = 'SELECT DISTINCT' + SQLStatement[6:] + ' ORDER BY ClonalPool'

		#''SELECT DISTINCT ClonalPool FROM vgenesDB WHERE Project = "Abx09-D0_IgA" ORDER BY ClonalPool''
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		for item in DataIs:
			SharedName = ''
			ClPl = item[0]
			SQLStatement = 'SELECT DISTINCT Project FROM vgenesDB WHERE ClonalPool = ' + ClPl + ' ORDER BY Project'
			DataIs2 = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for proj in DataIs2:
				if proj[0] in ProjDict:
					Pname = ProjDict[proj[0]]
				SharedName += Pname + '_'
			SharedName  = SharedName[:(len(SharedName)-1)]
			if len(SharedName) <= 3:
				SharedName = 'Exclusive'
			elif SharedName == 'D0_D1_D3_D7_Pre':
				SharedName = 'All'


			SQLStatement = 'UPDATE vgenesDB SET Blank7 = "' + SharedName + '" WHERE ClonalPool = "' + ClPl + '"'

			foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

		SQLStatement = 'UPDATE vgenesDB SET Blank7 = "Distinct" WHERE ClonalPool = 0'

		foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)
		print("Finished Shared analysis")

	#SQLCommand = 'UPDATE vgenesDB SET ' + Field + ' = "' + Value + '" WHERE ID = ' + ID





	@pyqtSlot()
	def on_actionFind_Clonal_triggered(self):
		# self.ShowProgressBar()
		from operator import itemgetter
		import itertools
		remove = False
		items = ('Clonal Pools', 'Annotate Duplicates', 'Remove Duplicates', 'Shared Clones')
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
		elif item[:3] == 'Sha':
			self.SharedClones()
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
			          'ID', 'GVend', 'GJbeg', 'GD1beg', 'GD1end']
		else:
			fields = ['SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence',
			          'ID','GVend', 'GJbeg', 'GD1beg', 'GD1end', fieldsearch]

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
		DataIs2 = []

		ProjName = data[75]

		ErLog = 'Clonal analysis for ' + ProjName + '\n'
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
			DataIs2.sort(key=itemgetter(14))
			for k, v in itertools.groupby(DataIs2, key=itemgetter(14)):  # first split out seperate clonal pools
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
						try:
							VGenesSQL.UpdateField(data[119], str(i), 'ClonalPool', DBFilename)
						except:

							print(item + ' caused error in finding clones at line 1798 and so was not annotated as a clone')

					else:
						if j == 1:
							SeqName = 'Duplicate of:  ' + item
							FirstOne = data[119]
						else:
							try:
								if remove == False:
									VGenesSQL.UpdateField(data[119], SeqName, 'Quality', DBFilename)
								else:
									VGenesSQL.UpdateField(data[119], 'Duplicate', 'Quality', DBFilename)
									VGenesSQL.UpdateField(data[119], 'Delete', 'Project', DBFilename)
								DupList += (item + ', ')
							except:
								print('problem line 1810 with: ' + item)
								

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

			self.LoadDB(DBFilename)
			self.ui.txtFieldSearch.setPlainText('Duplicate')
			self.ui.cboFindField.setCurrentText('Quality')
			done = self.on_btnFieldSearch_clicked()
			done = self.on_actionDelete_record_triggered()
			self.on_btnUpdateTree_clicked()

		if len(ErLog2) > 0:
			Erlog2 = ErLog2 + 'The following ' + str(
				Errs) + ' sequences could not be anaylzed for\nclonality because no CDR3s are indicated:\n' + ErLog
			ErlogFile = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'IgBlast', 'database',
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
		items = ('Intraclonal Diversity', 'HS-Summary', "Standard", "Hotspots", 'Mutation Frequencies', 'Heat Map', 'GC-Heat Map',
		         'AnalysisDB Heat Map', 'Reanalyze Base Mutations', 'Make all vairiants')
		title = 'Choose report type:'
		item = setItem(self, items, title)


		if item == "Cancel":
			return
		if item == 'Intraclonal Diversity':
			self.AnalyzeMutations()
		elif item == 'Reanalyze Base Mutations':
			self.UpdateMutationAnalysis()
		elif item == 'Make all vairiants':
			self.MakeAllVairants()
		else:
			VMapHotspots.MapHotspots(self, item, DBFilename, data[0])


	def MakeAllVairants(self):
		SeqName = data[0]
		SeqMutstr = data[97]
		VHGene = data[3]
		Species = data[78]
		beg = int(data[67])
		end = int(data[68]) #todo currently to vend...should be to jend?
		seq  =  data[79]
		seq = seq[:end]

		OutList = []
		SeqNameAp = SeqName + '-1'
		Entry = (SeqNameAp, seq)
		OutList.append(Entry)

		GSeq = data[80][:end]
		GseqName = SeqName + '-0'
		Entry = (GseqName, GSeq)
		OutList.append(Entry)

		SeqMutsR = []
		SeqMuts = VGenesSeq.RScaller(SeqMutstr, VHGene, Species, DBFilename)
		NumMuts = len(SeqMuts)
		for k in range(0, NumMuts - 1):
			Mutation = SeqMuts[k]
			MutDet = Mutation.split('-')
			RS = MutDet[3]
			MutPos = int(MutDet[1])

			if MutPos >= beg and MutPos <= end:
				if RS[0] == 'R':
					SeqMutsR.append(Mutation)
		lenSeq  =len(seq)
		NumMuts = len(SeqMutsR)
		i = 2
		AnoNum = '-' + str(i)

		for k in range(0, NumMuts):
			SeqList = list(seq)
			GSeqList = list(GSeq)
			ModSeqList = SeqList
			Mutation = SeqMutsR[k]
			MutDet = Mutation.split('-')
			print(MutDet)
			SeqNameAp = SeqName + AnoNum

			ModSeqList[int(MutDet[1])-1] = MutDet[0] #remove the mutation
			ModSeq = ''.join(ModSeqList)
			Entry = (SeqNameAp, ModSeq)
			OutList.append(Entry)

			i += 1
			AnoNum = '-' + str(i)
			SeqNameAp = SeqName + AnoNum

			ModSeqList = GSeqList

			ModSeqList[int(MutDet[1])-1] = MutDet[2] #have only that mutation
			ModSeq = ''.join(ModSeqList)
			Entry = (SeqNameAp, ModSeq)
			OutList.append(Entry)

			# ModSeqList = GSeqList
			for j in range(k+1, NumMuts):  #now do with all other mutations
				Mutation1 = SeqMutsR[j]
				MutDet1 = Mutation1.split('-')
				i += 1
				AnoNum = '-' + str(i)
				SeqNameAp = SeqName + AnoNum
				ModSeqList1 = ModSeqList
				ModSeqList1[int(MutDet1[1])-1] = MutDet1[2] #have first above and next only has 1 and 2

				ModSeq = ''.join(ModSeqList1)
				Entry = (SeqNameAp, ModSeq)
				OutList.append(Entry)



				# for l in range(j+1, NumMuts-1): #now do first and second mutate with each consecutively
				# 	Mutation1 = SeqMutsR[l]
				# 	MutDet1 = Mutation1.split('-')
				# 	i += 1
				# 	AnoNum = '-' + str(i)
				# 	SeqNameAp = SeqName + AnoNum
				# 	ModSeqList1 = ModSeqList
				# 	ModSeqList1[int(MutDet1[1])-1] = MutDet1[2]  # have first above and next only has 1 and 2
				#
				# 	ModSeq = ''.join(ModSeqList1)
				# 	Entry = (SeqNameAp, ModSeq)
				# 	OutList.append(Entry)
		print('done')

	@pyqtSlot()
	def AnalyzeMutations(self):

		if self.ui.cboFindField.currentText() == 'Name': self.ui.cboFindField.setCurrentText('Project')
		# SeqName, Sequence, ClonalPool, GermlineSequence, Mutations
		answer = questionMessage(self,
		                         'Use a field to delineate multiple subjects (default = "Project")?\n\n "No" will compare all within each project/group/subgroup.\n\n Press "Cancel" to choose field in the search panel before running analysis.',
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
			          'V1', 'Project', 'Grouping', 'Isotype', 'Quality' , 'Subspecificity', 'Blank8']

		else:
			fields = ['SeqName', 'Sequence', 'ClonalPool', 'GermlineSequence', 'Mutations', 'GVbeg', 'GVend', 'Species',
			          'V1', fieldsearch]

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first

		if answer == 'Yes':
			DataIn.sort(key=itemgetter(9))
		else:
			DataIn.sort(key=itemgetter(9,10,11))
		ClonalPool = []
		ClonalPools = []

		if answer == 'Yes':
			for k, v in itertools.groupby(DataIn, key=itemgetter(9)):  # first split out seperate clonal pools
				# i = int(k)

				if len(k) != 0:
					for item in v:
						ClonalPool.append(item)
					CurrentPool = tuple(ClonalPool)
					ClonalPools.append(CurrentPool)
					ClonalPool.clear()
		else:
			for k, v in itertools.groupby(DataIn, key=itemgetter(9,10,11)):  # first split out seperate clonal pools
				# i = int(k)

				if len(k) != 0:
					for item in v:
						ClonalPool.append(item)
					CurrentPool = tuple(ClonalPool)
					ClonalPools.append(CurrentPool)
					ClonalPool.clear()




			# ClonalPools.append(DataIn)

		f = saveFile(self, 'CSV')
		with open(f, 'w') as currentfile:
			doc = 'Comparison, Project, Subject, Strain, Clonotype, Sequence 1, Sequence 2, Activity, Differences, R-Differences, S-Differences, \
					Begin, End, Length,  Matches, Adjusted Matches,  \
			      Adjusted Differences, R-Matches, Adjusted R-Matches, Adjusted R-Diferences,  \
			      S-Matches, Adjusted S-Matches, Adjusted S-Differences, Warnings, R-Differences, R-Matches, S-Differences, S-Matches, Age\n'

			# CP, Name1, Name2, beg, end, lengthCompared, TotMatch, int(TotAdjMatch), TotDif, int(TotAdjDif), len(
			# 	RMatches),
			# int(RAdjMatch), len(RDifferences), int(RAdjDif), len(SMatches), int(SAdjMatch), len(SDifferences),
			# int(SAdjDif))


		for pool in ClonalPools:
			Pool = list(pool)

			result = VGenesSeq.Intraclonal(Pool, DBFilename)



			if len(result) > 0:

				with open(f, 'a') as currentfile:
					# if answer == 'Yes':
					# 	header = str(Pool[0][9]) + '\n'
					# 	doc += header
					# else:
					# 	header = str(Pool[0][9]) + '_' + str(Pool[0][10]) + '_' + str(Pool[0][11]) + '\n'
					# 	doc += header

					for item in result:
						for i in range(0, 29):
							doc += (str(item[i]) + ', ')
						doc += '\n'
					currentfile.write(doc)
					doc = ''

		self.ShowVGenesText(f)

	@pyqtSlot()
	def on_actionCreateAnalysisDB_triggered(self):
		filename = openFile(self, 'Nucleotide')
		VGenesSQL.CreateAnalysisDB(filename, DBFilename)

	@pyqtSlot()
	def on_actionMultiple_Alignement_triggered(self):
		self.AlignSequences('none')

	@pyqtSlot()
	def on_actionImport_Vgenes_database_triggered(self):
		filename = openFile(self, 'db')
		VGenesSQL.ImportVDB(filename, DBFilename)

	@pyqtSlot()
	def on_actionImportCluster_triggered(self):
		from operator import itemgetter  #from operator import itemgetter   #		SeqList.sort(key=itemgetter(0, 1, 2, 3))


		answerC = informationMessage(self,
		                            'Do you wish to import Seurat Clusters?',
		                            'YN')
		if answerC == 'Yes':
			typeOpen = 'csv'
			filename = openFile(self, typeOpen)

			CSubjects = []
			CClusters = []
			CBarcodes = []

			# LastName = ''
			with open(filename, 'r') as currentfile:
				for row in currentfile:  # imports data from file
					Rawentry = row.strip('\n')

					entryFields = Rawentry.split(',')
					if len(entryFields) == 2:
						if entryFields[0] != 'Barcode' and entryFields[1] != 'Cluster':
							CBarcodes.append(entryFields[0])
							CClusters.append(entryFields[1])
					elif len(entryFields) == 3:
						if entryFields[0] != 'Subject' and entryFields[1] != 'Barcode' and entryFields[2] != 'Cluster':
							CSubjects.append(entryFields[0])
							CBarcodes.append(entryFields[1])
							CClusters.append(entryFields[2])



			myClustersdict = {}

			if len(entryFields) == 2:
				for i in range(len(CBarcodes)):
					myClustersdict[CBarcodes[i]] = CClusters[i]
			if len(entryFields) == 3:

				for i in range(len(CBarcodes)):
					SjBarCode = CSubjects[i] + '-' + CBarcodes[i]
					myClustersdict[SjBarCode] = CClusters[i]




		fields = ['SeqName', 'Id', 'GeneType', 'Blank10']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
		SaidItAlready = 'No'
		DataIn.sort(key=itemgetter(2,0))
		LastName = ''
		DoubletsList = []
		for item in DataIn:
			# doublet = False
			SkipUpdate = False
			SeqName = item[0]
			SeqNameShort = SeqName[:3]
			ID = item[1]
			genetypeis = item[2]
			RealBarCode = item[3]
			TestLen = len(RealBarCode)-2
			TestChar = RealBarCode[TestLen]
			if TestChar =='-':
				RealBarCode = RealBarCode[:len(RealBarCode)-2]

			if len(entryFields) == 3:
				RealBarCode = SeqNameShort + '-' + RealBarCode


			if answerC == 'Yes':
				try:
					RealCluster = myClustersdict[RealBarCode]
				except:
					RealCluster = 'unknown'
			else:
				RealCluster = ''



			FieldName = 'Blank9'
			VGenesSQL.UpdateField(ID, RealCluster, FieldName, DBFilename)





		# import csv
		from operator import itemgetter
		CFilename = ''
		# typeOpen = 'csv'
		# CFilename = openFile(self, typeOpen)

		# QueryIS = 'Enter text to serve as the base name (i.e., subject number or name)'
		# DefaultText = data[75]  # data[77] + '-Expressed'
		# BaseName = setText(self, QueryIS, DefaultText)
		#
		# filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', '10x_barcodes.csv')
		# with open(filename, 'r') as currentfile:
		# 	# myCSVfile = csv.reader(currentfile)
		# 	myCSVfile  = []
		# 	# myCSVfile = currentfile.split(',')
		# 	for row in currentfile:
		# 		entry = row.strip('\n')
		# 		myCSVfile.append(entry)

		# with open(CFilename, 'r') as currentfile:
		#
		# 	myClusterfile  = []
		#
		# 	for row in currentfile:
		# 		entry = row.strip('\n').split(',')
		# 		myClusterfile.append(entry)
		#
		# fields = ['SeqName', 'Id', 'Comments']
		#
		# # checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		# SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		#
		# DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
		#
		# # DataIn.sort(key=itemgetter(2, 0))
		# Seurat = ''
		# CellRanger = ''
		# for item in DataIn:
		# 	SkipUpdate = False
		#
		# 	SeqName = item[0]
		# 	ID = item[1]
		# 	Comments = item[2]
		# 	barCode = Comments[:16]
		#
		#
		# 	try:
		# 		# BarCodeNum = myClusterfile.index(barCode)
		# 		Clusters = [item for item in myClusterfile if item[0] == barCode]
		#
		#
		# 		Seurat = Clusters[0][1]
		# 		CellRanger = Clusters[0][2]
		#
		#
		#
		#
		# 	except:
		# 		SkipUpdate = True
		#
		# 	if SkipUpdate == False:
		# 		FieldName = 'Specificity'
		# 		VGenesSQL.UpdateField(ID, Seurat, FieldName, DBFilename)
		# 		FieldName = 'Subspecificity'
		# 		VGenesSQL.UpdateField(ID, CellRanger, FieldName, DBFilename)


	@pyqtSlot()
	def on_actionRename10x_triggered(self):
		import csv
		from operator import itemgetter  #from operator import itemgetter   #		SeqList.sort(key=itemgetter(0, 1, 2, 3))


		QueryIS = 'Enter text to serve as the base name (i.e., subject number or name)'
		DefaultText = data[75]  # data[77] + '-Expressed'
		BaseName = setText(self, QueryIS, DefaultText)

		filename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', '10x_barcodes.csv')
		with open(filename, 'r') as currentfile:
			# myCSVfile = csv.reader(currentfile)
			myCSVfile  = []
			# myCSVfile = currentfile.split(',')
			for row in currentfile:
				entry = row.strip('\n')
				myCSVfile.append(entry)

		#Now import annotations to get barcodes:

		#filtered_contig_annotations.csv
		answer = informationMessage(self,
		                            'Select annotation file from same folder as consensus.fasta called "outs/filtered_contig_annotations.csv"',
		                            'OK')
		typeOpen = 'csv'
		filename =openFile(self, typeOpen)
		LastBar = ''
		LastName = ''
		with open(filename, 'r') as currentfile:

			barcodeis = []
			NameIs = []
			AllRows = []

			for row in currentfile:  #imports data from file
				Rawentry = row.strip('\n')

				entryFields = Rawentry.split(',')
				if entryFields[0] != 'barcode' and entryFields[17] != 'None':
					AllRows.append(entryFields)

			AllRows.sort(key=itemgetter(17))	#sorts by name so that seques with 2 barcodes can be identified



			for row in AllRows:

				if row[17] == LastName:
					multiBar = LastBar + ',' + row[0]
					LastBar = multiBar
					barcodeis.append(multiBar)

				else:
					barcodeis.append(row[0])
					LastBar = row[0]

				NameIs.append(row[17])


				LastName  = row[17]


			myAnnotationsCSVdict = {}

			for i in range(len(barcodeis)):
				myAnnotationsCSVdict[NameIs[i]] = barcodeis[i]
		# answerC = 'No'
		answerC = informationMessage(self,
		                            'Do you wish to import CellRanger Clusters from 5-prime RNAseq data filename "outs/analysis/clustering/graphclust/clusters.csv"?',
		                            'YN')
		if answerC == 'Yes':
			typeOpen = 'csv'
			filename = openFile(self, typeOpen)

			CClusters = []
			CBarcodes = []

			# LastName = ''
			with open(filename, 'r') as currentfile:
				for row in currentfile:  # imports data from file
					Rawentry = row.strip('\n')

					entryFields = Rawentry.split(',')
					if entryFields[0] != 'Barcode' and entryFields[1] != 'Cluster':
						CBarcodes.append(entryFields[0])
						CClusters.append(entryFields[1])


			myClustersdict = {}

			for i in range(len(CBarcodes)):
				myClustersdict[CBarcodes[i]] = CClusters[i]




		fields = ['SeqName', 'Id', 'GeneType']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
		SaidItAlready = 'No'
		DataIn.sort(key=itemgetter(2,0))
		LastName = ''
		DoubletsList = []
		for item in DataIn:
			doublet = False
			SkipUpdate = False
			SeqName = item[0]
			ID = item[1]
			genetypeis = item[2]
			if SeqName[:9]!= 'clonotype':
				barCode  = SeqName[:16]
				try:
					BarCodeNum = myCSVfile.index(barCode)
				except:
					if SaidItAlready == 'No':
						answer = informationMessage(self,
						                            'This function only renames barecoded sequences from 10x naming',
						                            'OK')
						SaidItAlready = 'yes'
					SkipUpdate = True


			else:
				RealBarCode = myAnnotationsCSVdict[SeqName]
				RealBarCode = RealBarCode[:18]
				if answerC == 'Yes':
					try:
						RealCluster = myClustersdict[RealBarCode]
					except:
						RealCluster = 'unknown'
				else:
					RealCluster = ''

				SeqNameParts = SeqName.split('_')
				barCode = SeqNameParts[0]
				barCode = str(barCode[9:])




			Contig = SeqName[len(SeqName)-1:]





			if SkipUpdate == False:
				if SeqName[:9] != 'clonotype':
					newName = BaseName + '_' + str(BarCodeNum) + '_' + genetypeis[0] + str(Contig)
				else:
					newName = BaseName + '_' + barCode + '_' + genetypeis[0] + str(Contig)

				if newName[:len(newName)-1] == LastName:
					doublet = True
					if genetypeis == 'Heavy':
						DoubletsList.append(barCode)
					FieldName = 'Quality'
					VGenesSQL.UpdateField(LastID, 'doublet', FieldName, DBFilename)
					FieldName = 'SubGroup'
					VGenesSQL.UpdateField(LastID, 'doublet', FieldName, DBFilename)


				if genetypeis != 'Heavy':
					try:
						DoubleHCIndex = DoubletsList.index(barCode)

					except ValueError:
						print("OK")
					else:
						# doublet = True
						FieldName = 'Quality'
						VGenesSQL.UpdateField(ID, 'doublet HC', FieldName, DBFilename)
						FieldName = 'SubGroup'
						VGenesSQL.UpdateField(ID, 'doublet HC', FieldName, DBFilename)


				LastName = newName[:len(newName)-1]





				FieldName = 'Comments'
				VGenesSQL.UpdateField(ID, SeqName, FieldName, DBFilename)

				FieldName = 'SeqName'
				VGenesSQL.UpdateField(ID, newName, FieldName, DBFilename)

				#RealBarCode
				FieldName = 'Blank10'

				VGenesSQL.UpdateField(ID, RealBarCode, FieldName, DBFilename)

				FieldName = 'Blank8'
				VGenesSQL.UpdateField(ID, RealCluster, FieldName, DBFilename)

				if doublet == True:
					FieldName = 'Quality'
					VGenesSQL.UpdateField(ID, 'doublet', FieldName, DBFilename)
					FieldName = 'SubGroup'
					VGenesSQL.UpdateField(ID, 'doublet', FieldName, DBFilename)



				# LastName = newName
				LastID = ID



				model = self.ui.tableView.model()

				model.refresh()

		answer = informationMessage(self, 'Close and restart database to see changes', 'OK')




			# print(newName)



	@pyqtSlot()
	def on_actionMergeMySeq_triggered(self):
		import shutil

		# WorkDir = '/Users/PCW-MacBookProRet/Applications/VGenes/FLASH-1.2.11/'


		filename = openfastq(self)

		read_1 = filename[0]
		read_2 = filename[1]

		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11', 'reads_1.fq')
		shutil.copy(read_1, WorkDir)
		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11', 'reads_2.fq')
		shutil.copy(read_2, WorkDir)

		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11')
		# (dirname, filename) = os.path.split(DBpathname)
		os.chdir(WorkDir)

		CommandLine = "./flash reads_1.fq reads_2.fq 2>&1 | tee flash.log"
		Result = os.popen(CommandLine)
		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11', 'out.extendedFrags.fastq')

		filename = saveFile(self, 'fastq')
		shutil.copy(WorkDir, filename)
		filename2 = filename + 'Unmerged-1.fastq'
		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11',
		                       'out.notCombined_1.fastq')
		shutil.copy(WorkDir, filename2)
		filename2 = filename + 'Unmerged-2.fastq'
		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11',
		                       'out.notCombined_2.fastq')
		shutil.copy(WorkDir, filename2)


		WorkDir = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'FLASH-1.2.11', 'flash.log')
		self.ShowVGenesText(WorkDir)


	#actionAnalyze_Isotypes   actionImport10XInfo
	@pyqtSlot()
	def on_actionImport10XInfo_triggered(self):
		print('10x code')



	@pyqtSlot()
	def on_actionAnalyze_Isotypes_triggered(self):

		fields = ['SeqName', 'Jend', 'Sequence', 'ID']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		for record in  DataIs:
			Sequence = record[2]
			Jend = int(record[1])
			SeqName = record[0]

			IsoSeq = (Sequence[(Jend):])
			# print(SeqName)
			IsoSeq = IsoSeq.strip('N')
			AGCTs = IsoSeq.count('A') + IsoSeq.count('G') + IsoSeq.count('C') + IsoSeq.count('T')
			Isot = ''
			if AGCTs > 5:  # todo decide if can determine isotype from < 5 or need more then
				Isot = VGenesSeq.CallIsotype(IsoSeq)
				print(Isot)

			SQLStatement = 'UPDATE vgenesDB SET Isotype = "' + Isot + '" WHERE SeqName = "' + SeqName + '"'
			# 'UPDATE vgenesDB SET SubGroup = "all" WHERE Project = "Heavy"'
			foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

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
		# import tempfile
		import os
		global GLMsg

		QApplication.setOverrideCursor(Qt.WaitCursor)

		if DataIn == 'none':
			fields = ['SeqName', 'Sequence']
			# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
			SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
			DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first

			if len(DataIs) == 1:
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

		# import subprocess


		# (fd, outfilename) = tempfile.mkstemp()
		try:


			outfilename = VGenesSeq.ClustalO(DataIs, 80, True)




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
		except:
			print('no')

		finally:
			os.remove(outfilename)

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
		SelectedItems =  self.getTreeSelected()
		NumSelected  = len(SelectedItems)
		if NumSelected > 1:
			NewHead  = str(NumSelected) + ' items selected'
			self.ui.label_Name.setText(NewHead)

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
			filename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'RecentPaths.vtx')
			with open(filename, 'r') as currentfile:
				vv = currentfile

		except:
			filename = 'RecentPaths.vtx'

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
			if DBFilename == None:
				return

		if DBFilename != None:
			if os.path.isfile(DBFilename):
				self.LoadDB(DBFilename)
			else:
				VGenesSQL.creatnewDB(DBFilename)

			self.UpdateRecentList(DBFilename, True)
			if GetName == True:
				self.SaveBackup
		else:
			self.hide()
			#self.ApplicationStarted()


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
			if currentRow < records:
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
		AASeq, ErMessage = VGenesSeq.Translator(data[80], 0)
		self.ui.txtAASeq.setText(AASeq)
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

	def on_textBarcode_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textBarcode.toPlainText()
		LastSelected = ('Blank10', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)


	# if data[106] != 'Blank8':
	# 	self.ui.textCluster.setText(data[106])
	# if data[107] != 'Blank9':
	# 	Population

	def on_textCluster_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textCluster.toPlainText()
		LastSelected = ('Blank8', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtLabel_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textEdit.toPlainText()
		LastSelected = ('Blank12', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtStatus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textEdit.toPlainText()
		LastSelected = ('Blank13', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_textEdit_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textEdit.toPlainText()
		LastSelected = ('Blank9', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtPopulation_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtPopulation.toPlainText()
		LastSelected = ('Blank11', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_textMutations_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textMutations.toPlainText()
		LastSelected = ('TotMut', valueis)
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
	def on_txtComments_textChanged(self):

		valueTo = self.ui.txtComments.toPlainText()
		FieldIS = 'Comments'
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
	def on_textBarcode_textChanged(self):

		valueTo = self.ui.textBarcode.toPlainText()
		FieldIS = 'Blank10'
		self.FieldChanger(valueTo, FieldIS)

	# if data[106] != 'Blank8':
	# 	self.ui.textCluster.setText(data[106])
	# if data[107] != 'Blank9':
	# 	self.ui.textEdit.setText(data[107])

	@pyqtSlot()
	def on_textCluster_textChanged(self):

		valueTo = self.ui.textCluster.toPlainText()
		FieldIS = 'Blank8'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtPopulation_textChanged(self):

		valueTo = self.ui.txtPopulation.toPlainText()
		FieldIS = 'Blank11'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_textMutations_textChanged(self):

		valueTo = self.ui.textMutations.toPlainText()
		FieldIS = 'TotMut'
		self.FieldChanger(valueTo, FieldIS)


	@pyqtSlot()
	def on_textEdit_textChanged(self):

		valueTo = self.ui.textEdit.toPlainText()
		FieldIS = 'Blank9'
		self.FieldChanger(valueTo, FieldIS)

	# @pyqtSlot()
	# def on_txtStatus_textChanged(self):
	#
	# 	valueTo = self.ui.txtStatus.toPlainText()
	# 	FieldIS = 'Blank13'
	# 	self.FieldChanger(valueTo, FieldIS)
	#
	# @pyqtSlot()
	# def on_txtLabel_textChanged(self):
	#
	# 	valueTo = self.ui.txtLabel.toPlainText()
	# 	FieldIS = 'Blank12'
	# 	self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtQuality_textChanged(self):

		valueTo = self.ui.txtQuality.toPlainText()
		FieldIS = 'Quality'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtLabel_textChanged(self):

		valueTo = self.ui.txtLabel.toPlainText()
		FieldIS = 'Blank12'
		self.FieldChanger(valueTo, FieldIS)

	@pyqtSlot()
	def on_txtStatus_textChanged(self):

		valueTo = self.ui.txtStatus.toPlainText()
		FieldIS = 'Blank13'
		self.FieldChanger(valueTo, FieldIS)

	# @pyqtSlot()
	def on_sbJend_valueChanged(self):

		valueTo = str(self.ui.sbJend.value())
		FieldIS = 'Jend'
		self.FieldChanger(valueTo, FieldIS)
		Jend = int(valueTo)
		VSeq = data[79]
		JendSeq = VSeq[Jend-11:Jend]
		JendAASeq, ErMessage = VGenesSeq.Translator(JendSeq, 0)
		JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + ' \n' + JendSeq[0:3] + ' ' + JendSeq[3:6] + ' ' + JendSeq[6:9]

		self.ui.txtJend_2.setText(JendDisplay)

	@pyqtSlot()
	def on_txtQuality_2_textChanged(self):

		valueTo = self.ui.txtQuality_2.toPlainText()
		FieldIS = 'Quality'
		self.FieldChanger(valueTo, FieldIS)

	def FieldChanger(self, valueTo, FieldIS):
		global FieldChanged, UpdateSpecific
		if JustMovedIt == False:
			FieldChanged = True
			self.UpdateFChanges(data[119], valueTo, FieldIS)
			if FieldIS == 'SeqName':
				self.NameChange(valueTo)
			if FieldIS == 'Specificity' or FieldIS == 'Subspecificity':
				UpdateSpecific = True

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


	@pyqtSlot()
	def on_actionFixNames_triggered(self):
		fields = ['SeqName', 'Id', 'GeneType']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first

		DataIn.sort(key=itemgetter(2,0))




		for item in DataIn:

			SkipUpdate = False
			SeqName = item[0]
			ID = item[1]
			IsCL = SeqName[4:6]
			if IsCL == 'Cl':
				NewName = SeqName[:4] + SeqName[6:]
				VGenesSQL.UpdateField(ID, NewName, 'SeqName', DBFilename)

		self.model.refresh()



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

	def on_txtIsotype_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtIsotype.toPlainText()
		LastSelected = ('Isotype', valueis)
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

	def on_txtStatus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtStatus.toPlainText()
		LastSelected = ('Blank13', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def on_txtLabel_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtLabel.toPlainText()
		LastSelected = ('Blank12', valueis)
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

	@pyqtSlot()
	def on_cboFindField_currentTextChanged(self):
		if self.ui.cboFindField.currentText() == 'Specificity':
			self.SpecSet()
		if self.ui.cboFindField.currentText() == 'Subspecificity':
			self.SubspecSet()

		if self.ui.cboFindField.currentText() == 'Autoreactivity':
			self.AutoRXSet()





	def SpecSet(self):
		global LastSelected


		# if JustMoved == False:
		valueis = self.ui.listViewSpecificity.currentText()
		LastSelected = ('Specificity', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def SubspecSet(self):
		global LastSelected


		# if JustMoved == False:
		valueis = self.ui.listViewSpecificity_2.currentText()
		LastSelected = ('Subspecificity', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def AutoRXSet(self):
		global LastSelected


		# if JustMoved == False:
		valueis = self.ui.Autoreactivity.currentText()
		LastSelected = ('Blank6', valueis)
		field = LastSelected[0]
		Fiedlvalue = self.TransLateFieldtoReal(field, False)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	# @pyqtSlot()
	def on_listViewSpecificity_currentTextChanged(self):
		valueTo = self.ui.listViewSpecificity.currentText()
		FieldIS = 'Specificity'
		self.FieldChanger(valueTo, FieldIS)
	#
	# @pyqtSlot()
	# def on_lblSpecificity_clicked(self):
		self.ui.cboFindField.setCurrentText('Specificity')

	# @pyqtSlot()
	def on_listViewSpecificity_2_currentTextChanged(self):
		valueTo = self.ui.listViewSpecificity_2.currentText()
		FieldIS = 'Subspecificity'
		self.FieldChanger(valueTo, FieldIS)


	def on_Autoreactivity_currentTextChanged(self):
		valueTo = self.ui.Autoreactivity.currentText()
		FieldIS = 'Blank6'
		self.FieldChanger(valueTo, FieldIS)

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
		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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

						outfilename = VGenesSeq.ClustalO(ToClustal, 1000, True)
						# ClustalOut = VGenesSeq.ClustalO(ToClustal, 1000, False)
						ToClustal.clear()
						# outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega',
						#                            'my-out-seqs.fa')
						Aligned = VGenesSeq.readClustalOutput(outfilename)

						os.remove(outfilename)

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
			# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
			SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
			foundRecs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for item in foundRecs:
				SeqName = item[0]
				Subgroup = item[1]
				NewSub = Subgroup + EditSubgroup

				# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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
		global wasClicked
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
		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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
		# self.ui.cboReportOptions.currentIndexChanged()  'FASTA Nucleotide Rename file'
		option = self.ui.cboReportOptions.currentText()
		if option == 'FASTA Nucleotide file':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)

		if option == '10x Synthesis report':
			self.ui.ckReportCSV.setDisabled(True)
			self.ui.ckReportCSV.setChecked(False)
			self.ui.ckReportDisplay.setDisabled(False)
			self.ui.ckReportDisplay.setChecked(True)
			self.ui.ckReportText.setDisabled(False)
			self.ui.ckReportText.setChecked(True)

		if option == 'FASTA Nucleotide Rename file':
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
		Backfilename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects',  'VGenes', 'BackUP.vdb')
		global DBFilename

		if DBFilename != None:
			shutil.copy(DBFilename, Backfilename)

	@pyqtSlot()
	def on_actionRevert_to_previous_triggered(self):
		import shutil

		buttons = 'OKC'
		answer = informationMessage(self, 'Any work done since opening this instance of VGenes will be lost.', buttons)

		if answer == 'Cancel':
			return

		Backfilename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'BackUP.vdb')
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
		global DBFilename
		IgBLASTAnalysis = []
		fields = ['*']
		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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
		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
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
				if data[108] != 'Blank10':
					self.ui.textBarcode.setText(data[108])
				if data[106] != 'Blank8':
					self.ui.textCluster.setText(data[106])

				if data[107] != 'Blank9':
					self.ui.textEdit.setText(data[107])

				if data[109] != 'Blank11':
					self.ui.txtPopulation.setText(data[109])
				else:
					self.ui.txtPopulation.setText("")

				VSeq = data[79]
				GVSeq = data[80]
				# AASeq, ErMessage = VGenesSeq.Translator(VGenesSeq, 0)

				try:
					Vbeg = int(data[67])

					GVbeg = int(data[59])
					self.ui.sbVbeg.setValue(Vbeg)
					Vend = int(data[68])
					self.ui.sbVend.setValue(Vend)

					VBegSeq = VSeq[0:9]
					VBegAASeq, ErMessage = VGenesSeq.Translator(VBegSeq, 0)
				except:
					print('none error')
				try:
					VbegDisplay = ' ' + VBegAASeq[0] + '   ' + VBegAASeq[1] + '   ' + VBegAASeq[2] + ' \n' + VBegSeq[0:3] + ' ' + VBegSeq[3:6] + ' ' + VBegSeq[6:9]

				except:
					print('oops')

				GVBegSeq = GVSeq[0:9]
				GVBegAASeq, ErMessage = VGenesSeq.Translator(GVBegSeq, 0)
				try:
					GVbegDisplay = ' ' + GVBegAASeq[0] + '   ' + GVBegAASeq[1] + '   ' + GVBegAASeq[2] + ' \n' + GVBegSeq[0:3] + ' ' + GVBegSeq[3:6] + ' ' + GVBegSeq[6:9]

				except:
					print('oops')


				try:
					self.ui.txtVbeg.setText(VbegDisplay)
					self.ui.txtVExp.setText(GVbegDisplay)


					Dbeg = int(data[69])
					self.ui.sbDbeg.setValue(Dbeg)
					Dend = int(data[70])
					self.ui.sbDend.setValue(Dend)

					Jbeg = int(data[73])

					self.ui.sbJbeg.setValue(Jbeg)
					Jend = int(data[74])
					GJend = int(data[66])
					self.ui.sbJend.setValue(Jend)

					# VSeq = data[79]
					JendSeq = VSeq[Jend-11:Jend]
					JendAASeq, ErMessage = VGenesSeq.Translator(JendSeq, 0)
					JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + ' \n' + JendSeq[0:3] + ' ' + JendSeq[3:6] + ' ' + JendSeq[6:9]

					self.ui.txtJend_2.setText(JendDisplay)
				except:
					print('none')
				# GJendSeq = GVSeq[GJend-11:GJend]
				# GJendAASeq, ErMessage = VGenesSeq.Translator(GJendSeq, 0)
				# GJendDisplay = ' ' + GJendAASeq[0] + '   ' + GJendAASeq[1] + '   ' + GJendAASeq[2] + ' \n' + GJendSeq[0:3] + ' ' + GJendSeq[3:6] + ' ' + GJendSeq[6:9]
				#
				# self.ui.txtJExp.setText(GJendDisplay)
				#


				self.ui.txtProject.setText(data[75])
				self.ui.textMutations.setText(data[57])


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
				try:
					if UpdateSpecific == True:
						self.PopulateSpec()
				except:
					print('NoTree')
				self.ui.listViewSpecificity.setCurrentText(data[86])
				self.ui.listViewSpecificity_2.setCurrentText(data[87])
				self.ui.Autoreactivity.setCurrentText(data[104])

				if self.ui.btnEditSeq.isChecked():
					msg = 'Sequence edit mode was activated, do you want to save changes and re-analyze this sequence before proceeding?'
					buttons = 'YN'
					answer = questionMessage(self, msg, buttons)
					if answer == 'Yes':
						self.ui.btnEditSeq.setChecked(False)
						self.UpdateSeqAnalysis()

						self.ui.btnEditSeq.setText("Edit Mode")
						msg = 'Press "Edit Mode" to edit sequence.'

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

				if data[110] != 'Blank12':
					self.ui.txtLabel.setText(data[110])
				else:
					self.ui.txtLabel.setText('')

				if data[111] != 'Blank13':
					self.ui.txtStatus.setText(data[111])
				else:
					self.ui.txtStatus.setText('')



				self.ui.txtComments.setText(data[94])
				self.ui.txtID.setText(data[98])
				self.ui.txtCDR3MW.setText(data[99])
				self.ui.txtCDR3pI.setText(data[100])
				self.ui.txtIsotype.setText(data[101])
				self.ui.txtIsotypeSeq.setText(data[101])

				valueToR = int(self.ui.sbPairRight.value())
				valueToL = int(self.ui.sbPairLeft.value())
				BarCode = data[0]
				BarCode = BarCode[:len(BarCode)- valueToR]
				BarCode = BarCode[valueToL:]
				self.ui.txtPairNames.setText(BarCode)
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
				self.on_cboFindField_currentTextChanged()



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

		Pathname = openFile(self, 'Nucleotide')
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
	def on_sbPairLeft_valueChanged(self):

		valueToL = int(self.ui.sbPairLeft.value())
		valueToR = int(self.ui.sbPairRight.value())
		BarCode = data[0]
		BarCode = BarCode[valueToL:]
		BarCode = BarCode[:len(BarCode)- valueToR]
		self.ui.txtPairNames.setText(BarCode)


	# @pyqtSlot()
	def on_sbPairRight_valueChanged(self):

		valueToR = int(self.ui.sbPairRight.value())
		valueToL = int(self.ui.sbPairLeft.value())
		BarCode = data[0]
		BarCode = BarCode[:len(BarCode)- valueToR]
		BarCode = BarCode[valueToL:]
		self.ui.txtPairNames.setText(BarCode)


	@pyqtSlot()
	def on_btnPairRename_clicked(self):
		msg = 'This function will pair Heavy/Light chain sequences, putting a common name in the barcode field so that antibody genes can be synthesized from both. For pairs with common portions of the names (i.e. Patient2C04H and Patient2C04K) use the the left and right spinboxes to truncate names to a common but distinct name (i.e., Patient2C04) for the heavy/light chains, ensuring only that single heavy/light pair will share the truncated name portion. All checked records will be searched to find matched pairs and the barcodes annotated. Continue (if a proper common name is indicated)?'
		buttons = 'YN'
		answer = informationMessage(self, msg, buttons)

		if answer == 'Yes':
			CurrName = self.ui.txtPairNames.toPlainText()
			valueToR = int(self.ui.sbPairRight.value())
			valueToL = int(self.ui.sbPairLeft.value())

			fields = ['SeqName', 'GeneType', 'Blank10']
			SequenceName = data[0]
			SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
			SQLStatement += ' ORDER BY Blank10'
			DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			NamedSets = []

			for Seq in DataIs:
				answer = 'Yes'
				Barcode = Seq[2]
				Seqname1 = Seq[0]


				NewBarCode = Seqname1
				NewBarCode = NewBarCode[:len(NewBarCode) - valueToR]
				NewBarCode = NewBarCode[valueToL:]
				search = NewBarCode
				for i in range(0, valueToR):
					NewBarCode = NewBarCode + '%'
				for i in range(0, valueToL):
					NewBarCode = '%' + NewBarCode

				SQLStatement = 'SELECT SeqName FROM vgenesDB WHERE SeqName LIKE "'+ NewBarCode + '"'
				HeavyLights = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				fieldsearch = "Blank10"

				for record in HeavyLights:


					SeqIsNamed = record[0]

					# self.ui.txtFieldSearch.setPlainText(SeqIsNamed)
					# self.ui.cboFindField.setCurrentText('Name')
					# done = self.on_btnFieldSearch_clicked()




					# fields = ['ID']
					# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
					# SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SeqIsNamed)
					SQLStatement = 'SELECT ID FROM vgenesDB WHERE SeqName = "'+ SeqIsNamed +'"'

					WhereStart = SQLStatement.find('WHERE')
					WhereState = SQLStatement[WhereStart - 1:]  # + ' AND '
					SQLStatement = 'UPDATE vgenesDB SET ' + fieldsearch + ' = "' + search + '"' + WhereState  # -7.000000' WHERE locId = 173567"
					# ' WHERE SeqName = "A116_1B04H-2" OR SeqName = "A116_1B04H-3"'
					foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

					model = self.ui.tableView.model()

					model.refresh()





	# self.ui.txtPairNames.setText(data[0])

	@pyqtSlot()
	def on_btn10xEditPreview_clicked(self):
		global PreviewHTcurrent
		global PreviewCurrentType
		global PreviewHTExp
		LSeqName = ''
		LGeneType = ''
		LVbeg = ''
		LVExp = ''
		LJend = ''
		LJexp = ''
		LCbeg = ''
		LCexp = ''
		LSequence = ''


		self.ui.tabWidget.setCurrentIndex(5)

		fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'Jend', 'Blank10']
		SequenceName = data[0]
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SequenceName)
		SQLStatement += ' ORDER BY Blank10'
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		JendMuts = []
		ExpectIgG = ' A   S   T   Q \n' + 'GCG TCG ACC AAG'
		ExpectKappa = ' R   T   V   A \n' + 'CGT ACG GTG GCA'
		ExpectLambda = ' G   Q   P   K \n' + 'GGT CAG CCC AAG'

		ExpIgGAA = 'ASTQ'
		ExpCKappaAA = 'RTVA'
		ExpCLambdaAA = 'GQPK'

		# GJseq = []
		# grab last 9 of J gene as going into report and last nine of germline J gene actual
		GermJH1 = '' + ''
		GermJH2 = '' + ''
		GermJH3 = ''+''
		GermJH4 = ''+''
		GermJH5 = ''+''
		GermJH6 = ''+''


		NameList = []
		# for item in fields:
		#     NameList.append(str(self.TransLateFieldtoReal(item, False)))
		CSVOut = 'Sequence name, Sequence \n'

		# preSeq = 'gcaactggtgtacattcc'  # 'CTG CAA CCG GTG TAC ATT CA'

		if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
			preSeq = 'gcaactggtgtacattcc'  # '
		elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':
			preSeq = 'atcctttttctagtagcaactgcaaccggtgtacattcc'  # '
		elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':
			preSeq = 'gcaaccggtgtacattcc'  # '

		LeaderSeq = ' G   V   H   S \nGGT GTA CAT TCC'



		PreviewHTExp.clear()
		PreviewEntry = ()
		#SeqName 0, HGeneType 1, HLeader 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,
		# LGeneType 10, LVbeg 11, LVexp 12, LJend 13, LJexp 14, LCbeg 15, LCexp 16, LSequence 17


		# todo change to app folder
		ErlogFile = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'IgBlast', 'database',
		                         'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
		ErLog = 'VGenes input beginning at: ' + time.strftime('%c') + '\n'
		with open(ErlogFile, 'w') as currentFile:  # using with for this automatically closes the file even if you crash
			currentFile.write(ErLog)

		try:
			DBpathname = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'VDJGenes.db')

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

		#  then need to create a cursor that lets you traverse the database
		cursor = conn.cursor()


		for Seq in DataIs:
			answer = 'Yes'
			Barcode = Seq[8]
			Seqname1 = Seq[0]
			# fields = ['SeqName', 'GeneType', 'V1', 'J1', 'productive', 'Sequence', 'Vbeg', 'GermlineSequence', 'Blank10']
			SQLStatement = 'SELECT SeqName, GeneType, V1, J1, productive, Sequence, Vbeg, GermlineSequence, Blank10, IDEvent, Mutations, Jend FROM vgenesDB WHERE Blank10 = "' + Barcode + '" ORDER BY GeneType'
			SeqClone = VGenesSQL.RunSQL(DBFilename, SQLStatement)


			if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':

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
				# SeqName 0, HGeneType 1, HLeader 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,
				# LGeneType 10, LVbeg 11, LVexp 12, LJend 13, LJexp 14, LCbeg 15, LCexp 16, LSequence 17
				for Record in SeqClone:
					done = False
					answer2 = 'Yes'
					answer3 = 'Yes'
					SeqName = Record[0]

					GeneType = Record[1]
					Productive = Record[4]
					Sequence = Record[5]
					Sequence = Sequence.upper()
					Sequence = Sequence.replace('-','')
					Vbeg = int(Record[6])
					GermSeq = Record[7]
					JendNumber = int(Record[11])

					VgeneName = Record[2]
					JgeneName = Record[3]



					SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + VgeneName + '"'
					cursor.execute(SqlStatementV)

					for row in cursor:
						GVgene = row[4]
						VStart = int(row[3])

					VExp = GVgene[VStart-1:12]


					VExpAASeq, ErMessage = VGenesSeq.Translator(VExp, 0)
					VExpDisplay = ' ' + VExpAASeq[0] + '   ' + VExpAASeq[1]  + '   ' + VExpAASeq[2] + '   ' + VExpAASeq[3] + ' \n' + VExp[0:3] + ' ' + VExp[3:6] + ' ' + VExp[6:9] + ' ' + VExp[9:12]

					# SeqName 0, HGeneType 1, LeaderSeq 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,
					# LGeneType 10, LVbeg 11, LVexp 12, LJend 13, LJexp 14, LCbeg 15, LCexp 16, LSequence 17



					Vend = Sequence[:12]

					VendAASeq, ErMessage = VGenesSeq.Translator(Vend, 0)
					VendDisplay = ' ' + VendAASeq[0] + '   ' + VendAASeq[1] + '   ' + VendAASeq[2] + '   ' + VendAASeq[
						3] + ' \n' + Vend[0:3] + ' ' + Vend[3:6] + ' ' + Vend[6:9] + ' ' + Vend[9:12]



					SqlStatementJ = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + JgeneName + '"'
					cursor.execute(SqlStatementJ)
					# GJseq.clear()

					for row in cursor:
						GJgene = row[4]

					if GeneType == 'Heavy' or GeneType == 'Lambda':
						if GJgene[len(GJgene) - 1] == 'G':
							GJgene = GJgene[:len(GJgene) - 1]

					if GeneType == 'Kappa':
						if GJgene[len(GJgene) - 1] == 'C':
							GJgene = GJgene[:len(GJgene) - 1]

					JExp = GJgene[len(GJgene)-12:len(GJgene)]
					JExpAASeq, ErMessage = VGenesSeq.Translator(JExp, 0)
					JExpDisplay = ' ' + JExpAASeq[0] + '   ' + JExpAASeq[1]  + '   ' + JExpAASeq[2] + '   ' + JExpAASeq[3] + ' \n' + JExp[0:3] + ' ' + JExp[3:6] + ' ' + JExp[6:9] + ' ' + JExp[9:12]





					# todo Figure GL for lambda and kappas and also fix Jend in database...
					# if GermSeq[len(GermSeq) - 1] == 'G':
					# 	Jend = len(GermSeq) - 1
					# else:
					# 	Jend = len(GermSeq)



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
					JendMut = False
					if answer2 == 'Yes' and answer3 == 'Yes':
						# if IDevent == 'Insertion':
							# msg = SeqName + 'Has an insertion and so needs to be prepared for synthesis by hand and will not be included.'
							# buttons = 'OK'
							# Sequence = ''
						Mutset = Mutations.split(',')
						for mutation in Mutset:
							MutDetails = mutation.split('-')
							for MutD in MutDetails:
								if MutD == 'Insertion':
									Insert = MutDetails[2]
									JendNumber += len(Insert)

							try:
								if int(MutDetails[1]) > JendNumber -15:
									JendMut = True
							except:
								print('hmmm')

						Sequence = Sequence.replace('-', '')
						Sequence = Sequence[:JendNumber]

						# SeqName 0, HGeneType 1, LeaderSeq 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,
						# LGeneType 10, LVbeg 11, LVexp 12, LJend 13, LJexp 14, LCbeg 15, LCexp 16, LSequence 17

						if GeneType == 'Heavy':

							if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
								postSeq = 'gcgtcgaccaagggcc'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':
								postSeq = 'gcgtcgaccaagggcccatcggtcttcc'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':
								postSeq = 'gcgtcgaccaagg'  # '




							# JendNumber -= 2
							Jend = Sequence[JendNumber - 12:]
							JendNumber2 = JendNumber
							if JendMut == False:
								if Jend != JExp:
									for i in range(0, 5):
										# Sequence = Sequence[:len(Sequence)-1]
										JendNumber2 = JendNumber2 - 1
										Jend = Sequence[JendNumber2 - 12:JendNumber2]
										if Jend == JExp:
											Sequence = Sequence[:JendNumber2]
											Jend = Sequence[JendNumber2 - 12:]
											break


							JendAASeq, ErMessage = VGenesSeq.Translator(Jend, 0)
							JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + '   ' + JendAASeq[3] + ' \n' + Jend[0:3] + ' ' + Jend[3:6] + ' ' + Jend[6:9] + ' ' + Jend[9:12]


							if len(Sequence) % 3 == 0:
								HCbeg = ExpectIgG
							else:
								HCbeg = 'Out-of-frame'



							HSeqName = SeqName
							HGeneType = 'Heavy'
							HCExp = ExpectIgG
							HVbeg = VendDisplay
							HVExp = VExpDisplay
							HJend = JendDisplay
							HJexp = JExpDisplay

							if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
								HSequence = Sequence + postSeq
							elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':
								HSequence = preSeq + Sequence + postSeq
							elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':
								HSequence = preSeq + Sequence + postSeq  # '






						elif GeneType == 'Kappa':
							# postSeq = 'cgtacggtggcacagaaccggtgtccattcc'  # CGT ACG gtg gca cag aAC CGG TGTcCATTCC 'gtacggtggctgcaccatctgtctt' # gtacggtggc'   #AA GAC AGA TGG TGC AGC CAC CGT ACG    CGTACGGTGGCTGCACCATCTGTCTT

							if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
								postSeq = 'cgtacggtggcacagaaccggtgtccattcc'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':
								postSeq = 'cgtacggtggctgcaccatctgtctt'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':
								postSeq = 'cgtacggtggctgc'  # '



							done = True



							Jend = Sequence[JendNumber - 12:]
							JendNumber2 = JendNumber
							if JendMut == False:
								if Jend != JExp:
									for i in range(0, 5):
										# Sequence = Sequence[:len(Sequence)-1]
										JendNumber2 = JendNumber2 - 1
										Jend = Sequence[JendNumber2 - 12:JendNumber2]
										if Jend == JExp:
											Sequence = Sequence[:JendNumber2]
											Jend = Sequence[JendNumber2 - 12:]
											break



							JendAASeq, ErMessage = VGenesSeq.Translator(Jend, 0)
							JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + '   ' + JendAASeq[3] + ' \n' + Jend[0:3] + ' ' + Jend[3:6] + ' ' + Jend[6:9] + ' ' + Jend[9:12]

							InFrameSeq = len(Sequence)/3

							if len(Sequence) % 3 == 0:
								LCbeg = ExpectKappa
							else:
								LCbeg = 'Out-of-frame'

							# SeqName 0, HGeneType 1, LeaderSeq 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,
							# LGeneType 10, LVbeg 11, LVexp 12, LJend 13, LJexp 14, LCbeg 15, LCexp 16, LSequence 17

							LSeqName = SeqName
							LGeneType = 'Kappa'
							LCexp = ExpectKappa
							LVbeg = VendDisplay
							LVExp = VExpDisplay
							LJend = JendDisplay
							LJexp = JExpDisplay



							LSequence = preSeq + Sequence + postSeq


						elif GeneType == 'Lambda':
							# postSeq = 'ggtcagcccaaggccaaccccactgtcactctgttcccgccctcgaggtggcacagaaccggtgtccattcc'  # 'ggtcagcccaaggctgccccctcggtcactctgttcccrccctcgagtgaggagcttcaagccaaca' #'ggtcagcccaaggctgccccctcggtcactctgttcccaccctcgagtgaggag'   #TG TTG GCT TGA AGC TCC TCA CTC GAG GGY GGG AAC AGA GTG  cactctgttcccrccctcgagtgaggagcttcaagccaaca

							if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
								postSeq = 'ggtcagcccaaggccaaccccactgtcactctgttcccgccctcgaggtggcacagaaccggtgtccattcc'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':
								postSeq = 'ggtcagcccaaggctgccccctcggtcactctgttcccrccctcgagtgaggagcttcaagccaaca'  # '
							elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':
								postSeq = 'ggtcagcccaaggctgccccctcggtcactctgttcccrccctcgagtgaggagc'  # '



							Jend = Sequence[JendNumber - 12:]
							JendNumber2 = JendNumber
							if JendMut == False:
								if Jend != JExp:
									for i in range(0, 5):
										# Sequence = Sequence[:len(Sequence)-1]
										JendNumber2 = JendNumber2 - 1
										Jend = Sequence[JendNumber2 - 12:JendNumber2]
										if Jend == JExp:
											Sequence = Sequence[:JendNumber2]
											Jend = Sequence[JendNumber2 - 12:]
											break



							JendAASeq, ErMessage = VGenesSeq.Translator(Jend, 0)
							JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + '   ' + JendAASeq[3] + ' \n' + Jend[0:3] + ' ' + Jend[3:6] + ' ' + Jend[6:9] + ' ' + Jend[9:12]


							if len(Sequence) % 3 == 0:
								LCbeg = ExpectLambda
							else:
								LCbeg = 'Out-of-frame'




							LSeqName = SeqName
							LGeneType = 'Lambda'
							LCexp = ExpectLambda
							LVbeg = VendDisplay
							LVExp = VExpDisplay
							LJend = JendDisplay
							LJexp = JExpDisplay

							LSequence = preSeq + Sequence + postSeq





			if answer == 'Yes':
				PreviewEntry = (HSeqName, HGeneType, LeaderSeq, HVbeg, HVExp, HJend, HJexp, HCbeg, HCExp, HSequence, LSeqName, LGeneType, LVbeg, LVExp, LJend, LJexp, LCbeg, LCexp, LSequence)
				PreviewHTExp.append(PreviewEntry)

					# SeqName 0, HGeneType 1, HLeader 2, HVbeg 3, HVexp 4, HJend 5, HJexp 6, HCbeg 7, HCexp 8, HSequence 9,

					# LSeqName 10, LGeneType 11, LVbeg 12, LVexp 13, LJend 14, LJexp 15, LCbeg 16, LCexp 17, LSequence 18

		#
		PreviewHTcurrent = 0
		PreviewCurrentType = 'H'
		self.PopulatePreviewHT('forward')


	@pyqtSlot()
	def on_btn10xEditForward_clicked(self):
		self.PopulatePreviewHT('forward')

	@pyqtSlot()
	def on_btn10xEditBack_clicked(self):
		self.PopulatePreviewHT('backward')

	def PopulatePreviewHT (self, direction):
		global PreviewHTcurrent
		global PreviewCurrentType
		global PreviewHTExp

		self.ui.lblVbegOff.setEnabled(False)
		self.ui.lblJendOff.setEnabled(False)
		# PreviewHTcurrent = 0
		# PreviewCurrentType == 'H'
		# PreviewEntry = (
		# HSeqName 0, HGeneType 1, LeaderSeq 2, HVbeg 3, HVExp 4, HJend 5, HJexp 6, HCbeg 7,
		# HCExp 8, HSequence 9, LSeqName 10, LGeneType 11, LVbeg 12,
		# LVExp 13, LJend 14, LJexp 15, LCbeg 16, LCexp 17, LSequence 18)

		NumSeqs = len(PreviewHTExp)

		CurrentRecord = PreviewHTExp[PreviewHTcurrent]


		HSeqName = CurrentRecord[0]
		LSeqName = CurrentRecord[10]
		HGeneType = CurrentRecord[1]
		LGeneType = CurrentRecord[11]

		LeaderSeq = CurrentRecord[2]
		HVBeg = CurrentRecord[3]
		HVExp = CurrentRecord[4]
		HJEnd = CurrentRecord[5]
		HJExp = CurrentRecord[6]
		HCBeg = CurrentRecord[7]
		HCExp = CurrentRecord[8]

		LVBeg = CurrentRecord[12]
		LVExp = CurrentRecord[13]
		LJEnd = CurrentRecord[14]
		LJExp = CurrentRecord[15]
		LCBeg = CurrentRecord[16]
		LCExp = CurrentRecord[17]

		self.ui.txtEditLeader_1.setText(LeaderSeq)
		self.ui.txtLeader_Expected.setText(LeaderSeq)


		if PreviewCurrentType == 'H':

			self.ui.txtFieldSearch.setPlainText(HSeqName)
			self.ui.cboFindField.setCurrentText('Name')
			done = self.on_btnFieldSearch_clicked()

			self.ui.lblSeqName.setText(HSeqName)
			self.ui.lblPrevType.setText(HGeneType)

			self.ui.txtEditVbeg_1.setText(HVBeg)
			self.ui.txtVbeg_Expected.setText(HVExp)

			if HVBeg != HVExp:
				self.ui.lblVbegOff.setEnabled(True)

			self.ui.txtEditJend_1.setText(HJEnd)
			self.ui.txtJendExp.setText(HJExp)

			if HJEnd != HJExp:
				self.ui.lblJendOff.setEnabled(True)

			self.ui.txtEditCend_1.setText(HCBeg)
			self.ui.txtCend_Expected.setText(HCExp)

		elif PreviewCurrentType == 'L':

			self.ui.txtFieldSearch.setPlainText(LSeqName)
			self.ui.cboFindField.setCurrentText('Name')
			done = self.on_btnFieldSearch_clicked()

			self.ui.lblSeqName.setText(LSeqName)
			self.ui.lblPrevType.setText(LGeneType)
			self.ui.txtEditVbeg_1.setText(LVBeg)
			self.ui.txtVbeg_Expected.setText(LVExp)
			if LVBeg != LVExp:
				self.ui.lblVbegOff.setEnabled(True)

			self.ui.txtEditJend_1.setText(LJEnd)
			self.ui.txtJendExp.setText(LJExp)
			if LJEnd != LJExp:
				self.ui.lblJendOff.setEnabled(True)

			self.ui.txtEditCend_1.setText(LCBeg)
			self.ui.txtCend_Expected.setText(LCExp)

		if PreviewCurrentType == 'H':

			PreviewCurrentType = 'L'
		else:
			PreviewCurrentType = 'H'

			if direction == 'forward':
				if PreviewHTcurrent < NumSeqs-1:
					PreviewHTcurrent += 1

				else:
					PreviewHTcurrent = 0
			else:
				if PreviewHTcurrent == 0:
					PreviewHTcurrent = NumSeqs-1

				else:
					PreviewHTcurrent -= 1

			# Pathname = saveFile(self, 'csv')
			#
			# with open(Pathname, 'w') as currentfile:
			# 	currentfile.write(CSVOut)

	@pyqtSlot()
	def on_btn10xEditFinal_clicked(self):

		global PreviewHTExp

		CSVOut = 'Name,Sequence\n'
		# CurrentRecord = PreviewHTExp[PreviewHTcurrent]

		if self.ui.cboAbCloningOptions.currentText() == 'HT-AbVec':
			for CurrentRecord in PreviewHTExp:
				HSeqName = CurrentRecord[0]
				Hsequence = CurrentRecord[9]
				LSeqName = CurrentRecord[10]
				Lsequence = CurrentRecord[18]
				Sequence = Lsequence + Hsequence

				CSVOut = CSVOut + HSeqName + ',' + Sequence + '\n'

		elif self.ui.cboAbCloningOptions.currentText() == 'Gibson AbVec':

			for CurrentRecord in PreviewHTExp:
				HSeqName = CurrentRecord[0]
				Hsequence = CurrentRecord[9]
				LSeqName = CurrentRecord[10]
				Lsequence = CurrentRecord[18]


				CSVOut = CSVOut + HSeqName + ',' + Hsequence + '\n' + LSeqName + ',' + Lsequence + '\n'
		elif self.ui.cboAbCloningOptions.currentText() == 'AbVec classic':

			for CurrentRecord in PreviewHTExp:
				HSeqName = CurrentRecord[0]
				Hsequence = CurrentRecord[9]
				LSeqName = CurrentRecord[10]
				Lsequence = CurrentRecord[18]


				CSVOut = CSVOut + HSeqName + ',' + Hsequence + '\n' + LSeqName + ',' + Lsequence + '\n'


		Pathname = saveFile(self, 'csv')

		with open(Pathname, 'w') as currentfile:
			currentfile.write(CSVOut)


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
			filename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'UpdateRecord.nt')
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
		datalist.append(0)
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

	def UpdateMutationAnalysis(self):
		datalist = []
		currentRow = self.ui.tableView.currentIndex().row()
		# todo change to app folder
		try:
			filename = os.path.join(os.path.expanduser('~'), 'Documents', 'Projects', 'VGenes', 'UpdateRecord.nt')
			# filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'UpdateRecord.nt')
			with open(filename, 'r') as currentfile:
				vv = currentfile
		except:
			filename = '/Volumes/Promise Pegasus/Dropbox/VGenes/UpdateRecord.nt'

		# WhichOnes = 'All'
		# if WhichOnes == 'One':
		# 	sequence = self.ui.txtDNASeq.toPlainText()
		# 	seqname = data[0]
		# 	updateseq = '>' + seqname + '\n' + sequence + '\n'
		#
		# 	with open(filename, 'w') as currentfile:
		# 		currentfile.write(updateseq)
		#
		# 	if filename == None:
		# 		return
		# 	project = self.ui.txtProject.toPlainText()
		# 	group = self.ui.txtGroup.toPlainText()
		# 	subgroup = self.ui.txtSubGroup.toPlainText()
		# 	species = data[78]
		# 	project = project.strip()
		# 	group = group.strip()
		# 	subgroup = subgroup.strip()
		# 	species = species.strip()
		#
		# 	datalist.append(project)
		# 	datalist.append(group)
		# 	datalist.append(subgroup)
		# 	datalist.append(species)
		# 	datalist.append(False)
		# # BlastIsDone = False
		#
		# elif WhichOnes == 'All':

		fields = ['SeqName', 'Sequence', 'Project', 'Grouping', 'SubGroup', 'Species', 'ID']
		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		for Seq in DataIs:
			seqname = Seq[0]
			sequence = Seq[1]
			sequence = sequence.replace('\n', '').replace('\r', '')
			sequence = sequence.upper()
			Newline = ''
			for nuc in sequence:
				if nuc == 'N' or nuc == 'A' or nuc == 'T' or nuc == 'G' or nuc == 'C':
					Newline += nuc
			updateseq = '>' + seqname + '\n' + Newline + '\n'
			with open(filename, 'w') as currentfile:
				currentfile.write(updateseq)

			project = Seq[2]
			group = Seq[3]
			subgroup = Seq[4]
			species  = Seq[5]
			DataRow = Seq[6]


			datalist.append(project)
			datalist.append(group)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(False)
			datalist.append(0)




			IgBLASTAnalysis = IgBLASTer.IgBLASTit(filename, datalist)


			i = 0
			model = self.ui.tableView.model()
			# DataRow = int(data[119])

			for record in IgBLASTAnalysis:
				for item in record:
					FieldName = FieldList[i]
					ItemValue = str(item)  # + ' '
					if i == 57 or i == 96 or i == 97 or i == 98:  #only SHM fields
						VGenesSQL.UpdateField(DataRow, ItemValue, FieldName, DBFilename)

					i += 1

			model.refresh()

		# global PreVID
		# NewID = int(data[119])+1
		# PreVID = NewID

		self.updateF(data[119])

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


if __name__ == '__main__':
	import sys

	app = QtWidgets.QApplication(sys.argv)
	Vgenes = VGenesForm()

	Vgenes.ApplicationStarted()
	sys.exit(app.exec_())
