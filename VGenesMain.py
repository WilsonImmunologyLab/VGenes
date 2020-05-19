__author__ = 'wilsonp, lei'
import os
import sys
import sqlite3 as db
import time
import re
import shutil
import math
import numpy

from PyQt5.QtCore import pyqtSlot, QTimer, Qt, QSortFilterProxyModel, pyqtSignal, QUrl, QObject, QThread, QEventLoop
from PyQt5 import QtWidgets
from PyQt5.QtPrintSupport import QPrintDialog, QPrinter
from PyQt5.QtGui import QTextCursor, QFont, QPixmap, QTextCharFormat, QBrush, QColor, QTextCursor, QCursor, QIcon
from PyQt5.QtWidgets import QApplication, QTableView, QGridLayout, QTableWidgetItem, QCheckBox, QAbstractItemView, QLabel, QLineEdit
from PyQt5.QtSql import QSqlQuery, QSqlQueryModel
from operator import itemgetter
from PyQt5.QtWebEngine import *
from PyQt5.QtWebEngineWidgets import *
from PyQt5.QtWebChannel import *
from pyecharts.charts import *
from pyecharts import options as opts
from pyecharts.globals import SymbolType
from weblogo import read_seq_data, LogoData, LogoOptions, LogoFormat, eps_formatter, svg_formatter
import itertools
from itertools import chain, groupby, zip_longest
import threading as thd
import seaborn as sns
import matplotlib
import random
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import VReports
from ui_VGenesMain import Ui_MainWindow
import IgBLASTer
from VgenesTextEdit import VGenesTextMain
import VGenesSQL
import VGenesSeq
from htmldialog import Ui_htmlDialog
from PyQt5.QtWidgets import QMainWindow

import tarfile
import zipfile
from tempfile import TemporaryDirectory

import VGenesCloneCaller
from ui_Import_Dialogue import Ui_DialogImport

global OldName
global UpdateSpecific
UpdateSpecific = True
global from_table
from_table = False
import csv
import copy
from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
	setText, openfastq
from ui_VGenesStartUpDialogue import Ui_VGenesStartUpDialog
from ui_import_data_dialog import Ui_ImportDataDialog
from ui_VGenesTextEdit import ui_TextEditor
from ui_alter_dialog import Ui_AlterDialog
from ui_annotatedialog import Ui_AnnoDialog
from ui_batch_dialog import Ui_BatchDialog
from ui_copydialog import Ui_CopyDialog
from ui_newfielddialog import Ui_NewFieldDialog
from VGenesProgressBar import ui_ProgressBar
# from VGenesPYQTSqL import EditableSqlModel, initializeModel , createConnection

global VGenesTextWindows
VGenesTextWindows = {}

from itertools import combinations
from collections import Counter
from Bio import SeqIO
from subprocess import call, Popen, PIPE

# import changeo
from changeo.IO import IMGTReader

from PyQt5.QtWidgets import QMessageBox, QInputDialog
from PyQt5.QtSql import QSqlDatabase, QSqlQuery

global LastPushed
LastPushed = ''

global working_prefix
global temp_folder
global js_folder
global clustal_path
global muscle_path
global raxml_path

working_prefix = os.path.dirname(os.path.realpath(sys.argv[0])) + '/'
temp_folder = os.path.join(working_prefix, 'Temp')
js_folder = os.path.join(working_prefix, 'JS')
clustal_path = os.path.join(working_prefix, 'Tools', 'clustalo')
muscle_path = os.path.join(working_prefix, 'Tools', 'muscle')
raxml_path = os.path.join(working_prefix, 'Tools', 'raxml')

ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
r = open(ErlogFile,'w')
r.write('')
r.close()
r = open(ErlogFile2, 'w')
r.write('')
r.close()

global IgBLASTAnalysis
IgBLASTAnalysis = []

global IMGTAnalysis
Analysis = []

global TreeSelected
TreeSelected = []

global timer
global data
global LastSelected
global DBFilename
global wasClicked
global LineLen
LastSelected = ()
global updateMarker
updateMarker = False
global DontFindTwice
DontFindTwice = False
wasClicked = False
data = []
global MoveNotChange
MoveNotChange = False
global PreviewHTExp
PreviewHTExp = []
global PreviewHTcurrent
PreviewHTcurrent = 0
global PreviewCurrentType
PreviewCurrentType = 'H'
global JustMovedIt
JustMovedIt = True
global JustMoved
JustMoved = False


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

global FieldTypeList
FieldTypeList = ["Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed",
                 "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed", "Fixed"]
global FieldCommentList
FieldCommentList = ["", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", "", "", "","", "", "",
                    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "",
                    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", "",
                    "", "", "", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", "", "", "", "",
                    "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", "", "", "", ""]

global StopCheckProgress
StopCheckProgress = True
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

def reName(ori_name, rep1, rep2, prefix):
	my_name = re.sub('clonotype',rep1, ori_name)
	my_name = re.sub('consensus_', rep2, my_name)
	if prefix != '':
		my_name = prefix + '_' + my_name
	return my_name

class MyFigure(FigureCanvas):
    def __init__(self,width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MyFigure,self).__init__(self.fig)
        self.axes = self.fig.add_subplot(111)

class NewFieldDialog(QtWidgets.QDialog, Ui_NewFieldDialog):
	NewFieldSignal = pyqtSignal(str, str)

	def __init__(self):
		super(NewFieldDialog, self).__init__()
		self.ui = Ui_NewFieldDialog()
		self.ui.setupUi(self)

		self.ui.labelTip.setHidden(True)
		self.ui.DisplayTip.setHidden(True)

		self.ui.Cancel.clicked.connect(self.reject)
		self.ui.OK.clicked.connect(self.accept)
		self.ui.comboBoxFrom.currentTextChanged.connect(self.StatFig)
		self.ui.radioButton.clicked.connect(self.StatFig)
		self.ui.DisplayTip.clicked.connect(self.StatFig)

	def accept(self):
		field_from = self.ui.comboBoxFrom.currentText()
		new_field = self.ui.lineEdit.text()

		if new_field == '':
			Msg = 'This field name can not be empty!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		if new_field in FieldList:
			Msg = 'This field name has been taken! Enter another field name!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		else:
			self.NewFieldSignal.emit(field_from, new_field)
			self.close()

	def StatFig(self):
		if DontFindTwice == True:
			print('skip')
			return

		sender_widget = self.sender()

		if self.ui.gridLayoutFig.count() > 0:
			for i in range(self.ui.gridLayoutFig.count()):
				self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()
		else:
			self.resize(800 + random.randint(0,9), 600 + random.randint(0,9))
		# numeric value
		try:
			if self.ui.radioButton.isChecked():
				field = re.sub(r'\(.+', '', self.ui.comboBoxFrom.currentText())
				if field == '':
					return

				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				value_list = []
				char_list = []
				non_number_count = 0
				for row in DataIn:
					try:
						value_list.append(float(row[0]))
					except:
						char_list.append(row[0])
						non_number_count += 1

				if len(value_list) == 0:
					Msg = 'No value can be converted to number!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					self.ui.radioButton.setChecked(False)
					self.StatFig()
					return

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.hist(value_list, bins=30)

				F.axes.tick_params(labelsize=7)
				F.fig.subplots_adjust(bottom=0.1)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)
				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
			# character value
			else:
				field = re.sub(r'\(.+', '', self.ui.comboBoxFrom.currentText())
				if field == '':
					return

				SQLStatement = 'SELECT DISTINCT(' + field + ') FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				value_list = [row[0] for row in DataIn]

				if sender_widget.objectName() == 'DisplayTip':
					pass
				else:
					if len(value_list) > 30:
						if self.ui.gridLayoutFig.count() > 0:
							for i in range(self.ui.gridLayoutFig.count()):
								self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()

						self.ui.labelTip.setHidden(False)
						self.ui.DisplayTip.setHidden(False)
						return

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				data = []
				for element in DataIn:
					data.append(element[0])
				result = Counter(data)
				labels = result.keys()
				values = result.values()
				colors = sns.color_palette("hls", len(values))

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.bar(labels, values, color=colors)
				F.axes.set_xticklabels(labels, rotation=-90)
				F.axes.tick_params(labelsize=7)

				# determine spacing
				lens = [len(lab) for lab in labels]
				max_len = max(lens)
				my_adjust = 0.1 + max_len / 50
				F.fig.subplots_adjust(bottom=my_adjust)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)
				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
		except:
			print('error')
			return

class CopyDialog(QtWidgets.QDialog, Ui_CopyDialog):
	CopySignal = pyqtSignal(str, str)

	def __init__(self):
		super(CopyDialog, self).__init__()
		self.ui = Ui_CopyDialog()
		self.ui.setupUi(self)
		self.ui.labelTip.setHidden(True)
		self.ui.DisplayTip.setHidden(True)

		self.ui.Cancel.clicked.connect(self.reject)
		self.ui.OK.clicked.connect(self.accept)
		self.ui.comboBoxFrom.currentTextChanged.connect(self.StatFig)
		self.ui.radioButton.clicked.connect(self.StatFig)
		self.ui.DisplayTip.clicked.connect(self.StatFig)

	def accept(self):
		field_from = self.ui.comboBoxFrom.currentText()
		field_to = self.ui.comboBoxTo.currentText()

		if field_from == '' or field_to == '':
			Msg = 'Both field names can not be empty!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		elif field_from == field_to:
			Msg = 'Two field names are the same!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		else:
			self.CopySignal.emit(field_from, field_to)
			self.close()

	def StatFig(self):
		sender_widget = self.sender()
		# numeric value
		if self.ui.radioButton.isChecked():
			field = re.sub(r'\(.+', '', self.ui.comboBoxFrom.currentText())
			if field == '':
				return
			
			SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			value_list = []
			char_list = []
			non_number_count = 0
			for row in DataIn:
				try:
					value_list.append(float(row[0]))
				except:
					char_list.append(row[0])
					non_number_count += 1

			if len(value_list) == 0:
				Msg = 'No value can be converted to number!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				self.ui.radioButton.setChecked(False)
				self.StatFig()
				return

			# update figure
			SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			F = MyFigure(width=3, height=3, dpi=160)
			F.axes.hist(value_list, bins=30)

			F.axes.tick_params(labelsize=7)
			F.fig.subplots_adjust(bottom=0.1)

			if self.ui.gridLayoutFig.count() > 0:
				for i in range(self.ui.gridLayoutFig.count()):
					self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()
			else:
				self.resize(1200 + random.randint(0,9), 700 + random.randint(0,9))
			self.ui.gridLayoutFig.addWidget(F, 0, 1)
			self.ui.labelTip.setHidden(True)
			self.ui.DisplayTip.setHidden(True)
		# character value
		else:
			field = re.sub(r'\(.+', '', self.ui.comboBoxFrom.currentText())
			if field == '':
				return
		
			SQLStatement = 'SELECT DISTINCT(' + field + ') FROM vgenesdb'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			value_list = [row[0] for row in DataIn]

			if sender_widget.objectName() == 'DisplayTip':
				pass
			else:
				if len(value_list) > 30:
					if self.ui.gridLayoutFig.count() > 0:
						for i in range(self.ui.gridLayoutFig.count()):
							self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()

					self.ui.labelTip.setHidden(False)
					self.ui.DisplayTip.setHidden(False)
					return

			# update figure
			SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			data = []
			for element in DataIn:
				data.append(element[0])
			result = Counter(data)
			labels = result.keys()
			values = result.values()
			colors = sns.color_palette("hls", len(values))

			F = MyFigure(width=3, height=3, dpi=160)
			F.axes.bar(labels, values, color=colors)
			F.axes.set_xticklabels(labels, rotation=-90)
			F.axes.tick_params(labelsize=7)

			# determine spacing
			lens = [len(lab) for lab in labels]
			max_len = max(lens)
			my_adjust = 0.1 + max_len / 50
			F.fig.subplots_adjust(bottom=my_adjust)

			if self.ui.gridLayoutFig.count() > 0:
				for i in range(self.ui.gridLayoutFig.count()):
					self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()
			else:
				self.resize(1200 + random.randint(0,9), 700 + random.randint(0,9))
			self.ui.gridLayoutFig.addWidget(F, 0, 1)
			self.ui.labelTip.setHidden(True)
			self.ui.DisplayTip.setHidden(True)

class BatchDialog(QtWidgets.QDialog, Ui_BatchDialog):
	BatchSignal = pyqtSignal(int, str, dict)

	def __init__(self):
		super(BatchDialog, self).__init__()
		self.ui = Ui_BatchDialog()
		self.ui.setupUi(self)
		self.initial = 0
		self.ui.LineEditCutoff.setHidden(True)
		self.ui.labelTip.setHidden(True)
		self.ui.DisplayTip.setHidden(True)
		self.ui.LineEditCutoff.min = 0
		self.ui.LineEditCutoff.max = 0
		self.ui.gridLayout.num_widget = 0
		self.ui.gridLayoutChar.num_widget = 0


		self.ui.pushButtonCancel.clicked.connect(self.reject)
		self.ui.pushButtonOK.clicked.connect(self.accept)
		self.ui.comboBox.currentTextChanged.connect(self.StatFig)
		self.ui.radioButton.clicked.connect(self.StatFig)
		self.ui.LineEditCutoff.textChanged.connect(self.StatFig)
		self.ui.DisplayTip.clicked.connect(self.StatFig)

	def updateNum(self):
		if self.ui.radioButton.isChecked():
			num_list = []

			if self.ui.LineEditCutoff.text() == '':
				self.load_data_num([], self.ui.LineEditCutoff.min, self.ui.LineEditCutoff.max)
			elif self.ui.LineEditCutoff.text()[0] == 'T':
				self.load_data_num([], self.ui.LineEditCutoff.min, self.ui.LineEditCutoff.max)
			else:
				temp_data = self.ui.LineEditCutoff.text().split(',')
				if len(temp_data) > 0:
					for ele in temp_data:
						try:
							num = float(ele)
							if num > self.ui.LineEditCutoff.min and num < self.ui.LineEditCutoff.max:
								num_list.append(num)
						except:
							return
					# remove redudant and sort
					num_list = list(set(num_list))
					num_list.sort()

					self.load_data_num(num_list, self.ui.LineEditCutoff.min, self.ui.LineEditCutoff.max)
		else:
			self.ui.LineEditCutoff.setHidden(True)

	def StatFig(self):
		sender_widget = self.sender()

		if self.ui.gridLayoutFig.count() > 0:
			for i in range(self.ui.gridLayoutFig.count()):
				self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()
		try:
			# numeric value
			if self.ui.radioButton.isChecked():
				try:
					sender = self.sender()
					if sender.objectName() == 'comboBox':
						self.ui.LineEditCutoff.setText('Type cutoff here, seprate by ,  (e.g.  500,600,700)')
					elif sender.objectName() == 'radioButton':
						self.ui.LineEditCutoff.setText('Type cutoff here, seprate by ,  (e.g.  500,600,700)')
				except:
					pass

				self.ui.LineEditCutoff.setHidden(False)
				field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				value_list = []
				char_list = []
				non_number_count = 0
				for row in DataIn:
					try:
						value_list.append(float(row[0]))
					except:
						char_list.append(row[0])
						non_number_count += 1

				if len(value_list) == 0:
					Msg = 'No value can be converted to number!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					self.ui.radioButton.setChecked(False)
					self.StatFig()
					return

				# do it later
				## number list
				self.ui.LineEditCutoff.min = min(value_list)
				self.ui.LineEditCutoff.max = max(value_list)
				self.updateNum()

				## char list
				char_list = list(set(char_list))
				self.load_data_char(char_list)

				num_list = []
				error = False
				if self.ui.LineEditCutoff.text() != '':
					temp_data = self.ui.LineEditCutoff.text().split(',')
					if len(temp_data) > 0:
						for ele in temp_data:
							try:
								num = float(ele)
								if num > self.ui.LineEditCutoff.min and num < self.ui.LineEditCutoff.max:
									num_list.append(num)
							except:
								error = True
						# remove redudant and sort
						num_list = list(set(num_list))
						num_list.sort()

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.hist(value_list, bins=30)

				if error == False:
					ymin, ymax = F.axes.get_ylim()
					for num in num_list:
						F.axes.plot([num,num], [ymin, ymax], color='r', linewidth = 1, label="Cutoff")

				F.axes.tick_params(labelsize=7)
				F.fig.subplots_adjust(bottom=0.1)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)

				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
			# character value
			else:
				self.ui.LineEditCutoff.setHidden(True)
				if self.initial == 0:
					return
				elif self.initial == 1:
					field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
					SQLStatement = 'SELECT DISTINCT(' + field + ') FROM vgenesdb'
					DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
					value_list = [row[0] for row in DataIn]
				elif self.initial == 2:
					field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
					SQLStatement = 'SELECT DISTINCT(' + field + ') FROM vgenesdb'
					DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
					value_list = [row[0] for row in DataIn]

					if sender_widget.objectName() == 'DisplayTip':
						pass
					else:
						if len(value_list) > 30:
							question = 'Distinct values of this field seems too many (number =  ' + str(
								len(value_list)) + ')\nAre you sure?'
							buttons = 'YN'
							answer = questionMessage(self, question, buttons)
							if answer == 'No':
								self.ui.labelTip.setHidden(False)
								self.ui.DisplayTip.setHidden(False)
								return

				self.load_data(value_list)

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				data = []
				for element in DataIn:
					data.append(element[0])
				result = Counter(data)
				labels = result.keys()
				values = result.values()
				colors = sns.color_palette("hls", len(values))

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.bar(labels, values, color=colors)
				F.axes.set_xticklabels(labels, rotation=-90)
				F.axes.tick_params(labelsize=7)

				# determine spacing
				lens = [len(lab) for lab in labels]
				max_len = max(lens)
				my_adjust = 0.1 + max_len/50
				F.fig.subplots_adjust(bottom=my_adjust)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)

				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
		except:
			return

	def accept(self):
		global DontFindTwice

		if self.ui.radioButton.isChecked():
			new_values = []
			# process num part
			layout = self.ui.gridLayout
			i = 1
			#print('row count = ' + str(layout.rowCount()))
			while i < layout.num_widget:
				str2 = layout.itemAtPosition(i, 1).widget().text()
				new_values.append(str2)
				i += 1

			num_list = [self.ui.LineEditCutoff.min, self.ui.LineEditCutoff.max]
			if self.ui.LineEditCutoff.text() == '':
				pass
			else:
				temp_data = self.ui.LineEditCutoff.text().split(',')
				if len(temp_data) > 0:
					for ele in temp_data:
						try:
							num = float(ele)
							if num > self.ui.LineEditCutoff.min and num < self.ui.LineEditCutoff.max:
								num_list.append(num)
						except:
							return
					# remove redudant and sort
					num_list = list(set(num_list))
					num_list.sort()

			# start update records
			field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
			SQLStatement = 'SELECT ' + field + ',ID FROM vgenesdb'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for ele in DataIn:
				try:
					cur_value = float(ele[0])
					cur_id = str(ele[1])
					for i in range(len(num_list)-1):
						if cur_value >= num_list[i] and cur_value < num_list[i+1]:
							new_value = new_values[i]
							SQLStatement = 'UPDATE vgenesdb SET ' + field + '= "' + new_value + '" WHERE ID=' + cur_id
							VGenesSQL.RunUpdateSQL(DBFilename, SQLStatement)
				except:
					pass
			# process char part
			Dict = {}
			layout = self.ui.gridLayoutChar
			i = 1
			while i < layout.num_widget:
				str1 = layout.itemAtPosition(i, 0).widget().text()
				str2 = layout.itemAtPosition(i, 1).widget().text()
				if str2 == '':
					pass
				else:
					Dict[str1] = str2
				i += 1

			if len(Dict) > 0:
				field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
				self.BatchSignal.emit(0, field, Dict)
			else:
				DontFindTwice = True
				Vgenes.refreshDB()
				value = Vgenes.ui.dial.value()
				Vgenes.updateF(value)
				DontFindTwice = False

				Msg = 'Update finished!'
				QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)
		else:
			Dict = {}
			layout = self.ui.gridLayout
			i = 1
			while i < layout.num_widget:
				str1 = layout.itemAtPosition(i, 0).widget().text()
				str2 = layout.itemAtPosition(i, 1).widget().text()
				if str2 == '':
					pass
				else:
					Dict[str1] = str2
				i += 1

			if len(Dict) > 0:
				field = re.sub(r'\(.+', '', self.ui.comboBox.currentText())
				self.BatchSignal.emit(0, field, Dict)
			else:
				self.BatchSignal.emit(1, '', Dict)

		self.close()

	def load_data(self, list):
		layout = self.ui.gridLayout
		if layout.count() > 0:
			for i in range(layout.count()):
				layout.itemAt(i).widget().deleteLater()
		layout.num_widget = 0

		layout.addWidget(QLabel("Original value"),0,0)
		layout.addWidget(QLabel("New value"), 0, 1)
		layout.num_widget += 1

		i = 1
		for item in list:
			f = QLineEdit(item)
			f.setReadOnly(True)
			layout.addWidget(f, i, 0)
			layout.addWidget(QLineEdit(""), i, 1)
			layout.num_widget += 1
			i += 1

		# delete everything in char layout
		layout = self.ui.gridLayoutChar
		if layout.count() > 0:
			for i in range(layout.count()):
				layout.itemAt(i).widget().deleteLater()

	def load_data_num(self, list, min, max):
		# clear old widgets
		layout = self.ui.gridLayout
		if layout.count() > 0:
			for i in range(layout.count()):
				layout.itemAt(i).widget().deleteLater()
		layout.num_widget = 0

		if len(list) == 0:
			layout.addWidget(QLabel("Data range:"), 0, 0)
			layout.addWidget(QLabel("New value"), 0, 1)
			layout.num_widget += 1

			item = str(min) + ' <= Value <= ' + str(max)
			f = QLineEdit(item)
			f.setReadOnly(True)
			layout.addWidget(f, 1, 0)
			layout.addWidget(QLineEdit(""), 1, 1)
			layout.num_widget += 1
		else:
			layout.addWidget(QLabel("Data range:"), 0, 0)
			layout.addWidget(QLabel("New value"), 0, 1)
			layout.num_widget += 1

			for i in range(len(list)):
				if i == 0:
					cur_range = str(min) + ' <= Value < ' + str(list[i])
				else:
					cur_range = str(list[i-1]) + ' <= Value < ' + str(list[i])

				f = QLineEdit(cur_range)
				f.setReadOnly(True)
				layout.addWidget(f, i+1, 0)
				layout.addWidget(QLineEdit(""), i+1, 1)
				layout.num_widget += 1

			cur_range = str(list[-1]) + ' <= Value <= ' + str(max)
			f = QLineEdit(cur_range)
			f.setReadOnly(True)
			layout.addWidget(f, len(list) + 1, 0)
			layout.addWidget(QLineEdit(""), len(list) + 1, 1)
			layout.num_widget += 1

	def load_data_char(self, list):
		layout = self.ui.gridLayoutChar
		if layout.count() > 0:
			for i in range(layout.count()):
				layout.itemAt(i).widget().deleteLater()
		layout.num_widget = 0
		
		if len(list) > 0:
			layout.addWidget(QLabel("Original value"),0,0)
			layout.addWidget(QLabel("New value"), 0, 1)
			layout.num_widget += 1
	
			i = 1
			for item in list:
				f = QLineEdit(item)
				f.setReadOnly(True)
				layout.addWidget(f, i, 0)
				layout.addWidget(QLineEdit(""), i, 1)
				layout.num_widget += 1
				i += 1

class htmlDialog(QtWidgets.QDialog):
	def __init__(self):
		super(htmlDialog, self).__init__()
		self.ui = Ui_htmlDialog()
		self.ui.setupUi(self)

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

class VGenesTextMain(QtWidgets.QMainWindow, ui_TextEditor):
	def __init__(self, parent=None):
		QtWidgets.QMainWindow.__init__(self, parent)
		# super(VGenesTextMain, self).__init__()
		self.setupUi()

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
			filename = os.path.join(working_prefix, 'Conf', 'RecentPaths.vtx')

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

class AnnoDielog(QtWidgets.QDialog, Ui_AnnoDialog):
	refreshDBSignal = pyqtSignal()

	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		#super(ImportDataDialogue, self).__init__()
		self.ui = Ui_AnnoDialog()
		self.ui.setupUi(self)

		self.ui.pushButtonCancel.clicked.connect(self.reject)
		self.ui.radioButton.clicked.connect(self.switchHeader)
		self.ui.pushButtonOK.clicked.connect(self.accept)

	def accept(self):
		global RealNameList
		global FieldCommentList
		global FieldTypeList
		global FieldList

		col_index = []
		col_fields = []

		header = self.header
		anchor_field = self.ui.comboBox.currentText()
		target_field = self.ui.comboBox2.currentText()
		target_field = re.sub(r'\s.+','',target_field)
		anchor_col_index = header.index(anchor_field)

		num_col = self.ui.tableWidget.columnCount()
		num_row = self.ui.tableWidget.rowCount()
		for col in range(num_col):
			if col == anchor_col_index:
				continue
			my_widget = self.ui.tableWidget.cellWidget(0, col)
			if my_widget.currentText() != "":
				col_index.append(col)
				col_fields.append(my_widget.currentText())

				# add new column if field name not exit in current col
				tmp_field_name = re.sub(r'\s.+', '', my_widget.currentText())
				if tmp_field_name in FieldList:
					if self.ui.radioButtonUpdateName.isChecked():
						VGenesSQL.UpdateFieldTable(tmp_field_name, header[col], 'FieldNickName', DBFilename)
						RealNameList[FieldList.index(tmp_field_name)] = header[col]
				else:
					# check if the new field name can be used:
					HEADERStatement = 'PRAGMA table_info(vgenesDB);'
					HeaderIn = VGenesSQL.RunSQL(DBFilename, HEADERStatement)
					ALL_Fields = [i[1] for i in HeaderIn]

					if tmp_field_name in ALL_Fields:
						pass
					else:
						# update vgene table
						SQLSTATEMENT1 = "ALTER TABLE vgenesDB ADD " + tmp_field_name + " text"

						try:
							VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
						except:
							msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
							QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
							                    QMessageBox.Ok)
							return

					# update field name table
					SQLSTATEMENT2 = 'INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment) ' \
					                'VALUES(' + str(len(FieldList) + 1) + ',"' + tmp_field_name + '", "' + \
					                header[col] + '", "Customized", "")'
					try:
						VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
					except:
						msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
						QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
						                    QMessageBox.Ok)
						return

					RealNameList.append(header[col])
					FieldCommentList.append('')
					FieldTypeList.append('Customized')
					FieldList.append(tmp_field_name)

		count = 0
		for row in range(num_row):
			if row == 0:
				continue
			else:
				SQLSTATEMENT = "UPDATE vgenesdb SET "
				for i in range(len(col_index)):
					col = col_index[i]
					field = col_fields[i]
					field = re.sub(r'\s.+','',field)
					value = self.ui.tableWidget.item(row, col).text()
					SQLSTATEMENT = SQLSTATEMENT + field + ' = "' + value + '",'
				current_anchor = self.ui.tableWidget.item(row, anchor_col_index).text()
				SQLSTATEMENT = SQLSTATEMENT.rstrip(',')
				SQLSTATEMENT = SQLSTATEMENT + " WHERE " + target_field + ' = "' + current_anchor + '"'
				try:
					count += VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT)
					#print(SQLSTATEMENT)
				except:
					Msg = 'SQL error! Current SQL statement is:\n' + SQLSTATEMENT
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok,
					                        QMessageBox.Ok)
					return
				#print(SQLSTATEMENT)
		if count == 0:
			Msg = "Total " + str(count) + ' records affected! Maybe check if you anchor column is correct?'
			QMessageBox.information(self, 'information', Msg, QMessageBox.Ok, QMessageBox.Ok)

		else:
			Msg = "Update successfully!\nTotal " + str(count) + ' records affected!'
			QMessageBox.information(self, 'information', Msg, QMessageBox.Ok, QMessageBox.Ok)

			self.refreshDBSignal.emit()
			self.hide()


	def reject(self):
		self.hide()

	def switchHeader(self):
		if self.ui.radioButton.isChecked():
			Content = list(self.Content)
			horizontalHeader = Content.pop(0)
			self.ui.comboBox.clear()
			self.ui.comboBox.addItems(horizontalHeader)

			HEADERStatement = 'SELECT Field, FieldNickName FROM fieldsname'
			HeaderIn = VGenesSQL.RunSQL(DBFilename, HEADERStatement)
			Fields = [i[0] + '  (' + i[1] + ')' for i in HeaderIn]
			self.ui.comboBox2.addItems(Fields)

			num_col = len(horizontalHeader)
			num_row = len(Content)

			self.ui.tableWidget.setRowCount(num_row)
			self.ui.tableWidget.setColumnCount(num_col)
			self.ui.tableWidget.setHorizontalHeaderLabels(horizontalHeader)
			self.header = horizontalHeader

			for col_index in range(num_col):
				cell_comBox = QtWidgets.QComboBox()
				cell_comBox.addItems([''] + Fields)
				cell_comBox.setMaximumSize(10086, 40)
				cell_comBox.setMinimumSize(50, 20)
				cell_comBox.setEditable(True)
				self.ui.tableWidget.setCellWidget(0, col_index, cell_comBox)

			for row_index in range(num_row):
				for col_index in range(num_col):
					self.ui.tableWidget.setItem(row_index + 1, col_index,
					                                       QTableWidgetItem(Content[row_index][col_index]))
		else:
			Content = list(self.Content)
			horizontalHeader = list(Content[0])
			for i in range(len(horizontalHeader)):
				horizontalHeader[i] = 'Column' + str(i+1)

			self.ui.comboBox.clear()
			self.ui.comboBox.addItems(horizontalHeader)

			HEADERStatement = 'SELECT Field, FieldNickName FROM fieldsname'
			HeaderIn = VGenesSQL.RunSQL(DBFilename, HEADERStatement)
			Fields = [i[0] + '  (' + i[1] + ')' for i in HeaderIn]
			self.ui.comboBox2.addItems(Fields)

			num_col = len(horizontalHeader)
			num_row = len(Content)

			self.ui.tableWidget.setRowCount(num_row)
			self.ui.tableWidget.setColumnCount(num_col)
			self.ui.tableWidget.setHorizontalHeaderLabels(horizontalHeader)
			self.header = horizontalHeader

			for col_index in range(num_col):
				cell_comBox = QtWidgets.QComboBox()
				cell_comBox.addItems([''] + Fields)
				cell_comBox.setMaximumSize(10086, 40)
				cell_comBox.setMinimumSize(50, 20)
				cell_comBox.setEditable(True)
				self.ui.tableWidget.setCellWidget(0, col_index, cell_comBox)

			for row_index in range(num_row):
				for col_index in range(num_col):
					self.ui.tableWidget.setItem(row_index + 1, col_index,
					                                       QTableWidgetItem(Content[row_index][col_index]))

class AlterDielog(QtWidgets.QDialog, Ui_AlterDialog):
	refreshDBSignal = pyqtSignal()

	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		#super(ImportDataDialogue, self).__init__()
		self.ui = Ui_AlterDialog()
		self.ui.setupUi(self)

		self.ui.labelTip.setHidden(True)
		self.ui.DisplayTip.setHidden(True)

		self.ui.pushButton.clicked.connect(self.reject)
		self.ui.pushButtonSave.clicked.connect(self.saveRecord)
		self.ui.lineEditNickName.textChanged.connect(self.valueChange)
		self.ui.textEditNote.textChanged.connect(self.valueChange)
		self.ui.listWidget.itemSelectionChanged.connect(self.ListItemChanged)
		self.ui.pushButtonNew.clicked.connect(self.addField)
		self.ui.pushButtonDelete.clicked.connect(self.deleteField)
		self.ui.lineEditNickName.setReadOnly(False)
		self.ui.pushButtonSave.setEnabled(False)
		self.ui.lineEditName.textChanged.connect(self.StatFig)
		self.ui.radioButton.clicked.connect(self.StatFig)
		self.ui.DisplayTip.clicked.connect(self.StatFig)

	def StatFig(self):
		sender_widget = self.sender()
		if self.ui.gridLayoutFig.count() > 0:
			for i in range(self.ui.gridLayoutFig.count()):
				self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()
		# numeric value
		try:
			if self.ui.radioButton.isChecked():
				field = re.sub(r'\(.+', '', self.ui.lineEditName.text())
				if field == '':
					return

				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				value_list = []
				char_list = []
				non_number_count = 0
				for row in DataIn:
					try:
						value_list.append(float(row[0]))
					except:
						char_list.append(row[0])
						non_number_count += 1

				if len(value_list) == 0:
					Msg = 'No value can be converted to number!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					self.ui.radioButton.setChecked(False)
					self.StatFig()
					return

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.hist(value_list, bins=30)

				F.axes.tick_params(labelsize=7)
				F.fig.subplots_adjust(bottom=0.1)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)
				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
			# character value
			else:
				field = re.sub(r'\(.+', '', self.ui.lineEditName.text())
				if field == '':
					return

				SQLStatement = 'SELECT DISTINCT(' + field + ') FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				value_list = [row[0] for row in DataIn]

				if sender_widget.objectName() == 'DisplayTip':
					pass
				else:
					if len(value_list) > 30:
						if self.ui.gridLayoutFig.count() > 0:
							for i in range(self.ui.gridLayoutFig.count()):
								self.ui.gridLayoutFig.itemAt(i).widget().deleteLater()

						self.ui.labelTip.setHidden(False)
						self.ui.DisplayTip.setHidden(False)
						return

				# update figure
				SQLStatement = 'SELECT ' + field + ' FROM vgenesdb'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				data = []
				for element in DataIn:
					data.append(element[0])
				result = Counter(data)
				labels = result.keys()
				values = result.values()
				colors = sns.color_palette("hls", len(values))

				F = MyFigure(width=3, height=3, dpi=160)
				F.axes.bar(labels, values, color=colors)
				F.axes.set_xticklabels(labels, rotation=-90)
				F.axes.tick_params(labelsize=7)

				# determine spacing
				lens = [len(lab) for lab in labels]
				max_len = max(lens)
				my_adjust = 0.1 + max_len / 50
				F.fig.subplots_adjust(bottom=my_adjust)

				self.ui.gridLayoutFig.addWidget(F, 0, 1)
				self.ui.labelTip.setHidden(True)
				self.ui.DisplayTip.setHidden(True)
		except:
			return

	def deleteField(self):
		global FieldList
		global RealNameList
		global FieldCommentList
		global FieldTypeList

		target_field = self.ui.listWidget.currentItem().text()
		if self.ui.lineEditType == "Fixed":
			Msg = 'You can not delete Fixed filed!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		Msg = "Are you sure to delete field: " + target_field + " ?"
		reply = QMessageBox.question(self, 'Information', Msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
		if reply == QMessageBox.No:
			return

		# update table
		#SQLSTATEMENT1 = "ALTER TABLE vgenesDB DROP COLUMN " + target_field
		#try:
		#	VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
		#except:
		#	msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
		#	QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
		#	                    QMessageBox.Ok)
		#	return

		SQLSTATEMENT2 = 'DELETE FROM fieldsname WHERE Field = "' + target_field + '"'
		try:
			VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
		except:
			msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
			                    QMessageBox.Ok)
			return



		index = FieldList.index(target_field)
		RealNameList = RealNameList[0:index] + RealNameList[index + 1:]
		FieldCommentList = FieldCommentList[0:index] + FieldCommentList[index + 1:]
		FieldTypeList = FieldTypeList[0:index] + FieldTypeList[index + 1:]
		FieldList = FieldList[0:index] + FieldList[index + 1:]

		# update widget
		self.ui.listWidget.clear()
		self.ui.listWidget.addItems(FieldList)
		self.ui.listWidget.setCurrentItem(self.ui.listWidget.item(0))
		self.ui.listWidget.scrollToTop()

		self.refreshDBSignal.emit()

	def addField(self):
		global DontFindTwice
		self.newFieldDialog = NewFieldDialog()
		field_list = [''] + [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		DontFindTwice = True
		self.newFieldDialog.ui.comboBoxFrom.addItems(field_list)
		DontFindTwice = False
		self.newFieldDialog.NewFieldSignal.connect(self.addFieldAndCopy)
		self.newFieldDialog.exec_()


	def addFieldAndCopy(self, field_from, new_field):
		print(field_from)
		print(new_field)
		global RealNameList
		global FieldCommentList
		global FieldTypeList
		global FieldList

		# update vgene table
		SQLSTATEMENT1 = "ALTER TABLE vgenesDB ADD " + new_field + " text"

		try:
			VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
		except:
			msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

		# update field name table
		SQLSTATEMENT2 = 'INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment) ' \
		                'VALUES(' + str(
			self.ui.listWidget.count() + 1) + ',"' + new_field + '", "' + new_field + '", "Customized", "")'
		try:
			VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
		except:
			msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

		# copy record
		if field_from != '':
			field_from = re.sub(r'\(.+', '', field_from)
			SQLStatement = 'UPDATE vgenesDB SET ' + new_field + ' = ' + field_from
			try:
				VGenesSQL.RunUpdateSQL(DBFilename, SQLStatement)
			except:
				Msg = 'Error occurs when updating the DB!\nCurrent SQL statement is:\n' + SQLStatement
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

		RealNameList.append(new_field)
		FieldCommentList.append('')
		FieldTypeList.append('Customized')
		FieldList.append(new_field)

		self.ui.listWidget.addItem(new_field)
		self.ui.listWidget.setCurrentItem(self.ui.listWidget.item(self.ui.listWidget.count() - 1))
		self.ui.listWidget.scrollToBottom()

		self.refreshDBSignal.emit()

	def addField1(self):
		global RealNameList
		global FieldCommentList
		global FieldTypeList
		global FieldList

		items = []
		for index in range(self.ui.listWidget.count()):
			items.append(self.ui.listWidget.item(index))
		cur_fields = [i.text() for i in items]

		tmp_field_name, ok = QInputDialog.getText(self, 'Input Dialog','Enter name of your new field:')
		while tmp_field_name in cur_fields and ok:
			tmp_field_name, ok = QInputDialog.getText(self, 'Input Dialog','This field name has been taken! '
			                                                               'Enter another name of your new field:')

		if ok:
			pass
		else:
			return

		# check if the new field name can be used:
		HEADERStatement = 'PRAGMA table_info(vgenesDB);'
		HeaderIn = VGenesSQL.RunSQL(DBFilename, HEADERStatement)
		ALL_Fields = [i[1] for i in HeaderIn]

		if tmp_field_name in ALL_Fields:
			pass
		else:
			# update vgene table
			SQLSTATEMENT1 = "ALTER TABLE vgenesDB ADD " + tmp_field_name + " text"

			try:
				VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
			except:
				msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
				QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
				                    QMessageBox.Ok)
				return

		# update field name table
		SQLSTATEMENT2 = 'INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment) ' \
		                'VALUES(' + str(
			self.ui.listWidget.count() + 1) + ',"' + tmp_field_name + '", "' + tmp_field_name + '", "Customized", "")'
		try:
			VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
		except:
			msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

		RealNameList.append(tmp_field_name)
		FieldCommentList.append('')
		FieldTypeList.append('Customized')
		FieldList.append(tmp_field_name)

		self.ui.listWidget.addItem(tmp_field_name)
		self.ui.listWidget.setCurrentItem(self.ui.listWidget.item(self.ui.listWidget.count()-1))
		self.ui.listWidget.scrollToBottom()

		self.refreshDBSignal.emit()

	def saveRecord(self):
		global RealNameList
		global FieldCommentList
		global FieldTypeList
		global FieldList

		items = self.ui.listWidget.selectedItems()
		for item in items:
			name = item.text()

		index = FieldList.index(name)
		RealNameList[index] = self.ui.lineEditNickName.text()
		FieldCommentList[index] = self.ui.textEditNote.toPlainText()

		# save to DB
		VGenesSQL.UpdateFieldTable(name, self.ui.textEditNote.toPlainText(), 'FieldComment', DBFilename)
		VGenesSQL.UpdateFieldTable(name, self.ui.lineEditNickName.text(), 'FieldNickName', DBFilename)

		self.ui.pushButtonSave.setEnabled(False)
		self.refreshDBSignal.emit()

	def valueChange(self):
		self.ui.pushButtonSave.setEnabled(True)

	def ListItemChanged(self):
		items = self.ui.listWidget.selectedItems()
		if len(items) == 0:
			return
		for item in items:
			name = item.text()

		index = FieldList.index(name)
		self.ui.lineEditName.setText(FieldList[index])
		self.ui.lineEditNickName.setText(RealNameList[index])
		self.ui.lineEditType.setText(FieldTypeList[index])
		self.ui.textEditNote.setText(FieldCommentList[index])
		self.ui.pushButtonSave.setEnabled(False)
		self.ui.lineEditNickName.setReadOnly(False)

	def reject(self):
		self.hide()

class ImportDataDialogue(QtWidgets.QDialog, Ui_DialogImport):
	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		#super(ImportDataDialogue, self).__init__()
		self.ui = Ui_ImportDataDialog()
		self.ui.setupUi(self)

		self.TextEdit = VGenesTextMain()

		self.ui.pushButtonOK.clicked.connect(self.accept)
		self.ui.pushButtonCancel.clicked.connect(self.reject)
		self.ui.tabWidget.currentChanged['int'].connect(self.switchTab)
		self.ui.browse10x.clicked.connect(self.browse10x)
		self.ui.browseCSV.clicked.connect(self.browseCSV)
		self.ui.browseFasta.clicked.connect(self.browseFasta)
		self.ui.browseVDB.clicked.connect(self.browseVDB)
		self.ui.browseIgBlast.clicked.connect(self.browseIgBlast)
		self.ui.browseIMGT.clicked.connect(self.browseIMGT)
		self.ui.radioButtonChain.clicked.connect(self.ChainClicked)
		self.ui.lineEditRep1.textChanged.connect(self.updateName)
		self.ui.lineEditRep2.textChanged.connect(self.updateName)
		self.ui.lineEditRep3.textChanged.connect(self.updateName)
		self.ui.rdoChoose.clicked.connect(self.updateGroupSetting)
		self.ui.rdoFunction.clicked.connect(self.updateGroupSetting)
		self.ui.checkBoxFileStruc.clicked.connect(self.updateGroupSetting)
		self.ui.toolButtonIgFasta.clicked.connect(self.browseIgBlastFasta)
		#self.ui.radioButtonAllcontig.clicked.connect(self.update10x)
		#self.ui.radioButtonConsensus.clicked.connect(self.update10x)
		#self.ui.radioButtonFiltercontig.clicked.connect(self.update10x)

		self.path10x = ''
		self.pathFasta = ''
		self.pathIMGT = ''
		self.pathIgBlast = ''
		self.pathCSV = ''
		self.pathVDB = ''

		global answer3
		answer3 = 'No'

	def update10x(self, directory):
		if isinstance(directory, str):
			pass
		else:
			if self.path10x == '':
				return
			else:
				pass

		choosed_seq = 'consensus.fasta'
		choose_annotate = 'filtered_contig_annotations.csv'

		'''
		choosed_seq = ''
		choose_annotate = ''
		if self.ui.radioButtonFiltercontig.isChecked():
			choosed_seq = 'filtered_contig.fasta'
			choose_annotate = 'filtered_contig_annotations.csv'
		elif self.ui.radioButtonConsensus.isChecked():
			choosed_seq = 'consensus.fasta'
			choose_annotate = 'filtered_contig_annotations.csv'
		elif self.ui.radioButtonAllcontig.isChecked():
			choosed_seq = 'all_contig.fasta'
			choose_annotate = 'all_contig_annotations.csv'
		'''
		if choosed_seq == '':
			return

		if isinstance(directory, str):
			self.path10x = directory

		seq_file = os.path.join(self.path10x, choosed_seq)
		anno_file = os.path.join(self.path10x, choose_annotate)
		if os.path.exists(seq_file):
			if os.path.exists(anno_file):
				self.ui.Seqpath.setText(seq_file)
				self.ui.Annopath.setText(anno_file)
			else:
				self.ui.Seqpath.setText(seq_file)
		else:
			seq_file = os.path.join(self.path10x, 'outs', choosed_seq)
			anno_file = os.path.join(self.path10x, 'outs', choose_annotate)
			if os.path.exists(seq_file):
				if os.path.exists(anno_file):
					self.ui.Seqpath.setText(seq_file)
					self.ui.Annopath.setText(anno_file)
				else:
					self.ui.Seqpath.setText(seq_file)
			else:
				Msg = 'Can not find consensus.fasta under your folder! Please check your input!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

	def updateGroupSetting(self):
		if self.ui.rdoChoose.isChecked():
			self.ui.comboBoxProject.clear()
			self.ui.comboBoxGroup.clear()
			self.ui.comboBoxSubgroup.clear()

			fields = ['Project']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			if len(DataIn) > 0:
				for item in DataIn:
					self.ui.comboBoxProject.addItem(item[0])
			DataIn.clear()

			fields = ['Grouping']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if len(DataIn) > 0:
				for item in DataIn:
					self.ui.comboBoxGroup.addItem(item[0])
			DataIn.clear()

			fields = ['SubGroup']  # , 'Grouping', 'SubGroup'
			SQLStatement1 = Vgenes.MakeSQLStatement(fields)

			SQLStatement = SQLStatement1[:7] + 'DISTINCT ' + SQLStatement1[7:]

			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

			if len(DataIn) > 0:
				for item in DataIn:
					self.ui.comboBoxSubgroup.addItem(item[0])

			self.ui.comboBoxProject.addItem('')
			self.ui.comboBoxGroup.addItem('')
			self.ui.comboBoxSubgroup.addItem('')

			self.ui.comboBoxProject.setCurrentText('')
			self.ui.comboBoxGroup.setCurrentText('')
			self.ui.comboBoxSubgroup.setCurrentText('')

			self.ui.comboBoxProject.setEditable(True)
			self.ui.comboBoxGroup.setEditable(True)
			self.ui.comboBoxSubgroup.setEditable(True)

			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)

		if self.ui.rdoFunction.isChecked():
			self.ui.comboBoxProject.clear()
			self.ui.comboBoxGroup.clear()
			self.ui.comboBoxSubgroup.clear()

			self.ui.comboBoxProject.setEnabled(False)
			self.ui.comboBoxGroup.setEnabled(False)
			self.ui.comboBoxSubgroup.setEnabled(False)

		if self.ui.checkBoxFileStruc.isChecked():
			self.ui.comboBoxProject.clear()
			self.ui.comboBoxGroup.clear()
			self.ui.comboBoxSubgroup.clear()

			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)

			if self.ui.tabWidget.currentIndex() == 0:
				if self.path10x == '':
					return
				else:
					dirs = self.path10x.split('/')
					if len(dirs) > 3:
						dirs = dirs[-3:]

					if len(dirs) == 3:
						project = dirs[0]
						group = dirs[1]
						subgroup = dirs[2]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
						self.ui.comboBoxSubgroup.addItem(subgroup)
						self.ui.comboBoxSubgroup.setCurrentText(subgroup)
					elif len(dirs) == 2:
						project = dirs[0]
						group = dirs[1]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
					elif len(dirs) == 1:
						project = dirs[0]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
					else:
						pass
			elif self.ui.tabWidget.currentIndex() == 1:
				if self.pathFasta == '':
					return
				else:
					dirs = self.pathFasta[0].split('/')
					if len(dirs) > 3:
						dirs = dirs[-3:]

					if len(dirs) == 3:
						project = dirs[0]
						group = dirs[1]
						subgroup = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
						self.ui.comboBoxSubgroup.addItem(subgroup)
						self.ui.comboBoxSubgroup.setCurrentText(subgroup)
					elif len(dirs) == 2:
						project = dirs[0]
						group = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
					elif len(dirs) == 1:
						project = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
					else:
						pass
			elif self.ui.tabWidget.currentIndex() == 4:
				if self.pathIgBlast == '':
					return
				else:
					dirs = self.pathIgBlast.split('/')
					if len(dirs) > 3:
						dirs = dirs[-3:]

					if len(dirs) == 3:
						project = dirs[0]
						group = dirs[1]
						subgroup = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
						self.ui.comboBoxSubgroup.addItem(subgroup)
						self.ui.comboBoxSubgroup.setCurrentText(subgroup)
					elif len(dirs) == 2:
						project = dirs[0]
						group = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
					elif len(dirs) == 1:
						project = '$FileName'
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
					else:
						pass
			elif self.ui.tabWidget.currentIndex() == 5:
				if self.pathIMGT == '':
					return
				else:
					dirs = self.pathIMGT.split('/')
					if len(dirs) > 3:
						dirs = dirs[-3:]

					dirs[-1] = re.sub(r'\W.+','',dirs[-1])

					if len(dirs) == 3:
						project = dirs[0]
						group = dirs[1]
						subgroup = dirs[2]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
						self.ui.comboBoxSubgroup.addItem(subgroup)
						self.ui.comboBoxSubgroup.setCurrentText(subgroup)
					elif len(dirs) == 2:
						project = dirs[0]
						group = dirs[1]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
						self.ui.comboBoxGroup.addItem(group)
						self.ui.comboBoxGroup.setCurrentText(group)
					elif len(dirs) == 1:
						project = dirs[0]
						self.ui.comboBoxProject.addItem(project)
						self.ui.comboBoxProject.setCurrentText(project)
					else:
						pass

	def updateName(self):
		ori_name = 'clonotype2_consensus_1'
		rep1 = self.ui.lineEditRep1.text()
		if self.ui.radioButtonChain.isChecked():
			rep2 = 'H'
		else:
			rep2 = self.ui.lineEditRep2.text()
		prefix = self.ui.lineEditRep3.text()

		new_name = reName(ori_name, rep1, rep2, prefix)
		self.ui.labelName.setText(new_name)

	def ChainClicked(self):
		if self.ui.radioButtonChain.isChecked():
			self.ui.lineEditRep2.setEnabled(False)
		else:
			self.ui.lineEditRep2.setEnabled(True)
		self.updateName()

	def browseVDB(self):
		files, filetype = QtWidgets.QFileDialog.getOpenFileNames(self, "getOpenFileNames", "~/Documents",
		                                                          "VGene DB Files (*.vdb);;All Files (*)")
		if len(files) == 0:
			return
		else:
			self.ui.listWidgetVDB.addItems(files)
			self.pathVDB = files

	def browse10x(self):
		directory = QtWidgets.QFileDialog.getExistingDirectory(self, "getExistingDirectory", "～/Documents")
		if directory == None or directory == '':
			return
		else:
			self.update10x(directory)
			self.updateGroupSetting()

	def browseFasta(self):
		files, filetype = QtWidgets.QFileDialog.getOpenFileNames(self, "getOpenFileNames", "~/Documents",
		                                                   "Fasta Files (*.fasta);Fasta Files (*.fas);Fasta Files (*.fa);All Files (*)")
		if len(files) == 0:
			return
		else:
			self.ui.listWidgetFasta.addItems(files)
			self.pathFasta = files
			self.updateGroupSetting()

	def browseCSV(self):
		files, filetype = QtWidgets.QFileDialog.getOpenFileNames(self, "getOpenFileNames", "~/Documents",
		                                                   "VGene exported CSV Files (*.csv);;All Files (*)")
		if len(files) == 0:
			return
		else:
			self.ui.listWidgetCSV.addItems(files)
			self.pathCSV = files

	def browseIgBlast(self):
		file, filetype = QtWidgets.QFileDialog.getOpenFileName(self, "getOpenFileName", "~/Documents",
		                                                   "igBlast output File (*.txt);;All Files (*)")
		if file == '' or file == None:
			return
		else:
			self.ui.lineEditIgOut.setText(file)
			self.pathIgBlast = file

	def browseIgBlastFasta(self):
		file, filetype = QtWidgets.QFileDialog.getOpenFileName(self, "getOpenFileName", "~/Documents",
		                                                       "Fasta Files (*.fasta);Fasta Files (*.fas);Fasta Files (*.fa);All Files (*)")
		if file == '' or file == None:
			return
		else:
			self.ui.lineEditIgFasta.setText(file)

	def browseIMGT(self):
		file, filetype = QtWidgets.QFileDialog.getOpenFileName(self, "getOpenFileName", "~/Documents",
		                                                   "IMGT output Zip File (*.zip);;All Files (*)")
		if len(file) == 0:
			return
		else:
			self.ui.lineEditIMGT.setText(file)
			self.pathIMGT = file

	def switchTab(self, num):
		self.updateGroupSetting()
		if num == 0:
			self.ui.radioButtonHuman.setEnabled(True)
			self.ui.radioButtonMouse.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.rdoVandJ.setEnabled(True)
			self.ui.rdoFunction.setEnabled(True)
			self.ui.rdoAll.setEnabled(True)
			self.ui.rdoChoose.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.checkBoxFileStruc.setEnabled(True)
			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)
			self.ui.txtComment.setEnabled(True)
		elif num == 1:
			self.ui.radioButtonHuman.setEnabled(True)
			self.ui.radioButtonMouse.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.rdoVandJ.setEnabled(True)
			self.ui.rdoFunction.setEnabled(True)
			self.ui.rdoAll.setEnabled(True)
			self.ui.rdoChoose.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.checkBoxFileStruc.setEnabled(True)
			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)
			self.ui.txtComment.setEnabled(True)
		elif num == 2:
			self.ui.radioButtonHuman.setEnabled(False)
			self.ui.radioButtonMouse.setEnabled(False)
			self.ui.rdoProductive.setEnabled(False)
			self.ui.rdoVandJ.setEnabled(False)
			self.ui.rdoFunction.setEnabled(False)
			self.ui.rdoAll.setEnabled(False)
			self.ui.rdoChoose.setEnabled(False)
			self.ui.rdoProductive.setEnabled(False)
			self.ui.checkBoxFileStruc.setEnabled(False)
			self.ui.comboBoxProject.setEnabled(False)
			self.ui.comboBoxGroup.setEnabled(False)
			self.ui.comboBoxSubgroup.setEnabled(False)
			self.ui.txtComment.setEnabled(False)
		elif num == 3:
			self.ui.radioButtonHuman.setEnabled(False)
			self.ui.radioButtonMouse.setEnabled(False)
			self.ui.rdoProductive.setEnabled(False)
			self.ui.rdoVandJ.setEnabled(False)
			self.ui.rdoFunction.setEnabled(False)
			self.ui.rdoAll.setEnabled(False)
			self.ui.rdoChoose.setEnabled(False)
			self.ui.rdoProductive.setEnabled(False)
			self.ui.checkBoxFileStruc.setEnabled(False)
			self.ui.comboBoxProject.setEnabled(False)
			self.ui.comboBoxGroup.setEnabled(False)
			self.ui.comboBoxSubgroup.setEnabled(False)
			self.ui.txtComment.setEnabled(False)
		elif num == 4:
			self.ui.radioButtonHuman.setEnabled(True)
			self.ui.radioButtonMouse.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.rdoVandJ.setEnabled(True)
			self.ui.rdoFunction.setEnabled(True)
			self.ui.rdoAll.setEnabled(True)
			self.ui.rdoChoose.setEnabled(True)
			self.ui.rdoProductive.setEnabled(True)
			self.ui.checkBoxFileStruc.setEnabled(True)
			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)
			self.ui.txtComment.setEnabled(True)
		elif num == 5:
			self.ui.radioButtonHuman.setEnabled(False)
			self.ui.radioButtonMouse.setEnabled(False)
			self.ui.rdoVandJ.setEnabled(False)
			self.ui.rdoFunction.setEnabled(True)
			self.ui.rdoAll.setEnabled(False)
			self.ui.rdoChoose.setEnabled(True)
			self.ui.rdoProductive.setEnabled(False)
			self.ui.checkBoxFileStruc.setEnabled(True)
			self.ui.comboBoxProject.setEnabled(True)
			self.ui.comboBoxGroup.setEnabled(True)
			self.ui.comboBoxSubgroup.setEnabled(True)
			self.ui.txtComment.setEnabled(True)
		else:
			pass

	def accept(self):
		num = self.ui.tabWidget.currentIndex()

		self.ui.progressBar.setValue(1)
		self.ui.labelpct.setText('Start processing and loading...')

		if num == 0:
			self.InitiateImportFrom10X('none', 0)
		elif num == 1:
			self.InitiateImportFromFasta('none', 0)
		elif num == 2:
			self.InitiateImportFromCSV('none', 0)
		elif num == 3:
			self.InitiateImportFromVDB('none', 0)
		elif num == 4:
			self.InitiateImportFromIgBlast('none', 0)
		elif num == 5:
			self.InitiateImportFromIMGT()
		else:
			pass

	def reject(self):
		print('close')
		self.hide()

	def disableWidgets(self):
		self.ui.tabWidget.setEnabled(False)

		self.ui.comboBoxGroup.setDisabled(True)
		self.ui.comboBoxProject.setDisabled(True)
		self.ui.comboBoxSubgroup.setDisabled(True)

		self.ui.radioButtonHuman.setDisabled(True)
		self.ui.radioButtonMouse.setDisabled(True)

		self.ui.rdoAll.setDisabled(True)
		self.ui.rdoProductive.setDisabled(True)
		self.ui.rdoVandJ.setDisabled(True)

		self.ui.rdoChoose.setDisabled(True)
		self.ui.rdoFunction.setDisabled(True)
		self.ui.checkBoxFileStruc.setDisabled(True)

		self.ui.txtComment.setDisabled(True)
		self.ui.pushButtonCancel.setEnabled(False)
		self.ui.pushButtonOK.setEnabled(False)

	def checkProgress(self):
		global timer
		progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
		file_handle = open(progressBarFile, 'r')
		text = file_handle.readline()
		text_list = text.split(',')

		try:
			progress = text_list[0]
			progress_text = text_list[1]
			progress = int(float(progress))
			self.ui.progressBar.setValue(progress)
			self.ui.labelpct.setText(progress_text + '(' + str(progress) + '%) records loaded')
		except:
			progress = self.ui.progressBar.value()
		if progress > 99:
			self.ui.labelpct.setText('Loading finished!')
			self.ui.progressBar.setValue(100)
			return

		if StopCheckProgress:
			return
		else:
			t = thd.Timer(1, self.checkProgress)
			t.start()

	def readBarcode(self, anno_path_name):
		# read annotation content
		barcode_dict = {}
		if anno_path_name != "":
			csvFile = open(anno_path_name, "r")
			reader = csv.reader(csvFile)
			for item in reader:
				# ignore header line
				if reader.line_num == 1:
					continue
				barcode_dict[item[17]] = item[0]
			csvFile.close()
		return barcode_dict

	def InitiateImportFrom10X(self, Filenamed, MaxNum):
		global StopCheckProgress
		self.calling = 1

		# need to transfer species grouping to IgBlaster
		answer = ''
		thetype = 'FASTA'
		species = ''
		datalist = []
		global answer3
		# answerTo = answer3

		if self.ui.radioButtonHuman.isChecked():
			species = 'Human'
		elif self.ui.radioButtonMouse.isChecked():
			species = 'Mouse'

		seq_pathname = self.ui.Seqpath.text()
		anno_path_name = self.ui.Annopath.text()
		self.anno_path_name = anno_path_name

		self.rep1 = self.ui.lineEditRep1.text()
		if self.ui.radioButtonChain.isChecked():
			self.rep2 = 'byChain'
		else:
			self.rep2 = self.ui.lineEditRep2.text()
		self.prefix = self.ui.lineEditRep3.text()

		if self.ui.rdoProductive.isChecked() == True:
			GetProductive = 0
		elif self.ui.rdoVandJ.isChecked() == True:
			GetProductive = 1
		else:
			GetProductive = 2

		if seq_pathname == None:
			return
		answer2 = ''

		ErlogFile = os.path.join(temp_folder,'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		header = "Began input at " + time.strftime('%c')
		with open(ErlogFile2, 'w') as currentFile:
			currentFile.write(header)
		# firstOne = True

		if self.ui.rdoChoose.isChecked() or self.ui.checkBoxFileStruc.isChecked():
			if Filenamed == 'none':
				project = self.ui.comboBoxProject.currentText()
				grouping = self.ui.comboBoxGroup.currentText()
				subgroup = self.ui.comboBoxSubgroup.currentText()
			else:
				project = Filenamed[1]
				grouping = Filenamed[2]
				subgroup = Filenamed[3]

			if project == '': project = 'none'
			if grouping == '': grouping = 'none'
			if subgroup == '': subgroup = 'none'

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread(self)
			workThread.item = seq_pathname
			workThread.datalist = datalist
			if self.ui.radioButtonFast1.isChecked():
				workThread.method = 'fast'
			else:
				workThread.method = 'slow'
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(seq_pathname)
			f.close()

			self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return
		elif self.ui.rdoFunction.isChecked():
			project = 'ByFunction'
			grouping = ''
			subgroup = ''

			multiProject = ''

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)
			datalist.append(multiProject)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread(self)
			workThread.item = seq_pathname
			workThread.datalist = datalist
			if self.ui.radioButtonFast1.isChecked():
				workThread.method = 'fast'
			else:
				workThread.method = 'slow'
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(seq_pathname)
			f.close()

			#self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return

	def InitiateImportFromFasta(self, Filenamed, MaxNum):
		self.calling = 2
		global StopCheckProgress

		# need to transfer species grouping to IgBlaster
		if self.pathFasta == "":
			return
		answer = ''
		thetype = 'FASTA'
		species = ''
		datalist = []
		global answer3
		# answerTo = answer3

		if self.ui.radioButtonHuman.isChecked():
			species = 'Human'
		elif self.ui.radioButtonMouse.isChecked():
			species = 'Mouse'

		# process sequence
		time_stamp = str(int(time.time() * 100)) + '.fasta'
		seq_pathname = os.path.join(temp_folder,time_stamp)
		fout = open(seq_pathname,'w')
		fasta_files = self.pathFasta
		for fasta_file in fasta_files:
			file_name = fasta_file.split('/')[-1]
			file_name = re.sub(r'\..+','',file_name)
			try:
				for record in SeqIO.parse(fasta_file, "fasta"):
					#print(record)
					seq_name = '>' + file_name + '_' + record.id
					fout.write(seq_name + '\n')
					fout.write(record.seq._data + '\n')
			except:
				Msg = 'Can not parse file\n ' + fasta_file + '\n as fasta file! Check your input!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return
		fout.close()

		# settings
		if self.ui.rdoProductive.isChecked() == True:
			GetProductive = 0
		elif self.ui.rdoVandJ.isChecked() == True:
			GetProductive = 1
		else:
			GetProductive = 2

		if seq_pathname == None:
			return
		answer2 = ''

		ErlogFile = os.path.join(temp_folder,'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		header = "Began input at " + time.strftime('%c')
		with open(ErlogFile2, 'w') as currentFile:
			currentFile.write(header)
		# firstOne = True

		if self.ui.checkBoxFileStruc.isChecked():
			project = self.ui.comboBoxProject.currentText()
			grouping = self.ui.comboBoxGroup.currentText()
			subgroup = self.ui.comboBoxSubgroup.currentText()

			datalist.clear()
			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread(self)
			workThread.item = seq_pathname
			workThread.datalist = datalist
			if self.ui.radioButtonFast2.isChecked():
				workThread.method = 'fast'
			else:
				workThread.method = 'slow'
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			file_list_str = '\n'.join(fasta_files)
			f.write(file_list_str)
			f.close()

			self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return
		elif self.ui.rdoChoose.isChecked():
			(dirname, filename) = os.path.split(seq_pathname)

			if Filenamed == 'none':
				project = self.ui.comboBoxProject.currentText()
				grouping = self.ui.comboBoxGroup.currentText()
				subgroup = self.ui.comboBoxSubgroup.currentText()
			else:
				project = Filenamed[1]
				grouping = Filenamed[2]
				subgroup = Filenamed[3]

			if project == '': project = 'none'
			if grouping == '': grouping = 'none'
			if subgroup == '': subgroup = 'none'

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread(self)
			workThread.item = seq_pathname
			workThread.datalist = datalist
			if self.ui.radioButtonFast2.isChecked():
				workThread.method = 'fast'
			else:
				workThread.method = 'slow'
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			file_list_str = '\n'.join(fasta_files)
			f.write(file_list_str)
			f.close()

			self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return
		elif self.ui.rdoFunction.isChecked():
			project = 'ByFunction'
			grouping = ''
			subgroup = ''

			multiProject = ''

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)
			datalist.append(multiProject)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread(self)
			workThread.item = seq_pathname
			workThread.datalist = datalist
			if self.ui.radioButtonFast2.isChecked():
				workThread.method = 'fast'
			else:
				workThread.method = 'slow'
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			file_list_str = '\n'.join(fasta_files)
			f.write(file_list_str)
			f.close()

			#self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return

	def InitiateImportFromCSV(self, Filenamed, MaxNu):
		global FieldList
		global FieldCommentList
		global FieldTypeList
		global RealNameList

		self.calling = 3
		files = self.pathCSV

		if os.path.isfile(DBFilename):
			# collect information, alter the DB structure
			Msg = 'We added new fileds from the follow VDBs:\n'
			for csv_file in files:
				csvFile = open(csv_file, "r")
				reader = csv.reader(csvFile)
				line = 0
				for item in reader:
					if line == 0:
						cur_field = item
					elif line == 1:
						cur_field_name = item
					else:
						break
					line += 1
				csvFile.close()

				new_cols = [i for i in cur_field if i not in FieldList]

				if len(new_cols) > 0:
					for new_col in new_cols:
						new_col_index = cur_field.index(new_col)
						new_col_name = cur_field_name[new_col_index]
						new_col_comment = ''

						# update table
						SQLSTATEMENT1 = "ALTER TABLE vgenesDB ADD " + new_col + " text"
						SQLSTATEMENT2 = 'INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment) ' \
						                'VALUES(' + str(len(FieldList) + 1) + ',"' + new_col + '", "' + new_col_name + \
						                '", "Customized", "' + new_col_comment + '")'
						try:
							VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
						except:
							msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
							QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
							                    QMessageBox.Ok)
							return

						try:
							VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
						except:
							msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
							QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
							                    QMessageBox.Ok)
							return

						RealNameList.append(new_col)
						FieldCommentList.append(new_col_comment)
						FieldTypeList.append('Customized')
						FieldList.append(new_col_name)

						Msg = Msg + new_col + ' from    ' + csv_file + '\n'

			QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)

			# import data from VDBs
			print('import data from VDBs')
			for csv_file in files:

				csvFile = open(csv_file, "r")
				reader = csv.reader(csvFile)
				line = 0
				DataIn = []
				for item in reader:
					if line == 0:
						cur_field = item
					elif line == 1:
						cur_field_name = item
					else:
						DataIn.append(item)
					line += 1
				csvFile.close()

				question_list = ['?' for n in range(len(cur_field))]
				field_str = ','.join(cur_field)
				question_str = ','.join(question_list)

				# update ID (ID should be unique)
				SQLSTATEMENT = 'SELECT COUNT(ID) FROM vgenesDB'
				Count = VGenesSQL.RunSQL(DBFilename, SQLSTATEMENT)
				Count = Count[0][0]

				InputData = []
				for records in DataIn:
					cur_line = list(records)
					cur_line[119] = str(int(records[119]) + Count)
					cur_line[58] = re.sub('#','\n',records[58])
					cur_line[97] = re.sub('|', ',', records[97])
					InputData.append(cur_line)

				# insert into DB
				print('insert into DB')
				conn = db.connect(DBFilename)
				cursor = conn.cursor()
				SQLSTATEMENT = "INSERT INTO vgenesDB(" + field_str + ") VALUES(" + question_str + ")"
				cursor.executemany(SQLSTATEMENT, InputData)
				conn.commit()
				conn.close()
				print(csv_file)
			Msg = 'Data import finished!'
			QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)

			Vgenes.LoadDB(DBFilename)
			self.hide()

	def InitiateImportFromVDB(self, Filenamed, MaxNu):
		global FieldList
		global FieldCommentList
		global FieldTypeList
		global RealNameList

		self.calling = 4
		files = self.pathVDB

		if os.path.isfile(DBFilename):
			# collect information, alter the DB structure
			Msg = 'We added new fileds from the follow VDBs:\n'
			for vdb_file in files:
				VGenesSQL.checkFieldTable(vdb_file)

				SQLSTATEMENT = 'SELECT Field,FieldNickName,FieldComment,FieldType FROM fieldsname'
				DataIn = VGenesSQL.RunSQL(vdb_file, SQLSTATEMENT)
				cur_field = [i[0] for i in DataIn]
				cur_field_name = [i[1] for i in DataIn]
				cur_field_comment = [i[2] for i in DataIn]
				cur_field_type = [i[3] for i in DataIn]

				new_cols = [i for i in cur_field if i not in FieldList]

				if len(new_cols) > 0:
					for	new_col in new_cols:
						new_col_index = cur_field.index(new_col)
						new_col_name = cur_field_name[new_col_index]
						new_col_comment = cur_field_comment[new_col_index]

						# update table
						SQLSTATEMENT1 = "ALTER TABLE vgenesDB ADD " + new_col + " text"
						SQLSTATEMENT2 = 'INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment) ' \
						                'VALUES(' + str(len(FieldList) + 1) + ',"' + new_col + '", "' + new_col_name + \
						                '", "Customized", "' + new_col_comment + '")'
						try:
							VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT1)
						except:
							msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT1
							QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
							                    QMessageBox.Ok)
							return

						try:
							VGenesSQL.RunUpdateSQL(DBFilename, SQLSTATEMENT2)
						except:
							msg = "DB operation Error! Current SQL statement is: \n" + SQLSTATEMENT2
							QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok,
							                    QMessageBox.Ok)
							return

						RealNameList.append(new_col)
						FieldCommentList.append(new_col_comment)
						FieldTypeList.append('Customized')
						FieldList.append(new_col_name)

						Msg = Msg + new_col + ' from    ' + vdb_file + '\n'

			QMessageBox.information(self, 'Information', Msg,QMessageBox.Ok, QMessageBox.Ok)

			# import data from VDBs
			print('import data from VDBs')
			for vdb_file in files:
				SQLSTATEMENT = 'SELECT Field FROM fieldsname'
				DataIn = VGenesSQL.RunSQL(vdb_file, SQLSTATEMENT)
				field_str = ''
				question_str = ''
				for item in DataIn:
					field_str = field_str + item[0] + ','
					question_str = question_str + '?,'
				field_str = field_str.rstrip(',')
				question_str = question_str.rstrip(',')

				SQLSTATEMENT = 'SELECT ' + field_str + ' FROM vgenesDB'
				DataIn = VGenesSQL.RunSQL(vdb_file, SQLSTATEMENT)

				# update ID (ID should be unique)
				SQLSTATEMENT = 'SELECT COUNT(ID) FROM vgenesDB'
				Count = VGenesSQL.RunSQL(DBFilename, SQLSTATEMENT)
				Count = Count[0][0]

				InputData = []
				for records in DataIn:
					cur_line = list(records)
					cur_line[119] = str(int(records[119]) + Count)
					InputData.append(cur_line)

				# insert into DB
				print('insert into DB')
				conn = db.connect(DBFilename)
				cursor = conn.cursor()
				SQLSTATEMENT = "INSERT INTO vgenesDB(" + field_str + ") VALUES(" + question_str + ")"
				cursor.executemany(SQLSTATEMENT, InputData)
				conn.commit()
				conn.close()
				print(vdb_file)
			Msg = 'Data import finished!'
			QMessageBox.information(self, 'Information', Msg,QMessageBox.Ok, QMessageBox.Ok)

			Vgenes.LoadDB(DBFilename)
			self.hide()
		else:
			return

	def InitiateImportFromIgBlast(self, Filenamed, MaxNum):
		if os.path.isfile(self.ui.lineEditIgOut.text()) and os.path.isfile(self.ui.lineEditIgFasta.text()):
			pass
		else:
			msg = 'Please setup IgOut file and Fasta file!'
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		self.calling = 4

		# need to transfer species grouping to IgBlaster
		answer = ''
		thetype = 'FASTA'
		species = ''
		datalist = []
		global answer3
		# answerTo = answer3

		if self.ui.radioButtonHuman.isChecked():
			species = 'Human'
		elif self.ui.radioButtonMouse.isChecked():
			species = 'Mouse'

		igOut = self.ui.lineEditIgOut.text()
		seq_pathname = self.ui.lineEditIgFasta.text()

		if self.ui.rdoProductive.isChecked() == True:
			GetProductive = 0
		elif self.ui.rdoVandJ.isChecked() == True:
			GetProductive = 1
		else:
			GetProductive = 2

		if seq_pathname == None:
			return
		if igOut == None:
			return
		answer2 = ''

		ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		header = "Began input at " + time.strftime('%c')
		with open(ErlogFile2, 'w') as currentFile:
			currentFile.write(header)
		# firstOne = True

		if self.ui.rdoChoose.isChecked() or self.ui.checkBoxFileStruc.isChecked():
			if Filenamed == 'none':
				project = self.ui.comboBoxProject.currentText()
				grouping = self.ui.comboBoxGroup.currentText()
				subgroup = self.ui.comboBoxSubgroup.currentText()
			else:
				project = Filenamed[1]
				grouping = Filenamed[2]
				subgroup = Filenamed[3]

			if project == '': project = 'none'
			if grouping == '': grouping = 'none'
			if subgroup == '': subgroup = 'none'

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread1(self)
			workThread.igOut = igOut
			workThread.item = seq_pathname
			workThread.datalist = datalist
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(seq_pathname)
			f.close()

			self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return
		elif self.ui.rdoFunction.isChecked():
			project = 'ByFunction'
			grouping = ''
			subgroup = ''

			multiProject = ''

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)
			datalist.append(species)
			datalist.append(GetProductive)
			datalist.append(MaxNum)
			datalist.append(multiProject)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThread1(self)
			workThread.igOut = igOut
			workThread.item = seq_pathname
			workThread.datalist = datalist
			workThread.start()
			workThread.trigger.connect(self.multi_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(seq_pathname)
			f.close()

			# self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return

	def InitiateImportFromIMGT(self):
		if os.path.isfile(self.ui.lineEditIMGT.text()):
			pass
		else:
			msg = 'Please setup IGMT output file!'
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		self.calling = 5

		# need to transfer species grouping to IgBlaster
		datalist = []
		global answer3
		# answerTo = answer3

		IMGT_out = self.ui.lineEditIMGT.text()

		answer2 = ''

		ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		header = "Began input at " + time.strftime('%c')
		with open(ErlogFile2, 'w') as currentFile:
			currentFile.write(header)
		# firstOne = True

		if self.ui.rdoChoose.isChecked() or self.ui.checkBoxFileStruc.isChecked():
			project = self.ui.comboBoxProject.currentText()
			grouping = self.ui.comboBoxGroup.currentText()
			subgroup = self.ui.comboBoxSubgroup.currentText()

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThreadIMGTparser(self)
			workThread.item = IMGT_out
			workThread.datalist = datalist
			workThread.start()
			workThread.trigger.connect(self.multiIMGT_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(IMGT_out)
			f.close()

			self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return
		elif self.ui.rdoFunction.isChecked():
			project = 'ByFunction'
			grouping = ''
			subgroup = ''

			multiProject = ''

			datalist.clear()

			datalist.append(project)
			datalist.append(grouping)
			datalist.append(subgroup)

			# try multi-thread
			progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
			file_handle = open(progressBarFile, 'w')
			file_handle.write('0')
			file_handle.close()
			workThread = WorkThreadIMGTparser(self)
			workThread.item = IMGT_out
			workThread.datalist = datalist
			workThread.start()
			workThread.trigger.connect(self.multiIMGT_callback)

			import_file = os.path.join(temp_folder, "import_file_name.txt")
			f = open(import_file, 'w')
			f.write(IMGT_out)
			f.close()

			# self.disableWidgets()
			StopCheckProgress = False
			self.checkProgress()
			return

	@pyqtSlot()
	def multi_callback(self):
		global IgBLASTAnalysis
		global StopCheckProgress

		StopCheckProgress = True

		Startprocessed = 0
		try:
			Startprocessed = len(IgBLASTAnalysis)
			self.close()
		except:
			if Startprocessed == 0:
				self.close()

		a = IgBLASTAnalysis

		if self.calling == 1:
			# match barcode
			barcodeDict = self.readBarcode(self.anno_path_name)

			for record in IgBLASTAnalysis:
				if record[0] in barcodeDict.keys():
					record[108] = barcodeDict[record[0]]

				if self.rep2 == "byChain":
					rep2 = record[2][0]
				else:
					rep2 = self.rep2
				record[0] = reName(record[0], self.rep1, rep2, self.prefix)
		elif self.calling == 2:
			for record in IgBLASTAnalysis:
				sampleName = record[0].split('_')[0]
				if len(record[75]) > 0:
					if record[75][0] == '$':
						record[75] = sampleName
						continue
				if len(record[76]) > 0:
					if record[76][0] == '$':
						record[76] = sampleName
						continue
				if len(record[77]) > 0:
					if record[77][0] == '$':
						record[77] = sampleName
						continue
		#a = IgBLASTAnalysis

		Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)

		import_file = os.path.join(temp_folder, "import_file_name.txt")
		file_handle = open(import_file, 'r')
		file_name = file_handle.readlines()
		file_name = ''.join(file_name)
		file_handle.close()

		i = 0
		newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: \n' + file_name + '\n'

		ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		with open(ErlogFile, 'r') as currentFile:  # using with for this automatically closes the file even if you crash
			for line in currentFile:
				if i > 0:
					newErLog += line
				i += 1

		with open(ErlogFile2, 'a') as currentFile:
			currentFile.write(newErLog)

		Vgenes.LoadDB(DBFilename)
		self.ShowVGenesText(ErlogFile2)
		self.hide()

	@pyqtSlot()
	def multiIMGT_callback(self):
		global IMGTAnalysis
		global StopCheckProgress

		StopCheckProgress = True

		Startprocessed = 0
		try:
			Startprocessed = len(IMGTAnalysis)
			self.close()
		except:
			if Startprocessed == 0:
				self.close()

		a = IMGTAnalysis

		Processed, answer = VGenesSQL.enterData(self, DBFilename, IMGTAnalysis, answer3)

		import_file = os.path.join(temp_folder, "import_file_name.txt")
		file_handle = open(import_file, 'r')
		file_name = file_handle.readlines()
		file_name = ''.join(file_name)
		file_handle.close()

		i = 0
		newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: \n' + file_name + '\n'

		ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
		with open(ErlogFile, 'r') as currentFile:  # using with for this automatically closes the file even if you crash
			for line in currentFile:
				if i > 0:
					newErLog += line
				i += 1

		with open(ErlogFile2, 'a') as currentFile:
			currentFile.write(newErLog)

		Vgenes.LoadDB(DBFilename)
		self.ShowVGenesText(ErlogFile2)
		self.hide()

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

			now = 'FASTA' + time.strftime('%c') + '.nt'
			FASTAFileName = os.path.join(temp_folder, now)
			# need to test

			with open(FASTAFileName,
			          'w') as currentFile:  # using with for this automatically closes the file even if you crash
				currentFile.write(FinalFASTA)

			return FASTAFileName

	def ShowVGenesText(self, filename):

		self.TextEdit.show()
		if filename != '':
			self.TextEdit.loadFile(filename)

class ImportDialogue(QtWidgets.QDialog, Ui_DialogImport):
	def __init__(self, parent=None):
		QtWidgets.QDialog.__init__(self, parent)
		self.setupUi(self)
		self.TextEdit = VGenesTextMain()

		self.toolButton.clicked.connect(self.browsedir)

		global answer3
		answer3 = 'No'

	def browsedir(self):  # browse and select path
		out_dir = openFiles(self, thetype)
		if out_dir == '' or out_dir == None:
			return
		self.lineEdit.setText(out_dir)

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

	def disableWidgets(self):
		self.comboBoxGroup.setDisabled(True)
		self.comboBoxProject.setDisabled(True)
		self.comboBoxSubgroup.setDisabled(True)

		self.radioButFASTA.setDisabled(True)
		self.radioButHuman.setDisabled(True)
		self.radioButIndSeq.setDisabled(True)
		self.radioButMouse.setDisabled(True)

		self.rdoAll.setDisabled(True)
		self.rdoProductive.setDisabled(True)
		self.rdoVandJ.setDisabled(True)

		self.rdoChoose.setDisabled(True)
		self.rdoFunction.setDisabled(True)
		self.checkBoxFileStruc.setDisabled(True)

		self.txtComment.setDisabled(True)
		self.MaxImport.setDisabled(True)

		#self.btnImportOldVGenes.setDisabled(True)
		self.buttonBox.setDisabled(True)

	def checkProgress(self):
		global timer
		progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
		file_handle = open(progressBarFile,'r')
		text = file_handle.readline()
		text_list = text.split(',')

		try:
			progress = text_list[0]
			progress_text = text_list[1]
			progress = int(float(progress))
			self.progressBar.setValue(progress)
			self.labelpct.setText(progress_text + '(' + str(progress) + '%) records loaded')
		except:
			progress = self.progressBar.value()
		if progress > 99:
			self.labelpct.setText('Loading finished!')
			self.progressBar.setValue(100)
			return
		t=thd.Timer(1, self.checkProgress)
		t.start()

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
			if filenames == None:
				return False
			pathname1 = self.ProcessSeqFiles(filenames)
			pathname = []
			pathname.append(pathname1)
		# filename = filenames[0]

		if self.ui.rdoProductive.isChecked() == True:
			GetProductive = 0
		elif self.ui.rdoVandJ.isChecked() == True:
			GetProductive = 1
		else:
			GetProductive = 2

		if pathname == None:
			return
		answer2 = ''

		ErlogFile = os.path.join(temp_folder, 'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs

		ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
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

				# try multi-thread
				progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
				file_handle = open(progressBarFile, 'w')
				file_handle.write('0')
				file_handle.close()
				workThread = WorkThread(self)
				workThread.item = item
				workThread.datalist = datalist
				workThread.start()
				workThread.trigger.connect(self.multi_callback)

				import_file = os.path.join(temp_folder, "import_file_name.txt")
				f = open(import_file, 'w')
				f.write(item)
				f.close()

				self.disableWidgets()
				self.checkProgress()
				return

				'''
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
				'''
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

				# try multi-thread
				progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
				file_handle = open(progressBarFile, 'w')
				file_handle.write('0')
				file_handle.close()
				workThread = WorkThread(self)
				workThread.item = item
				workThread.datalist = datalist
				workThread.start()
				workThread.trigger.connect(self.multi_callback)

				import_file = os.path.join(temp_folder, "import_file_name.txt")
				f = open(import_file, 'w')
				f.write(item)
				f.close()

				self.disableWidgets()
				self.checkProgress()
				return

				'''
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
				'''
		elif self.rdoFunction.isChecked():
			for item in pathname:
				self.lineEdit.setText(item)
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

				# try multi-thread
				progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
				file_handle = open(progressBarFile, 'w')
				file_handle.write('0')
				file_handle.close()
				workThread = WorkThread(self)
				workThread.item = item
				workThread.datalist = datalist
				workThread.start()
				workThread.trigger.connect(self.multi_callback)

				import_file = os.path.join(temp_folder, "import_file_name.txt")
				f = open(import_file, 'w')
				f.write(item)
				f.close()

				self.disableWidgets()
				self.checkProgress()
				return

				'''				
				start = time.time()
				IgBLASTAnalysis = IgBLASTer.IgBLASTit(item, datalist)
				end = time.time()
				print('Run time for IgBlast: ' +str(end - start))
				Startprocessed = 0
				try:

					Startprocessed = len(IgBLASTAnalysis)
				except:
					if Startprocessed == 0:
						self.close()

				start = time.time()
				Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)
				end = time.time()
				print('Run time for Importing DB: ' + str(end - start))

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
				'''
		Vgenes.LoadDB(DBFilename)
		self.ShowVGenesText(ErlogFile2)

	@pyqtSlot()
	def multi_callback(self):

		Startprocessed = 0
		try:
			Startprocessed = len(IgBLASTAnalysis)
			self.close()
		except:
			if Startprocessed == 0:
				self.close()

		Processed, answer = VGenesSQL.enterData(self, DBFilename, IgBLASTAnalysis, answer3)

		import_file = os.path.join(temp_folder, "import_file_name.txt")
		file_handle = open(import_file,'r')
		file_name = file_handle.readline()
		file_handle.close()

		i = 0
		newErLog = '\n' + str(Processed) + ' sequences were input by IgBLAST for file: ' + file_name + '\n'

		ErlogFile = os.path.join(temp_folder,'ErLog.txt')
		ErlogFile2 = os.path.join(temp_folder,'ErLog2.txt')
		with open(ErlogFile,'r') as currentFile:  # using with for this automatically closes the file even if you crash
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
		from operator import itemgetter  # SeqList.sort(key=itemgetter(0, 1, 2, 3))
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
		try:
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
			FileNamed = os.path.join(working_prefix, 'IgBlast', 'database',
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
		except:
			return

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

			now = 'FASTA' + time.strftime('%c') + '.nt'
			FASTAFileName = os.path.join(temp_folder, now)
			# need to test

			with open(FASTAFileName, 'w') as currentFile:  # using with for this automatically closes the file even if you crash
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
	resizeSignal = pyqtSignal(int, int)
	def __init__(self, parent=None):
		super(ResizeWidget, self).__init__()
		self.id = 0
		self.h = 0
		self.w = 0
		self.html = ''
		self._resize_timer = None

	def updateResizeTimer(self, interval=None):
		if self._resize_timer is not None:
			self.killTimer(self._resize_timer)
		if interval is not None:
			self._resize_timer = self.startTimer(interval)
		else:
			self._resize_timer = None

	def resizeEvent(self, event):
		w = event.size().width()
		h = event.size().height()
		self.h = h
		self.w = w
		self.updateResizeTimer(300)

	def timerEvent(self, event):
		if event.timerId() == self._resize_timer:
			self.updateResizeTimer()

			print(f' size now :{self.w, self.h, self.id}')
			self.resizeSignal.emit(self.w, self.h)

class MyObjectCls(QObject):
	downloadFigSignal = pyqtSignal(str)
	updateSelectionSignal = pyqtSignal(str)

	def __init__(self, parent=None):
		QObject.__init__(self, parent)

	@pyqtSlot(str)
	def consolePrint(self, msg):
		print(msg)

	@pyqtSlot(str)
	def updateSelection(self, msg):
		self.updateSelectionSignal.emit(msg)

	@pyqtSlot(str)
	def download(self, msg):
		self.downloadFigSignal.emit(msg)

class WorkThread(QThread):
	trigger = pyqtSignal(str)

	def __int__(self):
		super(WorkThread, self).__init__()
		self.parent = parent
		self.item = ''
		self.datalist = ''
		self.method = 'slow'

	def run(self):
		global IgBLASTAnalysis
		if self.method == 'fast':
			IgBLASTAnalysis = IgBlastParserFast(self.item, self.datalist)
		else:
			IgBLASTAnalysis = IgBLASTer.IgBLASTit(self.item, self.datalist)

		self.trigger.emit(self.item)

class WorkThread1(QThread):
	trigger = pyqtSignal(str)

	def __int__(self):
		super(WorkThread1, self).__init__()
		self.parent = parent
		self.item = ''
		self.datalist = ''
		self.igOut = ''

	def run(self):
		global IgBLASTAnalysis
		IgBLASTAnalysis = IgBLASTer.IgBLASTitResults(self.item, self.igOut, self.datalist)
		self.trigger.emit(self.item)

class WorkThreadIMGTparser(QThread):
	trigger = pyqtSignal(str)

	def __int__(self):
		super(WorkThread1, self).__init__()
		self.parent = parent
		self.item = ''
		self.datalist = ''

	def run(self):
		global IMGTAnalysis
		IMGTAnalysis = IMGTparser(self.item, self.datalist)
		a = IMGTAnalysis
		self.trigger.emit(self.item)

class VGenesForm(QtWidgets.QMainWindow):
	def __init__(self):  # , parent=None):
		super(VGenesForm, self).__init__()  # parent)
		global VGenesTextWindows
		global DBFilename

		self.ui = Ui_MainWindow()

		self.ui.setupUi(self)

		self.ui.comboBoxSpecies.currentTextChanged.connect(self.on_comboBoxSpecies_editTextChanged)
		self.ui.cboTreeOp1.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboTreeOp2.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboTreeOp3.currentTextChanged.connect(self.TreeviewOptions)
		self.ui.cboDecorate.currentTextChanged.connect(self.DecoratePeptide)
		self.ui.treeWidget.itemChanged.connect(self.handleChanged)
		self.ui.treeWidget.itemSelectionChanged.connect(self.TreeSelectChanged)
		self.ui.treeWidget.clicked['QModelIndex'].connect(self.treeWidgetClicked)
		self.ui.treeWidget.doubleClicked['QModelIndex'].connect(self.CheckMultiple)
		self.ui.cboReportOptions.currentTextChanged.connect(self.ReportOptions)
		self.ui.cboFindField.currentTextChanged.connect(self.on_cboFindField_currentTextChanged)
		self.ui.tabWidget.currentChanged['int'].connect(self.InitialGraphic)
		self.ui.pushButtonDraw.clicked.connect(self.GenerateFigure)
		self.ui.checkBoxFigLegend.clicked.connect(self.GenerateFigure)
		self.ui.checkBoxStack.clicked.connect(self.GenerateFigure)
		self.ui.radioButtonOriginal.clicked.connect(self.GenerateFigure)
		self.ui.radioButtonLog10.clicked.connect(self.GenerateFigure)
		self.ui.radioButtonLog2.clicked.connect(self.GenerateFigure)
		self.ui.checkBoxHideNull.clicked.connect(self.GenerateFigure)
		self.ui.pushButtonDownload.clicked.connect(self.downloadFig)
		self.ui.radioButtonTreeMap.clicked.connect(self.GenerateFigure)
		self.ui.radioButtonTree.clicked.connect(self.GenerateFigure)
		self.ui.pushButtonNT.clicked.connect(self.makeNTLogo)
		self.ui.pushButtonAA.clicked.connect(self.makeAALogo)
		self.ui.toolButtonIgphyml.clicked.connect(self.loadIgphyml)
		self.ui.toolButtonCloneRaxml.clicked.connect(self.buildCloneTree)
		self.ui.EditLock.clicked.connect(self.ChangeEditMode)
		self.ui.checkBoxAll.stateChanged.connect(self.checkAll)
		self.ui.checkBoxAll1.stateChanged.connect(self.checkAll1)
		self.ui.checkBoxRowSelection.stateChanged.connect(self.selectionMode)
		self.ui.pushButtonRefresh.clicked.connect(self.load_table)
		self.ui.SeqTable.clicked.connect(self.table_to_tree_selection)
		self.ui.pushButtonAlter.clicked.connect(self.openAlter)
		self.ui.listWidgetClone.itemSelectionChanged.connect(self.selectClone)
		self.ui.pushButtonDrawClone.clicked.connect(self.GenerateFigureClone)
		self.ui.checkBoxFigLegendClone.clicked.connect(self.GenerateFigureClone)
		self.ui.checkBoxStackClone.clicked.connect(self.GenerateFigureClone)
		self.ui.pushButtonDownloadClone.clicked.connect(self.downloadFigClone)
		self.ui.pushButtonCheckCloone.clicked.connect(self.checkClone)
		self.ui.listWidgetCloneMember.itemSelectionChanged.connect(self.updateClone)
		self.ui.pushButtonNewickTree.clicked.connect(self.loadNewickTree)
		self.ui.listWidgetAll.itemDoubleClicked.connect(self.addFieldsHeatmap)
		self.ui.listWidgetSelected.itemDoubleClicked.connect(self.delFieldsHeatmap)
		self.ui.pushButtonClear.clicked.connect(self.clearCheck)
		# self.ui.listViewSpecificity.highlighted['QString'].connect(self.SpecSet)
		# self.ui.listViewSpecificity.mouseDoubleClickEvent.connect(self.SpecSet)

		self.ui.cboTreeOp1.id = 'TreeOp1'
		self.ui.cboTreeOp2.id = 'TreeOp2'
		self.ui.cboTreeOp3.id = 'TreeOp3'

		self.ui.lcdNumber_max.display(self.ui.horizontalScrollBar.maximum())
		self.ui.dial.setMaximum(self.ui.horizontalScrollBar.maximum())

		self.TextEdit = VGenesTextMain()
		# self.VGProgress = VGenesProgressBar()
		# self.ImportOptions =ImportDialogue()
		self.ImportOptions = ImportDataDialogue()
		self.AlterWindow = AlterDielog()
		self.AlterWindow.refreshDBSignal.connect(self.refreshDB)
		self.AlterWindow.intial = True
		self.AlterWindow.ui.pushButtonSave.setEnabled(False)

		self.ui.HTMLview = ResizeWidget(self)
		self.ui.gridLayoutStat.addWidget(self.ui.HTMLview, 2, 0, 10, 0)
		self.ui.HTMLview.resizeSignal.connect(self.resizeHTML)

		self.ui.HTMLviewClone = ResizeWidget(self)
		self.ui.gridLayoutClone.addWidget(self.ui.HTMLviewClone)
		self.ui.HTMLviewClone.resizeSignal.connect(self.resizeHTMLClone)

		self.enableEdit = False

		self.HeatmapList = []

	@pyqtSlot()
	def on_actionSave_As_triggered(self):
		self.saveAS()

	@pyqtSlot()
	def on_btnSaveAs_clicked(self):
		self.saveAS()

	def saveAS(self):
		options = QtWidgets.QFileDialog.Options()
		newDB, _ = QtWidgets.QFileDialog.getSaveFileName(self,
		                                                      "New Fragment Database",
		                                                      "New Fragment database",
		                                                      "VGene database Files (*.vdb);;All Files (*)",
		                                                      options=options)
		if newDB != 'none' and newDB != '':
			shutil.copyfile(DBFilename, newDB)
			QMessageBox.information(self, 'Information', 'DB saved to new location!',
			                    QMessageBox.Ok, QMessageBox.Ok)

	def addFieldsHeatmap(self):
		listItems = self.ui.listWidgetAll.selectedItems()
		for item in listItems:
			label = item.text()
			if label in self.HeatmapList:
				pass
			else:
				self.HeatmapList.append(label)
				self.ui.listWidgetSelected.addItem(label)

	def delFieldsHeatmap(self):
		listItems = self.ui.listWidgetSelected.selectedItems()
		for item in listItems:
			label = item.text()
			if label in self.HeatmapList:
				self.HeatmapList.remove(label)
				self.ui.listWidgetSelected.clear()
				self.ui.listWidgetSelected.addItems(self.HeatmapList)
			else:
				pass

	def updateClone(self):
		items = self.ui.listWidgetCloneMember.selectedItems()
		seq_name = ''
		for item in items:
			seq_name = item.text()

		SQLStatement = 'SELECT GeneType,V1,D1,J1,CDR3Length FROM vgenesDB WHERE SeqName = "' + seq_name + '"'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		seq_num = len(DataIn)

		if seq_num == 0:
			return

		record = DataIn[0]

		self.ui.lineEditCloneType.setText(record[0])

		self.ui.lineEditCloneV.setText(record[1])
		self.ui.lineEditCloneD.setText(record[2])
		self.ui.lineEditCloneJ.setText(record[3])
		self.ui.lineEditCloneCDR3len.setText(record[4])

	def checkClone(self):
		global MoveNotChange

		member_names = []
		n_member = self.ui.listWidgetCloneMember.count()
		for i in range(n_member):
			item = self.ui.listWidgetCloneMember.item(i)
			member_names.append(item.text())

		MoveNotChange = True
		rows = self.ui.SeqTable.rowCount()
		for row in range(rows):
			cur_name = self.ui.SeqTable.item(row, 1).text()
			if cur_name in member_names:
				self.ui.SeqTable.cellWidget(row, 0).setChecked(True)
			else:
				self.ui.SeqTable.cellWidget(row, 0).setChecked(False)
		MoveNotChange = False

		self.match_table_to_tree()

	def initial_Clone(self):
		# identify if clones exist
		SQLStatement = 'SELECT ClonalPool FROM vgenesDB'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		list1 = []
		for ele in DataIn:
			list1.append(ele[0])
		list_unique = list(set(list1))
		try:
			list_unique.remove('0')
		except:
			return

		if len(list_unique) == 0:
			return

		list_unique = [int(i) for i in list_unique]
		list_unique.sort()
		list_unique = ['Clone' + str(i) for i in list_unique]

		self.ui.listWidgetClone.clear()
		self.ui.listWidgetClone.addItems(list_unique)
		msg = 'Total ' + str(len(list_unique)) + ' Clones identified'
		self.ui.titleClone.setText(msg)

		if DBFilename != '' and DBFilename != None and DBFilename != 'none':
			fields_name = [""] + [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
			self.ui.comboBoxPieClone.clear()
			self.ui.comboBoxPieClone.addItems(fields_name)
			self.ui.comboBoxCol1Clone.clear()
			self.ui.comboBoxCol1Clone.addItems(fields_name)
			self.ui.comboBoxCol2Clone.clear()
			self.ui.comboBoxCol2Clone.addItems(fields_name)
			self.ui.comboBoxBoxDataClone.clear()
			self.ui.comboBoxBoxDataClone.addItems(fields_name)
			self.ui.comboBoxBox1Clone.clear()
			self.ui.comboBoxBox1Clone.addItems(fields_name)
			self.ui.comboBoxBox2Clone.clear()
			self.ui.comboBoxBox2Clone.addItems(fields_name)

	def selectClone(self):
		items = self.ui.listWidgetClone.selectedItems()
		for item in items:
			clone_id = item.text()
		clone_id = re.sub('Clone','',clone_id)

		SQLStatement = 'SELECT GeneType,V1,D1,J1,CDR3Length,SeqName FROM vgenesDB WHERE ClonalPool = "' + clone_id + '"'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		seq_num = len(DataIn)

		if seq_num == 0:
			return

		record = DataIn[0]

		self.ui.lineEditCloneName.setText('Clone' + clone_id)
		self.ui.lineEditCloneType.setText(record[0])
		self.ui.lineEditCloneNum.setText(str(seq_num))

		self.ui.lineEditCloneV.setText(record[1])
		self.ui.lineEditCloneD.setText(record[2])
		self.ui.lineEditCloneJ.setText(record[3])
		self.ui.lineEditCloneCDR3len.setText(record[4])

		mab_names = []
		for record in DataIn:
			mab_names.append(record[5])
		self.ui.listWidgetCloneMember.clear()
		self.ui.listWidgetCloneMember.addItems(mab_names)
		self.ui.listWidgetCloneMember.setCurrentItem(self.ui.listWidgetCloneMember.item(0))

	@pyqtSlot()
	def refreshDB(self):
		self.load_table()

		# update combox
		f_txt = re.sub(r'\(.+', '', self.ui.cboFindField.currentText())
		f_txt1 = re.sub(r'\(.+', '', self.ui.cboFindField1.currentText())
		tree_txt1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
		tree_txt2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
		tree_txt3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())

		self.ui.cboFindField.clear()
		self.ui.cboFindField1.clear()
		self.ui.cboTreeOp1.clear()
		self.ui.cboTreeOp2.clear()
		self.ui.cboTreeOp3.clear()

		field_list = [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		self.ui.cboFindField.addItems(field_list)
		self.ui.cboFindField1.addItems(field_list)
		self.ui.cboTreeOp1.addItems(field_list)
		self.ui.cboTreeOp2.addItems(field_list)
		self.ui.cboTreeOp3.addItems(field_list)
		self.ui.cboTreeOp1.addItem('None')
		self.ui.cboTreeOp2.addItem('None')
		self.ui.cboTreeOp3.addItem('None')

		index = FieldList.index(f_txt)
		self.ui.cboFindField.setCurrentText(FieldList[index] + '(' + RealNameList[index] + ')')
		index = FieldList.index(f_txt1)
		self.ui.cboFindField1.setCurrentText(FieldList[index] + '(' + RealNameList[index] + ')')

		if tree_txt1 == 'None':
			self.ui.cboTreeOp1.setCurrentText('None')
		else:
			index = FieldList.index(tree_txt1)
			self.ui.cboTreeOp1.setCurrentText(FieldList[index] + '(' + RealNameList[index] + ')')

		if tree_txt2 == 'None':
			self.ui.cboTreeOp2.setCurrentText('None')
		else:
			index = FieldList.index(tree_txt2)
			self.ui.cboTreeOp2.setCurrentText(FieldList[index] + '(' + RealNameList[index] + ')')

		if tree_txt3 == 'None':
			self.ui.cboTreeOp3.setCurrentText('None')
		else:
			index = FieldList.index(tree_txt3)
			self.ui.cboTreeOp3.setCurrentText(FieldList[index] + '(' + RealNameList[index] + ')')

	@pyqtSlot()
	def on_actionImport_Annotate_triggered(self):
		anno_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "select annotation file", '~/',
		                                                  "Comma Separated Values (*.csv);;All Files (*)")
		if anno_file == '' or anno_file == None:
			return

		self.annoDialog = AnnoDielog()
		self.annoDialog.refreshDBSignal.connect(self.refreshDB)
		Content = []

		try:
			csvFile = open(anno_file, "r")
			reader = csv.reader(csvFile)
			for item in reader:
				Content.append(item)
			csvFile.close()
		except:
			Msg = 'Can not process your file as CSV! Please check your input!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		self.annoDialog.Content = copy.deepcopy(Content)

		horizontalHeader = Content.pop(0)
		self.annoDialog.ui.comboBox.addItems(horizontalHeader)

		HEADERStatement = 'SELECT Field, FieldNickName FROM fieldsname'
		HeaderIn = VGenesSQL.RunSQL(DBFilename, HEADERStatement)
		Fields = [i[0] + '  (' + i[1] + ')' for i in HeaderIn]
		self.annoDialog.ui.comboBox2.addItems(Fields)

		num_col = len(horizontalHeader)
		num_row = len(Content)

		self.annoDialog.ui.tableWidget.setRowCount(num_row)
		self.annoDialog.ui.tableWidget.setColumnCount(num_col)
		self.annoDialog.ui.tableWidget.setHorizontalHeaderLabels(horizontalHeader)
		self.annoDialog.header = horizontalHeader

		for col_index in range(num_col):
			cell_comBox = QtWidgets.QComboBox()
			cell_comBox.addItems([''] + Fields)
			cell_comBox.setMaximumSize(10086,40)
			cell_comBox.setMinimumSize(50,20)
			cell_comBox.setEditable(True)
			self.annoDialog.ui.tableWidget.setCellWidget(0, col_index, cell_comBox)

		for row_index in range(num_row):
			for col_index in range(num_col):
				self.annoDialog.ui.tableWidget.setItem(row_index + 1, col_index,
				                                       QTableWidgetItem(Content[row_index][col_index]))

		self.annoDialog.show()

	@pyqtSlot()
	def on_CopyValue_clicked(self):
		self.copy_dialog = CopyDialog()
		field_list = [''] + [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		self.copy_dialog.ui.comboBoxFrom.addItems(field_list)
		self.copy_dialog.ui.comboBoxTo.addItems(field_list)
		self.copy_dialog.CopySignal.connect(self.CopyRecord)
		self.copy_dialog.exec_()

	def CopyRecord(self, field_from, field_to):
		global DontFindTwice
		field_from = re.sub(r'\(.+', '', field_from)
		field_to = re.sub(r'\(.+', '', field_to)
		SQLStatement = 'UPDATE vgenesDB SET ' + field_to + ' = ' + field_from
		try:
			VGenesSQL.RunUpdateSQL(DBFilename, SQLStatement)
			Msg = 'Update finished!'
			DontFindTwice = True
			self.refreshDB()
			value = self.ui.dial.value()
			self.updateF(value)
			DontFindTwice = False

			QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		except:
			Msg = 'Error occurs when updating the DB!\nCurrent SQL statement is:\n' + SQLStatement
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

	def openAlter(self):
		self.AlterWindow.ui.pushButtonSave.setEnabled(False)
		if self.AlterWindow.intial:
			if DBFilename != "":
				self.AlterWindow.ui.listWidget.addItems(FieldList)
				self.AlterWindow.ui.listWidget.setCurrentItem(self.AlterWindow.ui.listWidget.item(0))
				self.AlterWindow.ui.lineEditName.setText(FieldList[0])
				self.AlterWindow.ui.lineEditNickName.setText(RealNameList[0])
				self.AlterWindow.ui.lineEditType.setText(FieldTypeList[0])
				self.AlterWindow.ui.textEditNote.setText(FieldCommentList[0])
				self.AlterWindow.intial = False
				self.AlterWindow.ui.pushButtonSave.setDisabled(True)
		self.AlterWindow.show()

	def selectionMode(self):
		if self.ui.checkBoxRowSelection.isChecked():
			self.ui.SeqTable.setSelectionBehavior(QAbstractItemView.SelectRows)
		else:
			self.ui.SeqTable.setSelectionBehavior(QAbstractItemView.SelectItems)

	@pyqtSlot()
	def on_actionAlignmentHTML_triggered(self):
		global VGenesTextWindows
		# load data
		AlignIn = []
		listItems = self.getTreeCheckedChild()
		listItems = listItems[3]

		WhereState = ''
		NumSeqs = len(listItems)
		i = 1
		if len(listItems) == 0:
			QMessageBox.warning(self, 'Warning', 'Please select sequence from active sequence panel!',
			                    QMessageBox.Ok,
			                    QMessageBox.Ok)
			return
		for item in listItems:
			WhereState += 'SeqName = "' + item + '"'
			if NumSeqs > i:
				WhereState += ' OR '
			i += 1

		SQLStatement = 'SELECT SeqName, Sequence FROM vgenesDB WHERE ' + WhereState
		DataIn =  VGenesSQL.RunSQL(DBFilename, SQLStatement)

		for item in DataIn:
			SeqName = item[0]
			Sequence = item[1]
			Sequence = Sequence.replace("-","")
			Sequence = Sequence.upper()
			EachIn = (SeqName, Sequence)
			AlignIn.append(EachIn)
		# make HTML
		html_file = AlignSequencesHTML(AlignIn, '')
		if html_file[0] == 'W':
			QMessageBox.warning(self, 'Warning', html_file, QMessageBox.Ok, QMessageBox.Ok)
			return
		# delete close window objects
		del_list = []
		for id, obj in VGenesTextWindows.items():
			if obj.isVisible() == False:
				del_list.append(id)
		for id in del_list:
			del_obj = VGenesTextWindows.pop(id)

		# display
		window_id = int(time.time() * 100)
		VGenesTextWindows[window_id] = htmlDialog()
		VGenesTextWindows[window_id].id = window_id
		layout = QGridLayout(VGenesTextWindows[window_id])
		view = QWebEngineView(self)
		view.load(QUrl("file://" + html_file))
		view.show()
		layout.addWidget(view)
		VGenesTextWindows[window_id].show()

	def loadNewickTree(self):
		treefile, _ = QtWidgets.QFileDialog.getOpenFileName(self, "select Newick file", '~/',
		                                                  "Newick File (*.nwk);;All Files (*)")
		if treefile == '' or treefile == None:
			return
		self.ui.lineEditNewick.setText(treefile)

		# generate html page
		treefile = treefile
		f = open(treefile, 'r')
		tree_str = f.readline()
		f.close()
		tree_str = 'var test_string = "' + tree_str.rstrip("\n") + '";\n'

		out_html_file = os.path.join(temp_folder, 'newick_tree.html')
		header_file = os.path.join(working_prefix, 'Data', 'template_newick_tree.html')
		shutil.copyfile(header_file, out_html_file)

		foot = 'var container_id = "#tree_container";\nvar svg = d3.select(container_id).append("svg")' \
		       '.attr("width", width).attr("height", height);\n$( document ).ready( function () {' \
		       'default_tree_settings();tree(test_string).svg (svg).layout();update_selection_names();' \
		       '});\n</script>\n</body>\n</html>'
		out_file_handle = open(out_html_file, 'a')
		out_file_handle.write(tree_str)
		out_file_handle.write(foot)
		out_file_handle.close()

		if self.ui.radioButtonNewick.isChecked():
			# display
			window_id = int(time.time() * 100)
			VGenesTextWindows[window_id] = htmlDialog()
			VGenesTextWindows[window_id].id = window_id
			layout = QGridLayout(VGenesTextWindows[window_id])
			view = QWebEngineView(self)
			view.load(QUrl("file://" + out_html_file))
			view.show()
			layout.addWidget(view)
			VGenesTextWindows[window_id].show()
		else:
			# display
			view = QWebEngineView()
			view.load(QUrl("file://" + out_html_file))
			view.show()

			layout = self.ui.groupBoxTree.layout()
			if layout == None:
				layout = QGridLayout(self.ui.groupBoxTree)
			else:
				for i in range(layout.count()):
					layout.removeWidget(layout.itemAt(i).widget())
			layout.addWidget(view)

	def buildCloneTree(self):
		clone_name = self.ui.comboBoxTree.currentText()
		clone_name = re.sub('Clone','',clone_name)

		WHEREStatement = 'WHERE ClonalPool = "' + clone_name + '"'
		SQLStatement = 'SELECT SeqName,Sequence,GermlineSequence FROM vgenesDB ' + WHEREStatement
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		if len(DataIn) < 3:
			Msg = 'Too few sequences in this clone! At least three sequences required!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		# write sequences into file
		time_stamp = str(int(time.time() * 100))
		this_folder = os.path.join(temp_folder, time_stamp)
		cmd = 'mkdir ' + this_folder
		try:
			os.system(cmd)
		except:
			QMessageBox.warning(self, 'Warning', 'Fail to make temp folder!', QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

		# alignment
		aafilename = this_folder + "/input.fas"
		outfilename = this_folder + "/alignment.fas"
		treefilename = 'tree'
		out_handle = open(aafilename, 'w')

		out_handle.write('>Germline\n')
		out_handle.write(DataIn[0][2] + '\n')
		for item in DataIn:
			SeqName = item[0]
			Sequence = item[1]

			# parse seq name
			SeqName = re.sub(r'[^\w\d\/\>]', '_', SeqName)
			SeqName = re.sub(r'_+', '_', SeqName)
			SeqName = SeqName.strip('_')

			out_handle.write('>' + SeqName + '\n')
			out_handle.write(Sequence + '\n')
		out_handle.close()

		# alignment
		cmd = muscle_path
		cmd += " -in " + aafilename + " -out " + outfilename
		try:
			os.system(cmd)
		except:
			QMessageBox.warning(self, 'Warning', 'Fail to run muscle! Check your muscle path!', QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

		# generate tree
		cmd = 'cd ' + this_folder + ';'
		cmd += raxml_path
		cmd += ' -m GTRGAMMA -p 12345 -T 2 -s ' + outfilename + ' -n ' + treefilename

		os.system(cmd)
		print("tree done!")

		# generate html page
		treefile = os.path.join(this_folder, 'RAxML_bestTree.tree')
		f = open(treefile, 'r')
		tree_str = f.readline()
		f.close()
		tree_str = 'var test_string = "' + tree_str.rstrip("\n") + '";\n'

		out_html_file = os.path.join(this_folder, 'tree.html')
		header_file = os.path.join(working_prefix, 'Data', 'template_raxml_tree.html')
		shutil.copyfile(header_file, out_html_file)

		foot = 'var container_id = "#tree_container";\nvar svg = d3.select(container_id).append("svg")' \
		       '.attr("width", width).attr("height", height);\n$( document ).ready( function () {' \
		       'default_tree_settings();tree(test_string).svg (svg).layout();update_selection_names();' \
		       '});\n</script>\n</body>\n</html>'
		out_file_handle = open(out_html_file, 'a')
		out_file_handle.write(tree_str)
		out_file_handle.write(foot)
		out_file_handle.close()

		if self.ui.radioButtonNewick.isChecked():
			# display
			window_id = int(time.time() * 100)
			VGenesTextWindows[window_id] = htmlDialog()
			VGenesTextWindows[window_id].id = window_id
			layout = QGridLayout(VGenesTextWindows[window_id])
			view = QWebEngineView(self)
			view.load(QUrl("file://" + out_html_file))
			view.show()
			layout.addWidget(view)
			VGenesTextWindows[window_id].show()
		else:
			# display
			view = QWebEngineView()
			view.load(QUrl("file://" + out_html_file))
			view.show()

			layout = self.ui.groupBoxTree.layout()
			if layout == None:
				layout = QGridLayout(self.ui.groupBoxTree)
			else:
				for i in range(layout.count()):
					layout.removeWidget(layout.itemAt(i).widget())
			layout.addWidget(view)

	def loadIgphyml(self):
		ig_out, _ = QtWidgets.QFileDialog.getOpenFileName(self, "select igphyml output", '~/',
		                                         "igphyml output File (*.tab);;All Files (*)")
		if ig_out == '' or ig_out == None:
			return
		self.ui.igphyml_line.setText(ig_out)

		trees = []
		f = open(ig_out, 'r')
		lines = f.readlines()
		f.close()
		if len(lines) < 3:
			QMessageBox.warning(self, 'Warning', 'No trees detected!',QMessageBox.Ok,QMessageBox.Ok)
			return
		lines = lines[2:]
		for line in lines:
			tmp_list = line.split('\t')
			this = ('Clone ' + tmp_list[0],tmp_list[14])
			trees.append(this)

		out_html_file = os.path.join(temp_folder, 'tree.html')
		header_file = os.path.join(working_prefix, 'Data', 'template_tree.html')
		shutil.copyfile(header_file, out_html_file)

		tree_str = 'var trees = new Array();\n'
		i = 0
		for tree in trees:
			tree_str = tree_str + '$("#mySelect").append("<option value=' + "'" + str(i) + "'" + '>' + tree[0] + '</option>");\n'
			tree_str = tree_str + 'trees[' + str(i) + ']="' + tree[1].strip('\n') + '";\n'
			i += 1

		tree_str = tree_str + 'var test_string = trees[0];\n'

		foot = 'var container_id = "#tree_container";\nvar svg = d3.select(container_id).append("svg")' \
		       '.attr("width", width).attr("height", height);\n$( document ).ready( function () {' \
		       'default_tree_settings();tree(test_string).svg (svg).layout();update_selection_names();' \
		       '});\n</script>\n</body>\n</html>'
		out_file_handle = open(out_html_file, 'a')
		out_file_handle.write(tree_str)
		out_file_handle.write(foot)
		out_file_handle.close()

		print("html done!")

		if self.ui.radioButtonNewick.isChecked():
			# display
			window_id = int(time.time() * 100)
			VGenesTextWindows[window_id] = htmlDialog()
			VGenesTextWindows[window_id].id = window_id
			layout = QGridLayout(VGenesTextWindows[window_id])
			view = QWebEngineView(self)
			view.load(QUrl("file://" + out_html_file))
			view.show()
			layout.addWidget(view)
			VGenesTextWindows[window_id].show()
		else:
			# display
			view = QWebEngineView()
			view.load(QUrl("file://" + out_html_file))
			view.show()

			layout = self.ui.groupBoxTree.layout()
			if layout == None:
				layout = QGridLayout(self.ui.groupBoxTree)
			else:
				for i in range(layout.count()):
					layout.removeWidget(layout.itemAt(i).widget())
			layout.addWidget(view)

	def makeNTLogo(self):
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		listItems = self.getTreeCheckedChild()
		listItems = listItems[3]
		WhereState = ''
		NumSeqs = len(listItems)
		# if not listItems: do nothing
		DataSet = []
		if NumSeqs < 1:
			msg = 'Please select(check) at least one sequence!'
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		else:
			i = 1
			for item in listItems:
				WhereState += 'SeqName = "' + item + '"'
				if NumSeqs > i:
					WhereState += ' OR '
				i += 1

			field = self.ui.comboBoxFieldLogo.currentText()
			SQLStatement = 'SELECT SeqName, ' + field + ' FROM vgenesDB WHERE ' + WhereState
			DataIn =  VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for item in DataIn:
				SeqName = item[0]
				Sequence = item[1]

				Sequence = Sequence.upper()
				EachIn = (SeqName, Sequence)
				DataSet.append(EachIn)

		# align selected sequences using ClustalOmega
		outfilename = ''
		try:
			if len(DataSet) == 1:
				time_stamp = str(int(time.time() * 100))
				outfilename = os.path.join(temp_folder, "out-" + time_stamp + ".fas")
				out_handle = open(outfilename, 'w')
				out_handle.write('>' + DataSet[0][0] + '\n')
				out_handle.write(DataSet[0][1])
				out_handle.close()
			else:
				if os.path.exists(clustal_path):
					outfilename = VGenesSeq.ClustalO_new(DataSet, 80, True, temp_folder, clustal_path)
				else:
					QMessageBox.warning(self, 'Warning',
					                    'The Clustal Omega does not exist! Check your path!', QMessageBox.Ok,
					                    QMessageBox.Ok)
					return
		except:
			return

		# start web logo
		f = open(outfilename)
		seqs = read_seq_data(f)
		data = LogoData.from_seqs(seqs)

		options = LogoOptions()
		options.fineprint = 'VGene Generated by WebLogo 3.7'
		format = LogoFormat(data, options)

		time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())

		if self.ui.radioButtonPop.isChecked():
			eps = eps_formatter(data, format)
			out_eps = os.path.join(temp_folder, "out-" + time_stamp + ".eps")
			with open(out_eps, 'wb') as f:
				f.write(eps)
			cmd = "open " + out_eps
			bot1 = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True,
			             env={"LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8"})
		else:
			try:
				svg = svg_formatter(data, format)
				svg = svg.decode("utf-8")

				out_svg = os.path.join(temp_folder, "out-" + time_stamp + ".html")
				with open(out_svg, 'w') as f:
					f.write('<!DOCTYPE html>\n<html>\n<body style="margin-left: 0px;\n">')
					f.write(svg)
					f.write('\n</body>\n</html>')

				# display
				view = QWebEngineView()
				view.load(QUrl("file://" + out_svg))
				view.show()

				layout = self.ui.groupBoxLogo.layout()
				if layout == None:
					layout = QGridLayout(self.ui.groupBoxLogo)
				else:
					for i in range(layout.count()):
						layout.removeWidget(layout.itemAt(i).widget())
				layout.addWidget(view)
			except:
				eps = eps_formatter(data, format)
				out_eps = os.path.join(temp_folder, "out-" + time_stamp + ".eps")
				with open(out_eps, 'wb') as f:
					f.write(eps)
				cmd = "open " + out_eps
				bot1 = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True,
				             env={"LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8"})

	def makeAALogo(self):
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		listItems = self.getTreeCheckedChild()
		listItems = listItems[3]
		WhereState = ''
		NumSeqs = len(listItems)
		# if not listItems: do nothing
		DataSet = []
		if NumSeqs < 1:
			msg = 'Please select(check) at least one sequence!'
			QMessageBox.warning(self, 'Warning', msg, QMessageBox.Ok, QMessageBox.Ok)
			return
		else:
			i = 1
			for item in listItems:
				WhereState += 'SeqName = "' + item + '"'
				if NumSeqs > i:
					WhereState += ' OR '
				i += 1

			field = self.ui.comboBoxFieldLogo.currentText()
			SQLStatement = 'SELECT SeqName, ' + field + ' FROM vgenesDB WHERE ' + WhereState
			DataIn =  VGenesSQL.RunSQL(DBFilename, SQLStatement)

			for item in DataIn:
				SeqName = item[0]
				Sequence = item[1]

				AAseq, msg = VGenesSeq.Translator(Sequence, 0)
				AAseq = AAseq.upper()
				EachIn = (SeqName, AAseq)
				DataSet.append(EachIn)

		# align selected sequences using ClustalOmega
		outfilename = ''
		try:
			if len(DataSet) == 1:
				time_stamp = str(int(time.time() * 100))
				outfilename = os.path.join(temp_folder, "out-" + time_stamp + ".fas")
				out_handle = open(outfilename, 'w')
				out_handle.write('>' + DataSet[0][0] + '\n')
				out_handle.write(DataSet[0][1])
				out_handle.close()
			else:
				if os.path.exists(clustal_path):
					outfilename = VGenesSeq.ClustalO_new(DataSet, 80, True, temp_folder, clustal_path)
				else:
					QMessageBox.warning(self, 'Warning',
					                    'The Clustal Omega does not exist! Check your path!', QMessageBox.Ok,
					                    QMessageBox.Ok)
					return
		except:
			return

		# start web logo
		f = open(outfilename)
		seqs = read_seq_data(f)
		data = LogoData.from_seqs(seqs)

		options = LogoOptions()
		options.fineprint = 'VGene Generated by WebLogo 3.7'
		format = LogoFormat(data, options)

		time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())

		if self.ui.radioButtonPop.isChecked():
			eps = eps_formatter(data, format)
			out_eps = os.path.join(temp_folder, "out-" + time_stamp + ".eps")
			with open(out_eps, 'wb') as f:
				f.write(eps)
			cmd = "open " + out_eps
			bot1 = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True,
			             env={"LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8"})
		else:

			try:
				svg = svg_formatter(data, format)
				svg = svg.decode("utf-8")

				out_svg = os.path.join(temp_folder, "out-" + time_stamp + ".html")
				with open(out_svg, 'w') as f:
					f.write('<!DOCTYPE html>\n<html>\n<body style="margin-left: 0px;\n">')
					f.write(svg)
					f.write('\n</body>\n</html>')

				# display
				view = QWebEngineView()
				view.load(QUrl("file://" + out_svg))
				view.show()

				layout = self.ui.groupBoxLogo.layout()
				if layout == None:
					layout = QGridLayout(self.ui.groupBoxLogo)
				else:
					for i in range(layout.count()):
						layout.removeWidget(layout.itemAt(i).widget())
				layout.addWidget(view)
			except:
				eps = eps_formatter(data, format)
				out_eps = os.path.join(temp_folder, "out-" + time_stamp + ".eps")
				with open(out_eps, 'wb') as f:
					f.write(eps)
				cmd = "open " + out_eps
				bot1 = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True,
				             env={"LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8"})

	def InitialGraphic(self):
		global DBFilename
		global updateMarker
		if DBFilename != '' and DBFilename != 'none' and DBFilename != None:
			#Msg = str(self.ui.tabWidget.currentIndex())
			#QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			if self.ui.tabWidget.currentIndex() == 6:
				if self.ui.comboBoxPie.count() == 0:
					# setup options for graph table
					fields_name = [''] + [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]

					self.ui.comboBoxPie.clear()
					self.ui.comboBoxPie.addItems(fields_name)
					self.ui.comboBoxCol1.clear()
					self.ui.comboBoxCol1.addItems(fields_name)
					self.ui.comboBoxCol2.clear()
					self.ui.comboBoxCol2.addItems(fields_name)
					self.ui.comboBoxBoxData.clear()
					self.ui.comboBoxBoxData.addItems(fields_name)
					self.ui.comboBoxBox1.clear()
					self.ui.comboBoxBox1.addItems(fields_name)
					self.ui.comboBoxBox2.clear()
					self.ui.comboBoxBox2.addItems(fields_name)
					self.ui.comboBoxRiver1.clear()
					self.ui.comboBoxRiver1.addItems(fields_name)
					self.ui.comboBoxRiver2.clear()
					self.ui.comboBoxRiver2.addItems(fields_name)
					self.ui.comboBoxWord.clear()
					self.ui.comboBoxWord.addItems(fields_name)
					self.ui.comboBoxScatterX.clear()
					self.ui.comboBoxScatterX.addItems(fields_name)
					self.ui.comboBoxScatterY.clear()
					self.ui.comboBoxScatterY.addItems(fields_name)
					self.ui.comboBoxScatterGroup.clear()
					self.ui.comboBoxScatterGroup.addItems(fields_name)
					self.ui.comboBoxTree1.clear()
					self.ui.comboBoxTree1.addItems(fields_name)
					self.ui.comboBoxTree2.clear()
					self.ui.comboBoxTree2.addItems(fields_name)
					self.ui.comboBoxTree3.clear()
					self.ui.comboBoxTree3.addItems(fields_name)
					self.ui.listWidgetAll.clear()
					self.ui.listWidgetAll.addItems(fields_name[1:])
					self.ui.comboBoxSortField.clear()
					self.ui.comboBoxSortField.addItems(fields_name)
			elif self.ui.tabWidget.currentIndex() == 8:
				SQLStatement = 'SELECT ClonalPool FROM vgenesDB'
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				if len(DataIn) > 1:
					list1 = []
					for ele in DataIn:
						list1.append(ele[0])
					list_unique = list(set(list1))
					if '0' in list_unique:
						list_unique.remove('0')
					if '' in list_unique:
						list_unique.remove('')

					if len(list_unique) > 0:
						list_unique = [int(i) for i in list_unique]
						list_unique.sort()
						list_unique = ['Clone' + str(i) for i in list_unique]

						self.ui.comboBoxTree.clear()
						self.ui.comboBoxTree.addItems(list_unique)
			elif self.ui.tabWidget.currentIndex() == 0:
				# if old table exists, clear table
				if self.ui.SeqTable.columnCount() > 0:
					if updateMarker == True:
						self.load_table()
						self.match_tree_to_table()
						self.tree_to_table_selection()
						updateMarker = False
				else:
					self.load_table()
					self.match_tree_to_table()
					self.tree_to_table_selection()

	def load_table(self):
		if self.ui.SeqTable.columnCount() > 0:
			self.ui.SeqTable.itemChanged.disconnect(self.EditTableItem)
		self.ui.SeqTable.setColumnCount(0)
		self.ui.SeqTable.setRowCount(0)

		a = DBFilename

		a = FieldList
		b = RealNameList
		c = FieldTypeList
		d = FieldCommentList



		# load data for new table
		if DBFilename != '' and DBFilename != 'none' and DBFilename != None:
			field1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
			field2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
			field3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())

			fields = ','.join(FieldList)
			SQLStatement = 'select '+ fields +' from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', ' + field3 + ', SeqName'
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)


			num_row = len(DataIn)
			num_col = len(RealNameList)
			self.ui.SeqTable.setRowCount(num_row)
			self.ui.SeqTable.setColumnCount(num_col + 1)

			horizontalHeader = [''] + RealNameList
			self.ui.SeqTable.setHorizontalHeaderLabels(horizontalHeader)
			self.ui.SeqTable.fields = horizontalHeader
			# re-size column size
			self.ui.SeqTable.horizontalHeader().resizeSection(0, 10)
			self.ui.SeqTable.setSelectionMode(QAbstractItemView.ExtendedSelection)
			if self.ui.checkBoxRowSelection.isChecked():
				self.ui.SeqTable.setSelectionBehavior(QAbstractItemView.SelectRows)
			else:
				self.ui.SeqTable.setSelectionBehavior(QAbstractItemView.SelectItems)

			for row_index in range(num_row):
				cell_checkBox = QCheckBox()
				# cell_checkBox.setText(DataIn[row_index][0])
				cell_checkBox.setChecked(False)
				cell_checkBox.stateChanged.connect(self.multipleSelection)
				self.ui.SeqTable.setCellWidget(row_index, 0, cell_checkBox)

				for col_index in range(num_col):
					unit = QTableWidgetItem(DataIn[row_index][col_index])
					unit.last_name = DataIn[row_index][col_index]
					self.ui.SeqTable.setItem(row_index, col_index + 1, unit)

			# disable edit
			self.ui.SeqTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
			# show sort indicator
			self.ui.SeqTable.horizontalHeader().setSortIndicatorShown(True)
			# connect sort indicator to slot function
			self.ui.SeqTable.horizontalHeader().sectionClicked.connect(self.sortTable)
			self.ui.SeqTable.itemChanged.connect(self.EditTableItem)

	def clearCheck(self):
		self.ui.checkBoxAll1.setChecked(False)
		self.ui.checkBoxAll.setChecked(False)

		self.checkAll()

	def checkAll(self):
		global MoveNotChange
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		MoveNotChange = True
		rows = self.ui.SeqTable.rowCount()
		if self.ui.checkBoxAll.isChecked():
			self.ui.checkBoxAll1.setChecked(True)
			for row in range(rows):
				self.ui.SeqTable.cellWidget(row, 0).setChecked(True)
		else:
			self.ui.checkBoxAll1.setChecked(False)
			for row in range(rows):
				self.ui.SeqTable.cellWidget(row, 0).setChecked(False)
		MoveNotChange = False
		self.match_table_to_tree()

	def checkAll1(self):
		global MoveNotChange
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		MoveNotChange = True
		rows = self.ui.SeqTable.rowCount()
		if self.ui.checkBoxAll1.isChecked():
			self.ui.checkBoxAll.setChecked(True)
			for row in range(rows):
				self.ui.SeqTable.cellWidget(row, 0).setChecked(True)
		else:
			self.ui.checkBoxAll.setChecked(False)
			for row in range(rows):
				self.ui.SeqTable.cellWidget(row, 0).setChecked(False)
		MoveNotChange = False
		self.match_table_to_tree()

	def multipleSelection(self, int_signal):
		global MoveNotChange
		if MoveNotChange:
			return

		MoveNotChange = True
		buttonClicked = self.sender()
		postitionOfWidget = buttonClicked.pos()
		index = self.ui.SeqTable.indexAt(postitionOfWidget)
		row = index.row()

		# get all selected rows
		rows = []
		item = self.ui.SeqTable.selectedItems()
		for i in item:
			if self.ui.SeqTable.indexFromItem(i).row() not in rows:
				rows.append(self.ui.SeqTable.indexFromItem(i).row())

		DataIn = []
		# if current row in selected rows, check all rows
		if row in rows:
			for cur_row in rows:
				if cur_row != row:
					if int_signal == 2:		# signal = 2 means check event
						self.ui.SeqTable.cellWidget(cur_row, 0).setChecked(True)
					else:
						self.ui.SeqTable.cellWidget(cur_row, 0).setChecked(False)
		MoveNotChange = False

		self.match_table_to_tree()

	def match_tree_to_table(self):
		global TreeSelected
		global MoveNotChange
		# check if checked item changed
		selected_list = self.getTreeCheckedChild()
		selected_list = selected_list[3]

		# if TreeSelected.sort() == selected_list.sort():
		#	print(",".join(selected_list))
		#	print(",".join(TreeSelected))
		#	return
		# else:
		MoveNotChange = True
		rows = self.ui.SeqTable.rowCount()
		for row in range(rows):
			cur_name = self.ui.SeqTable.item(row, 1).text()
			if cur_name in selected_list:
				self.ui.SeqTable.cellWidget(row, 0).setChecked(True)
			else:
				self.ui.SeqTable.cellWidget(row, 0).setChecked(False)
		MoveNotChange = False

	def match_table_to_tree(self):
		DataIn = []
		# get all checked table rows
		total_rows = self.ui.SeqTable.rowCount()
		for row in range(total_rows):
			if self.ui.SeqTable.cellWidget(row, 0).isChecked():
				DataIn.append(self.ui.SeqTable.item(row, 1).text())

		# update the selection
		self.clearTreeChecks()
		NumFound = len(DataIn)
		i = 0
		for Seqname in DataIn:
			found = self.ui.treeWidget.findItems(Seqname, Qt.MatchRecursive, 0)
			i += 1
			for record in found:
				if i == NumFound - 1:
					wasClicked = True
				record.setCheckState(0, Qt.Checked)

	def EditTableItem(self, item):
		global MoveNotChange
		global NameIndex

		if MoveNotChange:
			return

		row = item.row()
		col = item.column()
		CurVal = item.text()

		horizontalHeader = self.ui.SeqTable.fields
		col_name = horizontalHeader[col]
		index = RealNameList.index(col_name)
		col_name = FieldList[index]

		if col == 1:  # update sequence name
			SeqName = item.last_name
		else:
			SeqName = self.ui.SeqTable.item(row, 1).text()

		try:
			self.UpdateSeq(SeqName, CurVal, col_name)
			if col == 1:  # update sequence name
				item.last_name = CurVal

				# update name index
				#NameIndex[CurVal] = NameIndex[SeqName]
				#del NameIndex[SeqName]
				#refresh model
				#model = self.ui.tableView.model()
				#model.refresh()
				# update tree
				#SQLFields = (
				#self.ui.cboTreeOp1.currentText(), self.ui.cboTreeOp2.currentText(), self.ui.cboTreeOp3.currentText())
				'''
				index1 = RealNameList.index(self.ui.cboTreeOp1.currentText())
				index2 = RealNameList.index(self.ui.cboTreeOp2.currentText())
				index3 = RealNameList.index(self.ui.cboTreeOp3.currentText())

				SQLFields = (FieldList[index1], FieldList[index2], FieldList[index3])
				'''
				SQLFields = (
					re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText()),
					re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText()),
					re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
				)

				self.initializeTreeView(SQLFields)
				self.ui.treeWidget.expandAll()
		except:
			MoveNotChange = True
			col = item.column()
			self.ui.SeqTable.item(row, col).setText(SeqName)
			MoveNotChange = False
			QMessageBox.warning(self, 'Warning',
								'The name:\n' + CurVal + '\nhas been taken! Please choose another name!',
								QMessageBox.Ok, QMessageBox.Ok)

	@pyqtSlot()
	def UpdateSeq(self, ID, ItemValue, FieldName):
		global DBFilename
		# ID = item[0]
		VGenesSQL.UpdateFieldbySeqName(ID, ItemValue, FieldName, DBFilename)

	def sortTable(self, index):
		if self.ui.tabWidget.currentIndex() == 9:
			self.ui.SeqTable.sortByColumn(index, self.ui.SeqTable.horizontalHeader().sortIndicatorOrder())
		else:
			self.ui.SeqTable.sortByColumn(index, self.ui.SeqTable.horizontalHeader().sortIndicatorOrder())

	def ChangeEditMode(self):
		if self.ui.SeqTable.editTriggers() == QtWidgets.QAbstractItemView.NoEditTriggers:
			unlock_icon = QIcon()
			unlock_icon.addPixmap(QPixmap(":/PNG-Icons/unlocked.png"), QIcon.Normal, QIcon.Off)
			self.ui.EditLock.setIcon(unlock_icon)
			self.ui.EditLock.setText('Edit Lock: Unlock (Double click to edit)')
			self.ui.SeqTable.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked)
		else:
			lock_icon = QIcon()
			lock_icon.addPixmap(QPixmap(":/PNG-Icons/locked.png"), QIcon.Normal, QIcon.Off)
			self.ui.EditLock.setIcon(lock_icon)
			self.ui.EditLock.setText('Edit Lock: Locked')
			self.ui.SeqTable.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

	def GenerateFigureClone(self):
		msg = 'h:' + str(self.ui.HTMLviewClone.h) + ',w:' +str(self.ui.HTMLviewClone.w)
		print(msg)

		global DBFilename
		# global data

		# select data or not
		clone_id = self.ui.lineEditCloneName.text()
		clone_id = re.sub('Clone','',clone_id)
		try:
			int(clone_id)
		except:
			return
		where_statement = 'WHERE ClonalPool = "' + clone_id + '"'

		# pie chart
		if self.ui.tabWidgetClone.currentIndex() == 0:
			# get data
			field = re.sub(r'\(.+', '', self.ui.comboBoxPieClone.currentText())
			if field == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if len(DataIn) == 0:
				Msg = 'No records can be fetched!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

			data = []
			for element in DataIn:
				data.append(element[0])
			result = Counter(data)
			labels = result.keys()
			values = result.values()

			my_pyecharts = (
				Pie(init_opts=opts.InitOpts(width=str(self.ui.HTMLviewClone.w)+"px", height=str(self.ui.HTMLviewClone.h)+"px", renderer='svg'))\
				.add('', [list(z) for z in zip(labels, values)], radius=["40%", "75%"])\
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show=self.ui.checkBoxFigLegendClone.isChecked()
					),
				)
				.set_series_opts(label_opts=opts.LabelOpts(formatter=" {b}: {c} ({d}%)"))
			)  #
		# Bar chart
		elif self.ui.tabWidgetClone.currentIndex() == 1:
			# get data
			field1 = re.sub(r'\(.+', '', self.ui.comboBoxCol1Clone.currentText())
			field2 = re.sub(r'\(.+', '', self.ui.comboBoxCol2Clone.currentText())
			if field1 == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			if field2 == field1:
				QMessageBox.warning(self, 'Warning', 'Please select different group factors for field1 and field2!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

			multi_factor = False
			if field2 == "":
				field = field1
			else:
				field = field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if multi_factor == True:
				if self.ui.checkBoxStackClone.isChecked():
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

				my_bar = Bar(init_opts=opts.InitOpts(width=str(self.ui.HTMLviewClone.w)+"px", height=str(self.ui.HTMLviewClone.h)+"px", renderer='svg')) \
					.add_xaxis(labels) \
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
						),
					)

				for group in dic_keys:
					cur_data = data[group]
					group_data = []
					for ele in labels:
						group_data.append(cur_data.count(ele))
					my_bar.add_yaxis(group, group_data, stack=stack)
				my_bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False, formatter=" {b}: {c}"))

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
					Bar(init_opts=opts.InitOpts(width=str(self.ui.HTMLviewClone.w)+"px", height=str(self.ui.HTMLviewClone.h)+"px", renderer='svg'))
					.add_xaxis(labels)
					.add_yaxis(field1, values)
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegendClone.isChecked()
						),
					)
					.set_series_opts(label_opts=opts.LabelOpts(is_show=False, formatter=" {b}: {c}"))
				)
		# Box plot
		elif self.ui.tabWidgetClone.currentIndex() == 2:
			# get data
			data_field = re.sub(r'\(.+', '', self.ui.comboBoxBoxDataClone.currentText())
			field1 = re.sub(r'\(.+', '', self.ui.comboBoxBox1Clone.currentText())
			field2 = re.sub(r'\(.+', '', self.ui.comboBoxBox2Clone.currentText())
			if field1 == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			if field2 == field1:
				QMessageBox.warning(self, 'Warning', 'Please select different group factors for field1 and field2!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			multi_factor = False
			if field2 == "":
				field = data_field + "," + field1
			else:
				field = data_field + "," + field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			box_data = [i[0] for i in DataIn]
			try:
				box_data = list(map(float, box_data))
			except:
				QMessageBox.warning(self, 'Warning', 'The data field is not numerical! Check your input!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

			if min(box_data) >= 0:
				null_data = [0, 0, 0, 0, 0]
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

				my_bar = Boxplot(init_opts=opts.InitOpts(width=str(self.ui.HTMLviewClone.w)+"px", height=str(self.ui.HTMLviewClone.h)+"px", renderer='svg')) \
					.add_xaxis(labels) \
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
						),
						# toolbox_opts = opts.ToolboxOpts()
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
						cur_data = []
						for i in my_dict[ele]:
							cur_data.append(box_data[i])
						data_v1.append(cur_data)
					else:
						data_v1.append(null_data)

				my_pyecharts = (
					Boxplot(init_opts=opts.InitOpts(width=str(self.ui.HTMLviewClone.w)+"px", height=str(self.ui.HTMLviewClone.h)+"px", renderer='svg'))
					.add_xaxis(labels)
					.add_yaxis(field1, Boxplot.prepare_data(data_v1))
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
						),
					)
				)

		# load figure
		html_path = os.path.join(temp_folder, 'figure_clone.html')
		my_pyecharts.render(path=html_path)
		# adjust the window size seting
		file_handle = open(html_path, 'r')
		lines = file_handle.readlines()
		file_handle.close()
		# edit js line
		js_line = '<script type="text/javascript" src="' + \
		          os.path.join(js_folder, 'echarts.js') + '"></script>' + \
		          '<script src="' + os.path.join(js_folder, 'jquery.js') + '"></script>' + \
		          '<script src="qrc:///qtwebchannel/qwebchannel.js"></script>'
		lines[5] = js_line
		# edit style line
		style_line = lines[9]
		style_pos = style_line.find('style')
		style_line = style_line[
		             0:style_pos] + 'style="position: fixed; top: 0px; left: 5%;width:90%; height:' + str(
			self.ui.HTMLviewClone.h - 20) + 'px;"></div>'
		lines[9] = style_line
		insert_js = '<script type="text/javascript">$(document).ready(function() {' \
		            'new QWebChannel(qt.webChannelTransport, function(channel) {' \
		            'var my_object = channel.objects.connection;$("#download").click(function(){' \
		            'my_object.download(text);});$("#update").click(function(){' \
		            'my_object.updateSelection(text);});});});</script>'
		insert_btn = '<input id="download" type="button" value="" style="display:none;"/>' \
		             '<input id="update" type="button" value="" style="display:none;"/>'
		lines = lines[:6] + [insert_js] + lines[6:9] + [insert_btn] + lines[9:]

		if self.ui.tabWidgetFig.currentIndex() in [0, 1, 2, 6]:
			# insert click response function
			echart_init_line = lines[13]
			matchObj = re.match(r'.+var\s(\S+)\s=', echart_init_line)
			chart_id = matchObj.group(1)
			js_cmd = chart_id + ".on('click', function (params) {" \
			                    "if(params.data['0'] == null){text = params.name + ',' + params.seriesName + ',0,0';}" \
			                    "else{text = params.name + ',' + params.seriesName + ','+params.data['0']+','+params.data['1'];}" \
			                    "$('#update').click();});"
			lines = lines[:-3] + [js_cmd] + lines[-3:]

		content = '\n'.join(lines)
		file_handle = open(html_path, 'w')
		file_handle.write(content)
		file_handle.close()
		# show local HTML
		self.ui.HTMLviewClone.load(QUrl('file://' + html_path))
		self.ui.HTMLviewClone.show()
		self.ui.HTMLviewClone.html = "loaded"
		self.ui.HTMLviewClone.resizeSignal.connect(self.resizeHTML)

		# try to export figures
		# make_snapshot(snapshot, my_pyecharts.render(), "/Users/leil/Documents/Projects/VGenes/test.png")
		# self.ui.HTMLview.page().view().toPlainText(self.ui.HTMLview._callable)
		# print(self.ui.HTMLview.html)

		# build qweb channel
		channelClone = QWebChannel(self.ui.HTMLviewClone.page())
		my_objectClone = MyObjectCls(self.ui.HTMLviewClone)
		channelClone.registerObject('connection', my_objectClone)
		self.ui.HTMLviewClone.page().setWebChannel(channelClone)
		my_objectClone.downloadFigSignal.connect(self.downloadSVGClone)

	def downloadSVGClone(self, msg):
		options = QtWidgets.QFileDialog.Options()
		save_file_name, _ = QtWidgets.QFileDialog.getSaveFileName(self,
		                                                          "My_svg",
		                                                          "My_svg",
		                                                          "Scalable Vector Graphics (*.svg);;All Files (*)",
		                                                          options=options)

		try:
			file_handle = open(save_file_name, 'w')
			file_handle.write(msg)
			file_handle.close()
		except:
			return

	def GenerateFigure(self):
		global DBFilename
		#global data

		# select data or not
		if self.ui.checkBoxSelection.isChecked():
			where_statement = 'WHERE SeqName IN '
			selected_list = self.getTreeCheckedChild()
			selected_list = selected_list[3]
			if len(selected_list) == 0:
				Msg = 'No record has been selected! Will use the entire database!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				where_statement = 'WHERE 1'
			else:
				selected = "','".join(selected_list)
				where_statement = where_statement + "('" + selected + "')"
		else:
			where_statement = 'WHERE 1'

		# pie chart
		if self.ui.tabWidgetFig.currentIndex() == 0:
			# get data
			field = self.ui.comboBoxPie.currentText()
			if field == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			field = re.sub(r'\(.+', '', field)
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			if len(DataIn) == 0:
				Msg = 'No records can be fetched!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

			data = []
			for element in DataIn:
				data.append(element[0])
			result = Counter(data)
			labels = result.keys()
			values = result.values()

			my_pyecharts = (
				Pie(init_opts=opts.InitOpts(width= str(self.ui.HTMLview.w) + "px", height= str(self.ui.HTMLview.h) + "px", renderer='svg'))
				.add('', [list(z) for z in zip(labels, values)], radius=["40%", "75%"])
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show = self.ui.checkBoxFigLegend.isChecked()
					),
				)
				.set_series_opts(label_opts=opts.LabelOpts(formatter=" {b}: {c} ({d}%)"))
			)   #
		# Bar chart
		elif self.ui.tabWidgetFig.currentIndex() == 1:
			# get data
			field1 = self.ui.comboBoxCol1.currentText()
			field2 = self.ui.comboBoxCol2.currentText()
			if field1 == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			if field2 == field1:
				QMessageBox.warning(self, 'Warning', 'Please select different group factors for field1 and field2!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

			field1 = re.sub(r'\(.+', '', field1)
			field2 = re.sub(r'\(.+', '', field2)

			multi_factor = False
			if field2 == "":
				field = field1
			else:
				field = field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
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

				my_bar = Bar(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))\
					.add_xaxis(labels)\
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
						),
					)

				for group in dic_keys:
					cur_data = data[group]
					group_data = []
					for ele in labels:
						group_data.append(cur_data.count(ele))
					my_bar.add_yaxis(group, group_data,stack=stack)
				my_bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False, formatter=" {b}: {c}"))

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
					Bar(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))
						.add_xaxis(labels)
						.add_yaxis(field1, values)
						.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(
							is_show=self.ui.checkBoxFigLegend.isChecked()
						),
					)
					.set_series_opts(label_opts=opts.LabelOpts(is_show=False, formatter=" {b}: {c}"))
				)
		# Box plot
		elif self.ui.tabWidgetFig.currentIndex() == 2:
			# get data
			data_field = self.ui.comboBoxBoxData.currentText()
			field1 = self.ui.comboBoxBox1.currentText()
			field2 = self.ui.comboBoxBox2.currentText()

			data_field = re.sub(r'\(.+', '', data_field)
			field1 = re.sub(r'\(.+', '', field1)
			field2 = re.sub(r'\(.+', '', field2)

			if field1 == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			if field2 == field1:
				QMessageBox.warning(self, 'Warning', 'Please select different group factors for field1 and field2!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			multi_factor = False
			if field2 == "":
				field = data_field + "," + field1
			else:
				field = data_field + "," + field1 + "," + field2
				multi_factor = True
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			box_data = [i[0] for i in DataIn]
			try:
				box_data = list(map(float, box_data))
			except:
				QMessageBox.warning(self, 'Warning', 'The data field is not numerical! Check your input!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

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

				my_bar = Boxplot(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))\
					.add_xaxis(labels)\
					.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
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
					Boxplot(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))
						.add_xaxis(labels)
						.add_yaxis(field1, Boxplot.prepare_data(data_v1))
						.set_global_opts(
						title_opts=opts.TitleOpts(title=""),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						xaxis_opts=opts.AxisOpts(
							name=field1,
							name_location='center',
							name_gap=30,
						),
						yaxis_opts=opts.AxisOpts(
							name='Count',
							name_location='center',
							name_gap=30,
						),
					)
				)
		# Word Cloud
		elif self.ui.tabWidgetFig.currentIndex() == 3:
			# get data
			field = self.ui.comboBoxWord.currentText()

			field = re.sub(r'\(.+', '', field)

			if field == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
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
				WordCloud(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))
				.add("", word_data, word_size_range=[40, 200], shape=SymbolType.DIAMOND)
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(
						is_show=self.ui.checkBoxFigLegend.isChecked()
					),
				)
			)

			# load figure
			html_path = os.path.join(temp_folder,'figure.html')
			my_pyecharts.render(path=html_path)
			# adjust the window size seting
			file_handle = open(html_path, 'r')
			lines = file_handle.readlines()
			file_handle.close()
			# edit js line
			js_line = '<script type="text/javascript" src="' + \
			          os.path.join(js_folder,'echarts.js') + '"></script>' + \
			          '<script src="' + os.path.join(js_folder,'jquery.js') + '"></script>' + \
			          '<script src="qrc:///qtwebchannel/qwebchannel.js"></script>'
			lines[5] = js_line
			# edit style line
			style_line = lines[10]
			style_pos = style_line.find('style')
			style_line = style_line[
			             0:style_pos] + 'style="position: fixed; top: 0px; left: 5%;width:90%; height:' + str(
				self.ui.HTMLview.h - 20) + 'px;"></div>'
			lines[10] = style_line
			insert_js = '<script type="text/javascript">$(document).ready(function() {' \
			            'new QWebChannel(qt.webChannelTransport, function(channel) {' \
			            'var my_object = channel.objects.connection;$("#download").click(function(){' \
			            'my_object.download(text);});});});</script>'
			insert_btn = '<input id="download" type="button" value="" style="display:none;"/>'
			lines = lines[:6] + [insert_js] + lines[6:10] + [insert_btn] + lines[10:]
			content = '\n'.join(lines)
			file_handle = open(html_path, 'w')
			file_handle.write(content)
			file_handle.close()
			# show local HTML
			self.ui.HTMLview.load(QUrl('file://' + html_path))
			self.ui.HTMLview.show()
			self.ui.HTMLview.html = "loaded"
			self.ui.HTMLview.resizeSignal.connect(self.resizeHTML)

			# build qweb channel
			channel = QWebChannel(self.ui.HTMLview.page())
			my_object = MyObjectCls(self.ui.HTMLview)
			channel.registerObject('connection', my_object)
			self.ui.HTMLview.page().setWebChannel(channel)
			my_object.downloadFigSignal.connect(self.downloadSVG)
			my_object.updateSelectionSignal.connect(self.updateSelection)
			return
		# River chart
		elif self.ui.tabWidgetFig.currentIndex() == 4:
			# get data
			field1 = self.ui.comboBoxRiver1.currentText()
			field2 = self.ui.comboBoxRiver2.currentText()

			field1 = re.sub(r'\(.+', '', field1)
			field2 = re.sub(r'\(.+', '', field2)

			if field1 == "" or field2 == "":
				QMessageBox.warning(self, 'Warning', 'Your data field or group field is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

			field = field1 + "," + field2
			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
			DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			label_data = []
			time_data = []
			for element in DataIn:
				label_data.append(element[1])
				time_data.append(element[0])
			try:
				time_data = list(map(float, time_data))
			except:
				QMessageBox.warning(self, 'Warning', 'The data field is not numerical! Check your input!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

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
				ThemeRiver(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))
				.add(
					labels,
					data_river,
					singleaxis_opts=opts.SingleAxisOpts(type_='value', min_='dataMin', max_='dataMax'),
				)
				.set_global_opts(
					title_opts=opts.TitleOpts(title=""),
					legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
					xaxis_opts=opts.AxisOpts(
						name=field1,
						name_location='center',
						name_gap=30,
					),
					tooltip_opts=opts.TooltipOpts(trigger="axis", axis_pointer_type="line")
				)
			)
		# Tree Map
		elif self.ui.tabWidgetFig.currentIndex() == 5:
			group1 = self.ui.comboBoxTree1.currentText()
			group2 = self.ui.comboBoxTree2.currentText()
			group3 = self.ui.comboBoxTree3.currentText()

			group1 = re.sub(r'\(.+', '', group1)
			group2 = re.sub(r'\(.+', '', group2)
			group3 = re.sub(r'\(.+', '', group3)

			if group1 == "":
				QMessageBox.warning(self, 'Warning', 'Your Field1 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return
			if group2 == "":
				data = []

				field = group1
				SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				x_data = [d[0] for d in DataIn]

				x_data_res = Counter(x_data)
				for ele in x_data_res:
					unit = {"value": x_data_res[ele], "name": ele}
					data.append(unit)
			else:
				if group3 == "":
					if group2 == group1:
						QMessageBox.warning(self, 'Warning',
						                    'Please select different group factors for field1 and field2!',
						                    QMessageBox.Ok, QMessageBox.Ok)
						return
					data = []

					field = group1 + ',' + group2
					SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
					DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
					x_data = [d[0] for d in DataIn]
					y_data = [d[1] for d in DataIn]

					x_data_res = Counter(x_data)
					x_labels = x_data_res.keys()

					for x_label in x_labels:
						cur_x_label_data = []
						for i in range(0, len(x_data)):
							if x_data[i] == x_label:
								cur_x_label_data.append(y_data[i])
						cur_data = []
						cur_x_data_res = Counter(cur_x_label_data)
						for ele in cur_x_data_res:
							sub_unit = {"value": cur_x_data_res[ele], "name": ele}
							cur_data.append(sub_unit)

						unit = {"value": x_data_res[x_label], "name": x_label, "children": cur_data}
						data.append(unit)

				else:
					if group2 == group1 or group1 == group3 or group2 == group3:
						QMessageBox.warning(self, 'Warning',
						                    'Please select different group factors for field1, field2, and field3!',
						                    QMessageBox.Ok, QMessageBox.Ok)
						return
					data = []

					field = group1 + ',' + group2 + ',' + group3
					SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
					DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
					x_data = [d[0] for d in DataIn]
					y_data = [d[1] for d in DataIn]
					z_data = [d[2] for d in DataIn]

					x_data_res = Counter(x_data)
					x_labels = x_data_res.keys()

					# level 1
					for x_label in x_labels:
						cur_x_label_data = []
						cur_z_data = []
						for i in range(0, len(x_data)):
							if x_data[i] == x_label:
								cur_x_label_data.append(y_data[i])
								cur_z_data.append(z_data[i])
						cur_data = []
						cur_x_data_res = Counter(cur_x_label_data)
						y_labels = cur_x_data_res.keys()

						# level 2
						for y_label in y_labels:
							cur_y_label_data = []
							for i in range(0, len(cur_x_label_data)):
								if cur_x_label_data[i] == y_label:
									cur_y_label_data.append(cur_z_data[i])
							cur_sub_data = []
							cur_y_data_res = Counter(cur_y_label_data)
							for ele in cur_y_data_res:
								sub_sub_unit = {"value": cur_y_data_res[ele], "name": ele}
								cur_sub_data.append(sub_sub_unit)

							sub_unit = {"value": cur_x_data_res[y_label], "name": y_label, "children": cur_sub_data}
							cur_data.append(sub_unit)

						unit = {"value": x_data_res[x_label], "name": x_label, "children": cur_data}
						data.append(unit)

			if self.ui.radioButtonTree.isChecked():
				tree_data = [{"name":"MyData", "children":data}]
				my_pyecharts = (
					Tree(init_opts=opts.InitOpts(width="380px", height="380px", renderer='canvas'))
						.add("MyData", tree_data)
						.set_global_opts(
						title_opts=opts.TitleOpts(title="Tree"),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						#toolbox_opts=opts.ToolboxOpts()
					)
				)
			else:
				my_pyecharts = (
					TreeMap(init_opts=opts.InitOpts(width="380px", height="380px", renderer='canvas'))
					.add(
						series_name="MyData",
						data=data,
						visual_min=300,
						leaf_depth=1,
						label_opts=opts.LabelOpts(position="inside"),
					)
					.set_global_opts(
						title_opts=opts.TitleOpts(title="TreeMap"),
						legend_opts=opts.LegendOpts(is_show=self.ui.checkBoxFigLegend.isChecked()),
						#toolbox_opts=opts.ToolboxOpts()
					)
				)
		# Scatter Chart
		elif self.ui.tabWidgetFig.currentIndex() == 6:
			dim1 = self.ui.comboBoxScatterX.currentText()
			dim2 = self.ui.comboBoxScatterY.currentText()
			group = self.ui.comboBoxScatterGroup.currentText()

			dim1 = re.sub(r'\(.+', '', dim1)
			dim2 = re.sub(r'\(.+', '', dim2)
			group = re.sub(r'\(.+', '', group)

			if dim1 == "" or dim2 == "":
				QMessageBox.warning(self, 'Warning', 'Your dim1 or dim2 is empty!',
				                    QMessageBox.Ok, QMessageBox.Ok)
				return

			# create figure
			my_scatter = Scatter(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))\
				.set_series_opts(label_opts=opts.LabelOpts(is_show=False))\
				.set_global_opts(
					xaxis_opts=opts.AxisOpts(
						type_="value",
						splitline_opts=opts.SplitLineOpts(is_show=False),
						name=dim1,
						name_location='center',
						name_gap=30,
					),
					yaxis_opts=opts.AxisOpts(
						type_="value",
						name=dim2,
						name_location='center',
						name_gap=30,
						axistick_opts=opts.AxisTickOpts(is_show=True),
						splitline_opts=opts.SplitLineOpts(is_show=False),
					),
					tooltip_opts=opts.TooltipOpts(is_show=True, formatter="{c}, {a}"),
					legend_opts=opts.LegendOpts(
						is_show=self.ui.checkBoxFigLegend.isChecked()
					),
				)

			# load data
			if group == "":
				field = dim1 + "," + dim2
				SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				x_data = []
				y_data = []
				err = False
				for d in DataIn:
					try:
						x = float(d[0])
						y = float(d[1])
						x_data.append(x)
						y_data.append(y)
					except:
						err = True
				if err == True:
					pass
					#QMessageBox.information(self, 'Information', 'Non-numerical values/records have been removed from this figure!',
					#                    QMessageBox.Ok, QMessageBox.Ok)

				'''
				x_data = [d[0] for d in DataIn]
				y_data = [d[1] for d in DataIn]
				
				try:
					x_data = list(map(float,x_data))
					y_data = list(map(float, y_data))
				except:
					QMessageBox.warning(self, 'Warning', 'The dim1 or dim2 field is not numerical! Check your input!',
					                    QMessageBox.Ok, QMessageBox.Ok)
					return
				'''

				# attach data
				my_scatter.add_xaxis(xaxis_data=x_data)
				my_scatter.add_yaxis(series_name="Data", y_axis=y_data, label_opts=opts.LabelOpts(is_show=False))

			else:
				field = dim1 + "," + dim2 + "," + group
				SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

				x_data = []
				y_data = []
				group_data = []
				err = False
				for d in DataIn:
					try:
						x = float(d[0])
						y = float(d[1])
						x_data.append(x)
						y_data.append(y)
						group_data.append(d[2])
					except:
						err = True
				if err == True:
					pass
					#QMessageBox.information(self, 'Information',
					#                    'Non-numerical values/records have been removed from this figure!',
					#                    QMessageBox.Ok, QMessageBox.Ok)

				'''
				x_data = [d[0] for d in DataIn]
				y_data = [d[1] for d in DataIn]
				group_data = [d[2] for d in DataIn]
				try:
					x_data = list(map(float, x_data))
					y_data = list(map(float, y_data))
				except:
					QMessageBox.warning(self, 'Warning', 'The dim1 or dim2 field is not numerical! Check your input!',
					                    QMessageBox.Ok, QMessageBox.Ok)
					return
				'''
				group_result = Counter(group_data)
				groups = list(group_result.keys())

				for group in groups:
					sub_x_data = []
					sub_y_data = []
					for i in range(0, len(group_data)):
						if group_data[i] == group:
							sub_x_data.append(x_data[i])
							sub_y_data.append(y_data[i])

					# attach data
					my_scatter.add_xaxis(xaxis_data=sub_x_data)
					my_scatter.add_yaxis(series_name=group, y_axis=sub_y_data,label_opts=opts.LabelOpts(is_show=False))

			my_pyecharts = (
				my_scatter
			)
		# Heatmap
		elif self.ui.tabWidgetFig.currentIndex() == 7:
			if len(self.HeatmapList) == 0:
				Msg = 'Please select at least one field!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

			field_list = [re.sub(r'\(.+', '', item) for item in self.HeatmapList]

			field = ",".join(field_list)
			if self.ui.comboBoxSortField.currentText() == '':
				sort_statement = ''
				field = field + ',SeqName'
			else:
				sort_field = re.sub(r'\(.+', '', self.ui.comboBoxSortField.currentText())
				sort_statement = ' ORDER BY ' + sort_field + ' ASC'
				field = field + ',' + sort_field

			SQLStatement = 'SELECT ' + field + ' FROM vgenesDB ' + where_statement + sort_statement
			try:
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			except:
				Msg = 'SQL error! Current SQL statemernt is:\n' + SQLStatement
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

			# parse out record with all null value
			if self.ui.checkBoxHideNull.isChecked():
				index_remove = []
				index = 0
				for record in DataIn:
					values = record[:len(record)-1]
					not_null = False
					for value in values:
						try:
							int(value)
							not_null = True
							break
						except:
							pass

					if not_null == False:
						index_remove.append(index)

					index += 1

				# delete all null record
				cnt = 0
				for index in index_remove:
					del DataIn[index - cnt]
					cnt += 1

			num_row = len(DataIn)
			num_col = len(self.HeatmapList)

			if num_row == 0:
				Msg = 'No record fetched!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return
			'''
			try:
				value = [[i, j, int(DataIn[i][j])] for i in range(num_row) for j in range(num_col)]
				xaxis_data = [DataIn[i][-1] for i in range(num_row)]
			except:
				Msg = 'Some values of your selected field/record are not numbers!'
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return
				
			min_value = min(row[2] for row in value)
			max_value = max(row[2] for row in value)
			'''
			# replace all missing values by null
			value = []
			min_value = 9999
			max_value = 0
			for i in range(num_row):
				for j in range(num_col):
					try:
						value.append([i, j, int(DataIn[i][j])])
						if int(DataIn[i][j]) < min_value:
							min_value = int(DataIn[i][j])
						if int(DataIn[i][j]) > max_value:
							max_value = int(DataIn[i][j])
					except:
						value.append([i, j, 'null'])
			xaxis_data = [DataIn[i][-1] for i in range(num_row)]

			if self.ui.radioButtonLog10.isChecked():
				if min_value >= 0:
					for i in range(len(value)):
						try:
							value[i][2] = numpy.log1p(value[i][2])/numpy.log(10)
						except:
							pass

					min_value = numpy.log1p(min_value)/numpy.log(10)
					max_value = numpy.log1p(max_value) / numpy.log(10)
				else:
					Msg = 'Negative values detected in your data! Log scale is not appliable!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					return
			elif self.ui.radioButtonLog2.isChecked():
				if min_value >= 0:
					for i in range(len(value)):
						try:
							value[i][2] = numpy.log1p(value[i][2])/numpy.log(2)
						except:
							pass

					min_value = numpy.log1p(min_value) / numpy.log(2)
					max_value = numpy.log1p(max_value) / numpy.log(2)
				else:
					Msg = 'Negative values detected in your data! Log scale is not appliable!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					return

			my_pyecharts = (
				HeatMap(init_opts=opts.InitOpts(width="380px", height="380px", renderer='svg'))
					.add_xaxis(xaxis_data)
					.add_yaxis(
					"My data selection",
					self.HeatmapList,
					value,
					label_opts=opts.LabelOpts(is_show=False, position="inside"),
				)
					.set_global_opts(
					title_opts=opts.TitleOpts(title="HeatMap"),
					visualmap_opts=opts.VisualMapOpts(min_=min_value, max_=max_value),
				)
			)
		# load figure
		html_path = os.path.join(temp_folder,'figure.html')
		my_pyecharts.render(path=html_path)
		# adjust the window size seting
		file_handle = open(html_path, 'r')
		lines = file_handle.readlines()
		file_handle.close()
		# edit js line
		js_line = '<script type="text/javascript" src="' + \
		          os.path.join(js_folder, 'echarts.js') + '"></script>' + \
		          '<script src="' + os.path.join(js_folder, 'jquery.js') + '"></script>' + \
		          '<script src="qrc:///qtwebchannel/qwebchannel.js"></script>'
		lines[5] = js_line
		# edit style line
		style_line = lines[9]
		style_pos = style_line.find('style')
		style_line = style_line[
		             0:style_pos] + 'style="position: fixed; top: 0px; left: 5%;width:90%; height:' + str(
			self.ui.HTMLview.h - 20) + 'px;"></div>'
		lines[9] = style_line
		insert_js = '<script type="text/javascript">$(document).ready(function() {' \
		            'new QWebChannel(qt.webChannelTransport, function(channel) {' \
		            'var my_object = channel.objects.connection;$("#download").click(function(){' \
		            'my_object.download(text);});$("#update").click(function(){' \
		            'my_object.updateSelection(text);});});});</script>'
		insert_btn = '<input id="download" type="button" value="" style="display:none;"/>' \
		             '<input id="update" type="button" value="" style="display:none;"/>'
		lines = lines[:6] + [insert_js] + lines[6:9] + [insert_btn] + lines[9:]

		if self.ui.tabWidgetFig.currentIndex() in [0,1,2,6]:
			# insert click response function
			echart_init_line = lines[13]
			matchObj = re.match(r'.+var\s(\S+)\s=', echart_init_line)
			chart_id = matchObj.group(1)
			js_cmd = chart_id + ".on('click', function (params) {" \
			                    "if(params.data['0'] == null){text = params.name + ',' + params.seriesName + ',0,0';}" \
			                    "else{text = params.name + ',' + params.seriesName + ','+params.data['0']+','+params.data['1'];}" \
			                    "$('#update').click();});"
			lines = lines[:-3] + [js_cmd] + lines[-3:]

		content = '\n'.join(lines)
		file_handle = open(html_path, 'w')
		file_handle.write(content)
		file_handle.close()
		# show local HTML
		self.ui.HTMLview.load(QUrl('file://' + html_path))
		self.ui.HTMLview.show()
		self.ui.HTMLview.html = "loaded"
		self.ui.HTMLview.resizeSignal.connect(self.resizeHTML)

		# try to export figures
		#make_snapshot(snapshot, my_pyecharts.render(), "/Users/leil/Documents/Projects/VGenes/test.png")
		#self.ui.HTMLview.page().view().toPlainText(self.ui.HTMLview._callable)
		#print(self.ui.HTMLview.html)

		# build qweb channel
		channel = QWebChannel(self.ui.HTMLview.page())
		my_object = MyObjectCls(self.ui.HTMLview)
		channel.registerObject('connection', my_object)
		self.ui.HTMLview.page().setWebChannel(channel)
		my_object.downloadFigSignal.connect(self.downloadSVG)
		my_object.updateSelectionSignal.connect(self.updateSelection)

		# this section will add information for latest figure. The information will be used to update element selections on the left
		if self.ui.tabWidgetFig.currentIndex() == 0:
			self.ui.HTMLview.info = ['Pie', re.sub(r'\(.+', '', self.ui.comboBoxPie.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 1:
			self.ui.HTMLview.info = ['Bar', re.sub(r'\(.+', '', self.ui.comboBoxCol1.currentText()), re.sub(r'\(.+', '', self.ui.comboBoxCol2.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 2:
			self.ui.HTMLview.info = ['Box', re.sub(r'\(.+', '', self.ui.comboBoxBoxData.currentText()), re.sub(r'\(.+', '',self.ui.comboBoxBox1.currentText()), re.sub(r'\(.+', '', self.ui.comboBoxBox2.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 3:
			self.ui.HTMLview.info = ['Word', re.sub(r'\(.+', '', self.ui.comboBoxWord.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 4:
			self.ui.HTMLview.info = ['River', re.sub(r'\(.+', '', self.ui.comboBoxRiver1.currentText()), re.sub(r'\(.+', '', self.ui.comboBoxRiver2.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 5:
			self.ui.HTMLview.info = ['Tree', re.sub(r'\(.+', '', self.ui.comboBoxTree1.currentText()), re.sub(r'\(.+', '', self.ui.comboBoxTree2.currentText()) , re.sub(r'\(.+', '', self.ui.comboBoxTree3.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 6:
			self.ui.HTMLview.info = ['Scatter', re.sub(r'\(.+', '', self.ui.comboBoxScatterX.currentText()), re.sub(r'\(.+', '', self.ui.comboBoxScatterY.currentText()) ,re.sub(r'\(.+', '', self.ui.comboBoxScatterGroup.currentText())]
		elif self.ui.tabWidgetFig.currentIndex() == 7:
			pass

	def resizeHTML(self):
		if self.ui.HTMLview.html == '':
			return
		else:
			self.GenerateFigure()

	def resizeHTMLClone(self):
		if self.ui.HTMLviewClone.html == '':
			return
		else:
			self.GenerateFigureClone()

	def updateSelection(self, msg):
		if self.ui.checkBoxUpdateSelection.isChecked():
			# split msg
			messages = msg.split(',')
			# fetch data record
			info = self.ui.HTMLview.info
			if info[0] == "Pie":
				field = info[1]
				where_statement = 'WHERE ' + field + " = '" + messages[0] + "'"
				SQLStatement = 'SELECT SeqName FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			elif info[0] == "Bar":
				field1 = info[1]
				field2 = info[2]
				if field2 == "":
					where_statement = 'WHERE ' + field1 + " = '" + messages[0] + "'"
				else:
					where_statement = 'WHERE ' + field1 + " = '" + messages[0] + "'" + ' AND ' + field2 + " = '" + messages[1] + "'"
				SQLStatement = 'SELECT SeqName FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			elif info[0] == "Box":
				field1 = info[2]
				field2 = info[3]
				if field2 == "":
					where_statement = 'WHERE ' + field1 + " = '" + messages[0] + "'"
				else:
					where_statement = 'WHERE ' + field1 + " = '" + messages[0] + "'" + ' AND ' + field2 + " = '" + messages[1] + "'"
				SQLStatement = 'SELECT SeqName FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			elif info[0] == "Scatter":
				fieldx = info[1]
				fieldy = info[2]
				field = info[3]
				if field == "":
					where_statement = 'WHERE ' + fieldx + " = '" + messages[2] + "'" + \
					                  ' AND ' + fieldy + " = '" + messages[3] + "'"
				else:
					where_statement = 'WHERE ' + fieldx + " = '" + messages[2] + "'" + \
					                  ' AND ' + fieldy + " = '" + messages[3] + "'" + \
					                  ' AND ' + field + " = '" + messages[1] + "'"
				SQLStatement = 'SELECT SeqName FROM vgenesDB ' + where_statement
				DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
			else:
				return

			# update the selection
			self.clearTreeChecks()
			NumFound = len(DataIn)
			i = 0
			for item in DataIn:
				Seqname = item[0]
				found = self.ui.treeWidget.findItems(Seqname, Qt.MatchRecursive, 0)
				i += 1
				for record in found:
					if i == NumFound - 1:
						wasClicked = True
					record.setCheckState(0, Qt.Checked)

			self.match_tree_to_table()

	def downloadSVG(self, msg):
		options = QtWidgets.QFileDialog.Options()
		save_file_name, _ = QtWidgets.QFileDialog.getSaveFileName(self,
		                                                      "My_svg",
		                                                      "My_svg",
		                                                      "Scalable Vector Graphics (*.svg);;All Files (*)",
		                                                      options=options)

		try:
			file_handle = open(save_file_name, 'w')
			file_handle.write(msg)
			file_handle.close()
		except:
			return

	def downloadFig(self):
		js_cmd= 'text=document.getElementsByTagName("svg")[0].parentNode.innerHTML;$("#download").click();'

		#js_cmd= 'svg=document.getElementsByTagName("svg")[0];var a = document.createElement("a");a.href = svg.src;' \
		#        'a.download = "~/Downloads/test.svg";a.click();'
		self.ui.HTMLview.page().runJavaScript(js_cmd)

	def downloadFigClone(self):
		js_cmd= 'text=document.getElementsByTagName("svg")[0].parentNode.innerHTML;$("#download").click();'

		#js_cmd= 'svg=document.getElementsByTagName("svg")[0];var a = document.createElement("a");a.href = svg.src;' \
		#        'a.download = "~/Downloads/test.svg";a.click();'
		self.ui.HTMLviewClone.page().runJavaScript(js_cmd)

	def _callable(self, data):
		self.html = data

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
			#project = self.ui.cboTreeOp1.currentText()
			#fieldname = self.TransLateFieldtoReal(project, True)
			fieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
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
			#group = self.ui.cboTreeOp2.currentText()
			#fieldname = self.TransLateFieldtoReal(group, True)
			#project = self.ui.cboTreeOp1.currentText()
			#Projfieldname = self.TransLateFieldtoReal(project, True)
			fieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
			Projfieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())

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
			#Subgroup = self.ui.cboTreeOp3.currentText()
			#fieldname = self.TransLateFieldtoReal(Subgroup, True)
			#group = self.ui.cboTreeOp2.currentText()
			#Groupfieldname = self.TransLateFieldtoReal(group, True)
			#project = self.ui.cboTreeOp1.currentText()
			#Projfieldname = self.TransLateFieldtoReal(project, True)
			fieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
			Groupfieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
			Projfieldname = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())

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
		# clear tree
		self.ui.treeWidget.clear()
		# add all children
		self.TreeAddItems(self.ui.treeWidget.invisibleRootItem(), SQLFields)
		# rebuild name index
		self.buildNameidnex(SQLFields)

		global wasClicked
		wasClicked = True

	def buildNameidnex(self,SQLFields):
		global NameIndex

		field1, field2, field3 = SQLFields
		NameIndex.clear()
		if field1 != 'None' and field2 != 'None' and field3 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', ' + field3 + ', SeqName'
		elif field1 != 'None' and field2 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', SeqName'
		elif field1 != 'None' and field3 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field1 + ', ' + field3 + ', SeqName'
		elif field2 != 'None' and field3 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field2 + ', ' + field3 + ', SeqName'
		elif field1 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field1 + ', SeqName'
		elif field2 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field2 + ', SeqName'
		elif field3 != 'None':
			SQLStatement = 'select SeqName from vgenesdb ORDER BY ' + field3 + ', SeqName'
		else:
			SQLStatement = 'select SeqName from vgenesdb ORDER BY SeqName'

		DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		Maxi = len(DataIs)
		for i in range(0, Maxi):
			NameIs = DataIs[i][0]
			NameIndex[NameIs] = i

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
		else:
			if item.checkState(column) == Qt.Checked:
				self.CheckBelow(item, True)
				self.CheckUp(item, True)
				# childs = item.childCount()
			if item.checkState(column) == Qt.Unchecked:
				self.CheckBelow(item, False)
				self.CheckUp(item, False)
				print("unchecked")

		self.match_tree_to_table()

	def tree_to_table_selection(self):
		Selected = self.ui.treeWidget.selectedItems()
		try:
			Selected = Selected[-1]
		except:
			Selected = Selected[0]
		name = Selected.text(0)

		print(name)

		rows = self.ui.SeqTable.rowCount()
		for row in range(rows):
			cur_name = self.ui.SeqTable.item(row, 1).text()
			if cur_name == name:
				self.ui.SeqTable.setCurrentCell(row,1)
				return

	def table_to_tree_selection(self):
		global from_table
		try:
			items = self.ui.SeqTable.selectedItems()
			item = items[-1]
			name = self.ui.SeqTable.item(self.ui.SeqTable.indexFromItem(item).row(), 1).text()

			print(name)

			found = self.ui.treeWidget.findItems(name, Qt.MatchRecursive, 0)
			if len(found) > 0:
				found = found[0]
				from_table = True
				self.ui.treeWidget.setCurrentItem(found)
				from_table = False
		except:
			return

	def select_tree_by_name(self, name):
		found = self.ui.treeWidget.findItems(name, Qt.MatchRecursive, 0)
		if len(found) > 0:
			found = found[0]
			self.ui.treeWidget.setCurrentItem(found)

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

	def getTreeCheckedChild(self):
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
								for SGindex in range(SubGroup.childCount()):
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
		elif self.ui.tabWidget.currentIndex() == 2:
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
		else:
			Msg = 'Current page has nothing to print!\n ' \
			      'Print functions ony available for Tablulated Data, Sequences and Alignment page!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

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
	def on_actionDuplicate_current_triggered(self):
		SQLStatement = 'SELECT * FROM vgenesDB WHERE SeqName = "' + data[0] + '"'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		copy_data = list(DataIn[0])
		new_name = copy_data[0] + '_copy'
		copy_data[0] = new_name

		SQLStatement = 'SELECT ID FROM vgenesDB'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		id_list = [int(i[0]) for i in DataIn]
		copy_data[119] = str(max(id_list) + 1)
		my_list = ['?']*len(FieldList)

		CMD = 'INSERT INTO vgenesDB(' + ','.join(FieldList) + ') VALUES(' + ','.join(my_list) + ')'
		conn = db.connect(DBFilename)
		cursor = conn.cursor()
		try:
			cursor.execute(CMD, copy_data)
			conn.commit()  # saves data into file

			self.refreshDB()
			self.on_btnUpdateTree_clicked()

			found = self.ui.treeWidget.findItems(new_name, Qt.MatchRecursive, 0)
			if len(found) > 0:
				found = found[0]
				from_table = True
				self.ui.treeWidget.setCurrentItem(found)
				from_table = False
				self.tree_to_table_selection()
		except:
			QMessageBox.warning(self, 'Warning', 'Failed to make a copy of current record!\nMost likely you already have a copy of current record!\nTo make a new copy, please rename the old copy first!',
			                        QMessageBox.Ok, QMessageBox.Ok)

		conn.close()
		
		# focus on the copy

	@pyqtSlot()
	def on_actionDelete_record_triggered(self):
		fields = ['SeqName', 'ID']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])

		sel_names = self.getTreeCheckedChild()
		sel_names = sel_names[3]

		if len(sel_names) == 0:
			question = 'You did not check any sequences, do you want delete current selected sequence?\n' + data[0]
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)

			if answer == 'Yes':
				SQLStatement = ' FROM vgenesDB WHERE SeqName = "' + data[0] + '"'
				VGenesSQL.deleterecords(DBFilename, SQLStatement)

				QMessageBox.information(self, 'Information', 'The selected record has been deleted!',
				                        QMessageBox.Ok, QMessageBox.Ok)
				if DBFilename != None:
					if os.path.isfile(DBFilename):
						self.LoadDB(DBFilename)
		else:
			question = 'Do you want delete current selected sequences?\n' + '\n'.join(sel_names)
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)

			if answer == 'Yes':
				SQLStatement = ' FROM vgenesDB WHERE SeqName IN ("' + '","'.join(sel_names) + '")'
				VGenesSQL.deleterecords(DBFilename, SQLStatement)

				QMessageBox.information(self, 'Information', 'The selected record has been deleted!',
				                        QMessageBox.Ok, QMessageBox.Ok)
				if DBFilename != None:
					if os.path.isfile(DBFilename):
						self.LoadDB(DBFilename)
						self.on_btnUpdateTree_clicked()
						


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

	@pyqtSlot()
	def on_actionFind_Clonal_triggered(self):
		global updateMarker
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
			#field = self.ui.cboFindField.currentText()
			#fieldsearch = self.TransLateFieldtoReal(field, True)
			fieldsearch = re.sub(r'\(.+', '', self.ui.cboFindField.currentText())
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

		#model = self.ui.tableView.model()
		#model.refresh()

		if answer == 'Yes':
			if DBFilename != None:
				if os.path.isfile(DBFilename):
					self.LoadDB(DBFilename)
		else:
			if re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText()) == 'Clonal Pool' \
					or re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText()) == 'Clonal Pool' \
					or re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText()) == 'Clonal Pool':
				self.on_btnUpdateTree_clicked()

		self.findTreeItem(Currentrecord)
		ErLog2 = str(CPs) + ' clonal pools containing ' + str(CPseqs) + ' sequences were identified from ' + str(
			TotSeqs) + ' total sequences analyzed.\n'

		SQLStatement = 'SELECT ClonalPool FROM vgenesDB'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		list1 = []
		for ele in DataIn:
			list1.append(ele[0])
		list_unique = list(set(list1))
		try:
			list_unique.remove('0')
		except:
			pass

		list_unique = [int(i) for i in list_unique]
		list_unique.sort()
		list_unique = ['Clone' + str(i) for i in list_unique]

		self.ui.comboBoxTree.clear()
		self.ui.comboBoxTree.addItems(list_unique)

		if remove == True:

			self.LoadDB(DBFilename)
			self.ui.txtFieldSearch.setText('Duplicate')
			self.ui.cboFindField.setCurrentText('Quality')
			done = self.on_btnFieldSearch_clicked()
			done = self.on_actionDelete_record_triggered()
			self.on_btnUpdateTree_clicked()

		if len(ErLog2) > 0:
			Erlog2 = ErLog2 + 'The following ' + str(
				Errs) + ' sequences could not be anaylzed for\nclonality because no CDR3s are indicated:\n' + ErLog
			ErlogFile = os.path.join(temp_folder, 'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs

			with open(ErlogFile, 'w') as currentFile:
				currentFile.write(Erlog2)

			self.ShowVGenesText(ErlogFile)

		if self.ui.tabWidget.currentIndex() == 0:
			self.load_table()
			self.match_tree_to_table()
			self.tree_to_table_selection()
		else:
			updateMarker == True

		self.initial_Clone()



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
			#field = self.ui.cboFindField.currentText()
			#fieldsearch = self.TransLateFieldtoReal(field, True)
			fieldsearch = re.sub(r'\(.+', '', self.ui.cboFindField.currentText())
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
		if f == '' or f == None:
			return

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

			if filename == None:
				return

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

	@pyqtSlot()
	def on_actionRename10x_triggered(self):
		global NameIndex

		from operator import itemgetter  #from operator import itemgetter   #		SeqList.sort(key=itemgetter(0, 1, 2, 3))

		QueryIS = 'Enter text to serve as the base name (i.e., subject number or name)'
		DefaultText = data[75]  # data[77] + '-Expressed'
		BaseName = setText(self, QueryIS, DefaultText)
		if BaseName == "Cancelled Action":
			return

		filename = os.path.join(working_prefix, '10x_barcodes.csv')
		with open(filename, 'r') as currentfile:
			# myCSVfile = csv.reader(currentfile)
			myCSVfile  = []
			# myCSVfile = currentfile.split(',')
			for row in currentfile:
				entry = row.strip('\n')
				myCSVfile.append(entry)

		# if not check any sequence, will import for all
		checked_list = self.getTreeCheckedChild()
		checked_list = checked_list[3]

		if len(checked_list) == 0:
			root = self.ui.treeWidget.invisibleRootItem()
			self.CheckBelow(root, True)

		#filtered_contig_annotations.csv
		answer = informationMessage(self,
		                            'Select annotation file from same folder as consensus.fasta called "outs/filtered_contig_annotations.csv"',
		                            'OK')
		typeOpen = 'csv'
		filename =openFile(self, typeOpen)
		if filename == None:
			return
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
			if filename != None:
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
			else:
				answerC == 'No'

		fields = ['SeqName', 'Id', 'GeneType']

		# checkedProjects, checkedGroups, checkedSubGroups, checkedkids = getTreeChecked()
		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, data[0])
		a = data

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
				#model = self.ui.tableView.model()
				#model.refresh()

		# update tree
		#SQLFields = (self.ui.cboTreeOp1.currentText(), self.ui.cboTreeOp2.currentText(),self.ui.cboTreeOp3.currentText())
		'''
		index1 = RealNameList.index(self.ui.cboTreeOp1.currentText())
		index2 = RealNameList.index(self.ui.cboTreeOp2.currentText())
		index3 = RealNameList.index(self.ui.cboTreeOp3.currentText())

		SQLFields = (FieldList[index1], FieldList[index2], FieldList[index3])
		'''
		SQLFields = (
			re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText()),
			re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText()),
			re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
		)

		self.initializeTreeView(SQLFields)
		self.ui.treeWidget.expandAll()

		# update table
		self.load_table()
		#answer = informationMessage(self, 'Close and restart database to see changes', 'OK')

	

	@pyqtSlot()
	def on_actionclearTrash_triggered(self):
		question = 'Clean all TEMP files?'
		buttons = 'YN'
		answer = questionMessage(self, question, buttons)
		if answer == 'No':
			return
		else:
			cmd = 'cd ' + temp_folder + '; rm -rf ' + temp_folder + '/*'
			try:
				os.system(cmd)
			except:
				QMessageBox.warning(self, 'Warning', 'Fail to clear TEMP folder!', QMessageBox.Ok,
				                    QMessageBox.Ok)
				return

			ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
			ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
			r = open(ErlogFile,'w')
			r.write('')
			r.close()
			r = open(ErlogFile2, 'w')
			r.write('')
			r.close()



	@pyqtSlot()
	def on_actionMergeMySeq_triggered(self):
		import shutil

		# WorkDir = '/Users/PCW-MacBookProRet/Applications/VGenes/FLASH-1.2.11/'

		try:
			filename = openfastq(self)

			read_1 = filename[0]
			read_2 = filename[1]

			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11', 'reads_1.fq')
			shutil.copy(read_1, WorkDir)
			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11', 'reads_2.fq')
			shutil.copy(read_2, WorkDir)

			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11')
			# (dirname, filename) = os.path.split(DBpathname)
			os.chdir(WorkDir)

			CommandLine = "./flash reads_1.fq reads_2.fq 2>&1 | tee flash.log"
			Result = os.popen(CommandLine)
			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11', 'out.extendedFrags.fastq')

			filename = saveFile(self, 'fastq')
			shutil.copy(WorkDir, filename)
			filename2 = filename + 'Unmerged-1.fastq'
			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11',
			                       'out.notCombined_1.fastq')
			shutil.copy(WorkDir, filename2)
			filename2 = filename + 'Unmerged-2.fastq'
			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11',
			                       'out.notCombined_2.fastq')
			shutil.copy(WorkDir, filename2)


			WorkDir = os.path.join(working_prefix, 'FLASH-1.2.11', 'flash.log')
			self.ShowVGenesText(WorkDir)
		except:
			return

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

	def getTreeSelected(self):
		root = self.ui.treeWidget.invisibleRootItem()
		ListSelected = []
		for item in self.ui.treeWidget.selectedItems():
			ListSelected.append(str(item.text(0)))
			# (item.parent() or root).removeChild(item)

		return ListSelected

	@pyqtSlot()
	def TreeSelectChanged(self):
		global from_table

		value = self.ui.treeWidget.selectedItems()

		name = ''
		# name2 = ''
		for item in value:
			name = item.text(0)
		#FieldCheck = self.FieldChangeCheck()
		try:
			MatchingIndex = NameIndex[name]
			self.DialScroll(MatchingIndex, False)
		except:
			print('wrong')
			return
			# self.findTableViewRecord(name)
			# print(role)
		SelectedItems =  self.getTreeSelected()
		NumSelected  = len(SelectedItems)
		if NumSelected > 1:
			NewHead  = str(NumSelected) + ' items selected'
			self.ui.label_Name.setText(NewHead)

		if from_table:
			return
		else:
			self.tree_to_table_selection()
			self.match_tree_to_table()

	def MatchingValue(self, IndexIs):
		try:
			return list(NameIndex.keys())[list(NameIndex.values()).index(int(IndexIs))]
		except:
			return 'None'

	def TreeviewOptions(self):
		global DontFindTwice
		if DontFindTwice == True:
			return
		try:
			sender = self.sender()
			if sender.id == 'TreeOp1':
				Option1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
				self.ui.cboTreeOp1.setToolTip('Press update tree to implement changes')
				if Option1 == 'None':
					DontFindTwice = True
					self.ui.cboTreeOp2.setCurrentText('None')
					self.ui.cboTreeOp3.setCurrentText('None')
					DontFindTwice = False
			elif sender.id == 'TreeOp2':
				Option1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
				self.ui.cboTreeOp1.setToolTip('Press update tree to implement changes')
				if Option1 == 'None':
					DontFindTwice = True
					self.ui.cboTreeOp2.setCurrentText('None')
					self.ui.cboTreeOp3.setCurrentText('None')
					DontFindTwice = False
					Msg = 'Upper level grouping factor is empty!\nPlease determine grouping factors from the upper level!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					return
				self.ui.cboTreeOp2.setToolTip('Press update tree to implement changes')
				Option2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
				if Option2 == 'None':
					DontFindTwice = True
					self.ui.cboTreeOp3.setCurrentText('None')
					DontFindTwice = False
			else:
				Option2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
				self.ui.cboTreeOp2.setToolTip('Press update tree to implement changes')
				Option3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
				if Option3 == 'None':
					pass
				elif Option2 == 'None':
					DontFindTwice = True
					self.ui.cboTreeOp3.setCurrentText('None')
					DontFindTwice = False
					Msg = 'Upper level grouping factor is empty!\nPlease determine grouping factors from the upper level!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					return
				self.ui.cboTreeOp3.setToolTip('Press update tree to implement changes')
		except:
			pass

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

		a = DBFilename
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
			pass
			#self.hide()
			#self.ApplicationStarted()
			# self.EditableSqlModel.refresh()

	def UpdateRecentList(self, DBFilename, AddOne):
		# todo need to make this filename fall in VGenes directory upon deployment
		# todo may need to switch this to configparser which is python modle to save ini files
		# todo change to app folder
		try:
			filename = os.path.join(working_prefix, 'Conf', 'RecentPaths.vtx')
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
		#self.ImportOptions = ImportDialogue()
		self.ImportOptions = ImportDataDialogue()
		self.ImportOptions.show()

	@pyqtSlot()
	def on_actionGenerate_from_file_triggered(self):  # how to activate menu and toolbar actions!!!
		fileNames = openFiles(self, 'seq')
		if fileNames != '' and fileNames != None:
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
		self.DialScroll(value, True)

	@pyqtSlot(int)
	def DialScroll(self, value, update):
		global from_table

		records = len(NameIndex)
		if value < records and value > -1:
			self.updateF(value)
			if update:
				if from_table:
					return
				else:
					NameIs = self.MatchingValue(value)
					self.select_tree_by_name(NameIs)
					self.tree_to_table_selection()

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
			pass
			#FontIs = self.ui.tableView.font()
			#font = QFont(FontIs)

			#FontSize = int(font.pointSize())
			#if FontSize > 7:
			#	FontSize += 1
			#font.setPointSize(FontSize)
			#font.setFamily('Lucida Grande')

			#self.ui.tableView.setFont(font)
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
			pass
			#FontIs = self.ui.tableView.font()
			#font = QFont(FontIs)

			#FontSize = int(font.pointSize())
			#if FontSize > 7:
			#	FontSize -= 1
			#font.setPointSize(FontSize)
			#font.setFamily('Lucida Grande')

			#self.ui.tableView.setFont(font)

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

		Selected = self.ui.treeWidget.selectedItems()
		Selected = Selected[-1]
		name = Selected.text(0)
		currentRow = NameIndex[name]

		records = len(NameIndex)

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

		if currentRow == -1:
			currentRow = 0

		name = list(NameIndex.keys())[list(NameIndex.values()).index(int(currentRow))]
		found = self.ui.treeWidget.findItems(name, Qt.MatchRecursive, 0)
		if len(found) > 0:
			found = found[0]
			self.ui.treeWidget.setCurrentItem(found)
		self.tree_to_table_selection()
		self.ui.radioButtonSeqView.setChecked(True)

		try:
			self.SeqButton(LastPushed)
		except:
			self.SeqButton('v')
		JustMoved = False

	def createView(self, model):
		# def createView(title, model):
		global views
		views = []
		view = self.ui.tableView
		# self.ui.tableView.keyPressEvent()
		views.append(view)
		view.setModel(model)

	def LoadDB(self, DBFilename):
		global FieldList
		global FieldCommentList
		global FieldTypeList
		global RealNameList

		if not self.createConnection(DBFilename):
			sys.exit(1)

		VGenesSQL.checkFieldTable(DBFilename)
		SQLSATEMENT = 'SELECT * FROM fieldsname ORDER BY ID'
		Records = VGenesSQL.RunSQL(DBFilename, SQLSATEMENT)

		FieldList = [i[1] for i in Records]
		RealNameList = [i[2] for i in Records]
		FieldTypeList = [i[3] for i in Records]
		FieldCommentList = [i[4] for i in Records]

		a = FieldList
		b = RealNameList
		c = FieldTypeList
		d = FieldCommentList

		self.OnOpen()
		titletext = 'VGenes - ' + DBFilename
		self.setWindowTitle(titletext)

		self.load_table()
		self.initial_Clone()

	def GenerateNameIndex(self):
		global NameIndex
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
		#self.ui.tableView.resizeColumnsToContents()

		self.ui.lcdNumber_current.display(1)
		# self.ui.txtGotoRecord.setPlainText('1')
		#model = self.ui.tableView.model()

		#index = model.index(0, 0)
		#self.ui.tableView.setCurrentIndex(index)
		#self.GenerateNameIndex()

		self.ui.cboFindField.clear()
		self.ui.cboFindField1.clear()
		self.ui.cboTreeOp1.clear()
		self.ui.cboTreeOp2.clear()
		self.ui.cboTreeOp3.clear()

		'''
		self.ui.cboFindField.addItem("Project")
		self.ui.cboFindField.addItem("Grouping")
		self.ui.cboFindField.addItem("Subgroup")

		self.ui.cboFindField1.addItem("Project")
		self.ui.cboFindField1.addItem("Grouping")
		self.ui.cboFindField1.addItem("Subgroup")

		for item in RealNameList:
			self.ui.cboFindField.addItem(item)
			self.ui.cboFindField1.addItem(item)
			self.ui.cboTreeOp1.addItem(item)
			self.ui.cboTreeOp2.addItem(item)
			self.ui.cboTreeOp3.addItem(item)
		'''

		field_list = [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		self.ui.cboFindField.addItems(field_list)
		self.ui.cboFindField1.addItems(field_list)
		self.ui.cboTreeOp1.addItems(field_list)
		self.ui.cboTreeOp2.addItems(field_list)
		self.ui.cboTreeOp3.addItems(field_list)

		self.ui.cboTreeOp1.addItem('None')
		self.ui.cboTreeOp2.addItem('None')
		self.ui.cboTreeOp3.addItem('None')
		self.ui.cboTreeOp1.setCurrentText(FieldList[75] + '(' + RealNameList[75] + ')')
		self.ui.cboTreeOp2.setCurrentText(FieldList[76] + '(' + RealNameList[76] + ')')
		self.ui.cboTreeOp3.setCurrentText(FieldList[77] + '(' + RealNameList[77] + ')')

		try:
			SQLFields = ('Project', 'Grouping', 'SubGroup')
			self.initializeTreeView(SQLFields)
			self.updateF(-2)
			self.findTreeItem(data[0])
		except:
			return

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

	def makeFieldName(self, field):
		index = FieldList.index(field)
		cur_field = FieldList[index] + '(' + RealNameList[index] + ')'
		return cur_field

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
		#Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('V gene 1')
		self.ui.valueLine.setText(valueis)

	def on_txtDgene_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDgene.toPlainText()
		LastSelected = ('D1', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('D gene 1')
		self.ui.valueLine.setText(valueis)

	def on_txtName_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtName.toPlainText()
		LastSelected = ('SeqName', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('D gene 1')
		self.ui.valueLine.setText(valueis)

	def on_txtJgene_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJgene.toPlainText()
		LastSelected = ('J1', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('J gene 1')
		self.ui.valueLine.setText(valueis)

	# todo need to put LastSelected for rest of fields

	def on_txtVLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVLocus.toPlainText()
		LastSelected = ('VLocus', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('VLocus')
		self.ui.valueLine.setText(valueis)

	def on_txtDLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDLocus.toPlainText()
		LastSelected = ('DLocus', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('DLocus')
		self.ui.valueLine.setText(valueis)

	def on_txtJLocus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJLocus.toPlainText()
		LastSelected = ('JLocus', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('JLocus')
		self.ui.valueLine.setText(valueis)

	def on_textBarcode_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textBarcode.toPlainText()
		LastSelected = ('Blank10', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Barcode')
		self.ui.valueLine.setText(valueis)


	# if data[106] != 'Blank8':
	# 	self.ui.textCluster.setText(data[106])
	# if data[107] != 'Blank9':
	# 	Population

	def on_textCluster_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textCluster.toPlainText()
		LastSelected = ('Blank8', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('10X cluster')
		self.ui.valueLine.setText(valueis)

	def on_textEdit_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textEdit.toPlainText()
		LastSelected = ('Blank9', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Seurat cluster')
		self.ui.valueLine.setText(valueis)

	def on_txtPopulation_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtPopulation.toPlainText()
		LastSelected = ('Blank11', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Population')
		self.ui.valueLine.setText(valueis)

	def on_textMutations_selectionChanged(self):
		global LastSelected
		valueis = self.ui.textMutations.toPlainText()
		LastSelected = ('TotMut', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('TotMut')
		self.ui.valueLine.setText(valueis)

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
		self.on_btnSaveChange_clicked()

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

		#self.model.refresh()



	def on_txtProject_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtProject.toPlainText()
		LastSelected = ('Project', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Project')
		self.ui.valueLine.setText(valueis)

	def on_txtGroup_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtGroup.toPlainText()
		LastSelected = ('Grouping', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Grouping')
		self.ui.valueLine.setText(valueis)

	def on_txtGeneType_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtGeneType.toPlainText()
		LastSelected = ('GeneType', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('GeneType')
		self.ui.valueLine.setText(valueis)

	def on_txtSubGroup_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtSubGroup.toPlainText()
		LastSelected = ('SubGroup', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('SubGroup')
		self.ui.valueLine.setText(valueis)

	def on_txtVend_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVend.toPlainText()
		LastSelected = ('VSeqend', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('VSeqend')
		self.ui.valueLine.setText(valueis)

	def on_txtD_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtD.toPlainText()
		LastSelected = ('Dregion', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Dregion')
		self.ui.valueLine.setText(valueis)

	def on_txtVD_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtVD.toPlainText()
		LastSelected = ('VDJunction', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('VDJunction')
		self.ui.valueLine.setText(valueis)

	def on_txtDJ_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtDJ.toPlainText()
		LastSelected = ('DJJunction', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('DJJunction')
		self.ui.valueLine.setText(valueis)

	def on_txtJend_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtJend.toPlainText()
		LastSelected = ('begJ', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('begJ')
		self.ui.valueLine.setText(valueis)

	def on_txtCDR3DNA_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3DNA.toPlainText()
		LastSelected = ('CDR3DNA', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('CDR3DNA')
		self.ui.valueLine.setText(valueis)

	def on_txtIsotype_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtIsotype.toPlainText()
		LastSelected = ('Isotype', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Isotype')
		self.ui.valueLine.setText(valueis)


	def on_txtCDR3AA_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3AA.toPlainText()
		LastSelected = ('CDR3AA', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('CDR3AA')
		self.ui.valueLine.setText(valueis)

	def on_txtCDR3Length_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3Length.toPlainText()
		LastSelected = ('CDR3Length', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('CDR3Length')
		self.ui.valueLine.setText(valueis)

	def on_txtCDR3pI_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3pI.toPlainText()
		LastSelected = ('CDR3pI', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('CDR3pI')
		self.ui.valueLine.setText(valueis)

	def on_txtCDR3MW_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtCDR3MW.toPlainText()
		LastSelected = ('CDR3MW', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('CDR3MW')
		self.ui.valueLine.setText(valueis)

	def on_txtProductive_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtProductive.toPlainText()
		LastSelected = ('productive', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('productive')
		self.ui.valueLine.setText(valueis)

	def on_txtReadingFrame_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtReadingFrame.toPlainText()
		LastSelected = ('ReadingFrame', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('ReadingFrame')
		self.ui.valueLine.setText(valueis)

	def on_txtStop_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtStop.toPlainText()
		LastSelected = ('StopCodon', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('StopCodon')
		self.ui.valueLine.setText(valueis)

	def on_txtQuality_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtQuality.toPlainText()
		LastSelected = ('Quality', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Quality')
		self.ui.valueLine.setText(valueis)

	def on_txtQuality_2_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtQuality_2.toPlainText()
		LastSelected = ('Quality', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Quality')
		self.ui.valueLine.setText(valueis)

	def on_txtStatus_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtStatus.toPlainText()
		LastSelected = ('Blank13', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Status')
		self.ui.valueLine.setText(valueis)

	def on_txtLabel_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtLabel.toPlainText()
		LastSelected = ('Blank12', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('Label')
		self.ui.valueLine.setText(valueis)

	def on_txtID_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtID.toPlainText()
		LastSelected = ('IDEvent', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('IDEvent')
		self.ui.valueLine.setText(valueis)

	def on_txtClonalPool_selectionChanged(self):
		global LastSelected
		valueis = self.ui.txtClonalPool.toPlainText()
		LastSelected = ('ClonalPool', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

		self.ui.fieldLine.setText('ClonalPool')
		self.ui.valueLine.setText(valueis)

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

		self.ui.cboFindField1.setCurrentText(self.ui.cboFindField.currentText())

	def SpecSet(self):
		global LastSelected

		# if JustMoved == False:
		valueis = self.ui.listViewSpecificity.currentText()
		LastSelected = ('Specificity', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def SubspecSet(self):
		global LastSelected


		# if JustMoved == False:
		valueis = self.ui.listViewSpecificity_2.currentText()
		LastSelected = ('Subspecificity', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
		self.ui.cboFindField.setCurrentText(Fiedlvalue)

	def AutoRXSet(self):
		global LastSelected


		# if JustMoved == False:
		valueis = self.ui.Autoreactivity.currentText()
		LastSelected = ('Blank6', valueis)
		field = LastSelected[0]
		# Fiedlvalue = self.TransLateFieldtoReal(field, False)
		Fiedlvalue = self.makeFieldName(field)
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
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		self.clearTreeChecks()
		self.match_tree_to_table()
		self.ui.checkBoxAll.setChecked(False)

	@pyqtSlot()
	def on_actionSuggestCanonical_triggered(self):
		self.ui.treeWidget.expandAll()

		# fields = self.ui.cboTreeOp1.currentText()
		# field1 = self.TransLateFieldtoReal(fields, True)
		field1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
		i = 0
		for item in FieldList:
			if field1 == item:
				field1Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp2.currentText()
		# field2 = self.TransLateFieldtoReal(fields, True)
		field2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
		i = 0
		for item in FieldList:
			if field2 == item:
				field2Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp3.currentText()
		# field3 = self.TransLateFieldtoReal(fields, True)
		field3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
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

			#model = self.ui.tableView.model()
			#model.refresh()

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

		# fields = self.ui.cboTreeOp1.currentText()
		# field1 = self.TransLateFieldtoReal(fields, True)
		field1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
		i = 0
		for item in FieldList:
			if field1 == item:
				field1Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp2.currentText()
		# field2 = self.TransLateFieldtoReal(fields, True)
		field2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
		i = 0
		for item in FieldList:
			if field2 == item:
				field2Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp3.currentText()
		# field3 = self.TransLateFieldtoReal(fields, True)
		field3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
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
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		global wasClicked
		search = self.ui.txtFieldSearch.text()
		field = self.ui.cboFindField.currentText()
		fieldsearch = re.sub(r'\(.+', '', field)

		#     if search != '':
		#         SQLstatement = 'select * from vgenesdb WHERE ' + fieldsearch + ' LIKE ' + search + ' ORDER BY ' + fieldsearch + ', Project, Grouping, SubGroup, SeqName'
		# #         todo select one if exact or if no then all similar items in tree and table...best if tree checkable

		if self.ui.rdoLocal.isChecked():
			# fields = self.ui.cboTreeOp1.currentText()
			# field1 = self.TransLateFieldtoReal(fields, True)
			field1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
			i = 0
			for item in FieldList:
				if field1 == item:
					field1Value = data[i]
				i += 1

			# fields = self.ui.cboTreeOp2.currentText()
			# field2 = self.TransLateFieldtoReal(fields, True)
			field2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
			i = 0
			for item in FieldList:
				if field2 == item:
					field2Value = data[i]
				i += 1

			# fields = self.ui.cboTreeOp3.currentText()
			# field3 = self.TransLateFieldtoReal(fields, True)
			field3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
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
	def on_pushButtonBatch_clicked(self):
		global DontFindTwice
		cur_field = re.sub(r'\(.+', '', self.ui.cboFindField1.currentText())
		'''
		if cur_field in RealNameList:
			index = RealNameList.index(cur_field)
			cur_field = FieldList[index]
		else:
			answer = informationMessage(self, 'The field name you determined is not exist!', 'OK')
			return
		'''
		SQLStatement = 'SELECT DISTINCT(' + cur_field + ') FROM vgenesdb'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		value_list = [row[0] for row in DataIn]

		if len(value_list) > 30:
			question = 'Distinct values of this field seems too many (number =  ' + str(len(value_list)) + ')\nAre you sure?'
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)
			if answer == 'No':
				return

		self.myBatchDialog = BatchDialog()
		self.myBatchDialog.initial = 0
		#self.myBatchDialog.load_data(value_list)
		field_list = [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		self.myBatchDialog.ui.comboBox.addItems(field_list)
		self.myBatchDialog.initial = 1
		self.myBatchDialog.StatFig()
		self.myBatchDialog.ui.comboBox.setCurrentText(self.ui.cboFindField1.currentText())
		self.myBatchDialog.BatchSignal.connect(self.updateFieldBatch)
		self.myBatchDialog.show()
		self.myBatchDialog.initial = 2
		self.myBatchDialog.resize(1200, 700)
	
	@pyqtSlot()
	def on_StatUpdate_clicked(self):
		self.myBatchDialog = BatchDialog()
		self.myBatchDialog.initial = 0
		# self.myBatchDialog.load_data(value_list)
		field_list = [FieldList[i] + '(' + RealNameList[i] + ')' for i in range(len(FieldList))]
		self.myBatchDialog.ui.comboBox.addItems(field_list)
		self.myBatchDialog.initial = 1
		self.myBatchDialog.StatFig()
		self.myBatchDialog.ui.comboBox.setCurrentText(self.ui.cboFindField1.currentText())
		self.myBatchDialog.BatchSignal.connect(self.updateFieldBatch)
		self.myBatchDialog.show()
		self.myBatchDialog.initial = 2
		self.myBatchDialog.resize(1200, 700)

	def updateFieldBatch(self, indicator, field, dict):
		if indicator == 1:
			Msg = 'You did not specific any value! Do nothing!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
		else:
			for key in dict:
				SQLStatement = 'UPDATE vgenesdb SET ' + field + ' = "' + dict[key] + '" WHERE ' + field + ' = "' + key + '"'
				print(SQLStatement)
				VGenesSQL.RunUpdateSQL(DBFilename, SQLStatement)
			Msg = 'Update finished!'
			QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)

			self.refreshDB()
			value = self.ui.dial.value()
			self.updateF(value)

	@pyqtSlot()
	def on_btnFieldBulk_clicked(self):
		# find ID for all selected items
		selected_name = self.getTreeCheckedChild()
		selected_name = selected_name[3]

		if len(selected_name) == 0:
			Msg = 'Please check at least one record you want to edit!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

		SQLStatement = 'SELECT ID FROM vgenesdb WHERE SeqName IN ("' + '","'.join(selected_name) + '")'
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)

		cur_field = re.sub(r'\(.+', '', self.ui.cboFindField1.currentText())
		'''
		if cur_field in RealNameList:
			index = RealNameList.index(cur_field)
			cur_field = FieldList[index]
		else:
			answer = informationMessage(self, 'The field name you determined is not exist!', 'OK')
			return
		'''

		str_target = self.ui.lineEditTarget.text()
		if str_target == "":
			question = "Your target field is empty, will update the whole cell, OK?"
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)
			if answer == 'No':
				return

		str_replace = self.ui.lineEditReplace.text()
		if str_replace == '':
			question = "Your replace field is empty, will clear everything match target field, OK?"
			buttons = 'YN'
			answer = questionMessage(self, question, buttons)
			if answer == 'No':
				return

		# update each records
		for record in DataIn:
			id = record[0]
			SQLStatement = ""
			if str_target == "":
				SQLStatement = "UPDATE vgenesdb SET " + cur_field + " = '" + str_replace + "' WHERE ID = " + str(id)
			else:
				SQLStatement = 'SELECT ' + cur_field + ' FROM vgenesdb WHERE ID = ' + str(id)
				thisData = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				thisData = thisData[0][0]

				if self.ui.radioButtonReg.isChecked():
					regx = re.compile(str_target)
					new_str = re.sub(regx, str_replace, thisData)
				else:
					new_str = re.sub(str_target, str_replace, thisData)

				if new_str != thisData:
					SQLStatement = "UPDATE vgenesdb SET " + cur_field + " = '" + new_str + "' WHERE ID = " + str(id)

			try:
				VGenesSQL.RunUpdateSQL(DBFilename, SQLStatement)
			except:
				Msg = 'SQL ERROR! Current CMD:\n' + SQLStatement
				QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
				return

		#SQLFields = (
		#	self.ui.cboTreeOp1.currentText(), self.ui.cboTreeOp2.currentText(), self.ui.cboTreeOp3.currentText())
		'''
		index1 = RealNameList.index(self.ui.cboTreeOp1.currentText())
		index2 = RealNameList.index(self.ui.cboTreeOp2.currentText())
		index3 = RealNameList.index(self.ui.cboTreeOp3.currentText())

		SQLFields = (FieldList[index1], FieldList[index2], FieldList[index3])
		'''
		SQLFields = (
			re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText()),
			re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText()),
			re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
		)

		self.initializeTreeView(SQLFields)
		self.ui.treeWidget.expandAll()
		self.refreshDB()

	@pyqtSlot()
	def on_btnEditLock_clicked(self):
		if self.enableEdit == True:
			lock_icon = QIcon()
			lock_icon.addPixmap(QPixmap(":/PNG-Icons/locked.png"), QIcon.Normal, QIcon.Off)
			self.ui.btnEditLock.setIcon(lock_icon)

			self.ui.btnSaveChange.setEnabled(False)

			self.ui.txtDateTime.setReadOnly(True)
			self.ui.txtGroup.setReadOnly(True)
			self.ui.txtLabel.setReadOnly(True)
			self.ui.txtName.setReadOnly(True)
			self.ui.txtProject.setReadOnly(True)
			self.ui.txtQuality.setReadOnly(True)
			self.ui.txtStatus.setReadOnly(True)
			self.ui.txtSubGroup.setReadOnly(True)
			self.ui.comboBoxSpecies.setEnabled(False)
			
			self.ui.txtDLocus.setReadOnly(True)
			self.ui.txtDgene.setReadOnly(True)
			self.ui.txtVLocus.setReadOnly(True)
			self.ui.txtVgene.setReadOnly(True)
			self.ui.txtJLocus.setReadOnly(True)
			self.ui.txtJgene.setReadOnly(True)
			self.ui.txtIsotype.setReadOnly(True)
			self.ui.txtID.setReadOnly(True)
			self.ui.txtStop.setReadOnly(True)

			self.ui.textBarcode.setReadOnly(True)
			self.ui.textCluster.setReadOnly(True)
			self.ui.textMutations.setReadOnly(True)
			self.ui.textEdit.setReadOnly(True)
			self.ui.txtPopulation.setReadOnly(True)

			self.ui.txtD.setReadOnly(True)
			self.ui.txtDJ.setReadOnly(True)
			self.ui.txtJend.setReadOnly(True)
			self.ui.txtVD.setReadOnly(True)
			self.ui.txtVend.setReadOnly(True)
			self.ui.txtClonalPool.setReadOnly(True)
			self.ui.txtClonalRank.setReadOnly(True)
			self.ui.txtProductive.setReadOnly(True)
			self.ui.txtReadingFrame.setReadOnly(True)
			self.ui.txtCDR3AA.setReadOnly(True)
			self.ui.txtCDR3DNA.setReadOnly(True)
			self.ui.txtCDR3Length.setReadOnly(True)
			self.ui.txtCDR3MW.setReadOnly(True)
			self.ui.txtCDR3pI.setReadOnly(True)

			self.ui.radioButton_21.setEnabled(False)
			self.ui.radioButton_22.setEnabled(False)
			self.ui.radioButton_23.setEnabled(False)
			self.ui.Autoreactivity.setEnabled(False)
			self.ui.listViewSpecificity.setEnabled(False)
			self.ui.listViewSpecificity_2.setEnabled(False)

			self.ui.txtComments.setReadOnly(True)

			self.enableEdit = False

		else:
			unlock_icon = QIcon()
			unlock_icon.addPixmap(QPixmap(":/PNG-Icons/unlocked.png"), QIcon.Normal, QIcon.Off)
			self.ui.btnEditLock.setIcon(unlock_icon)

			self.ui.btnSaveChange.setEnabled(True)

			self.ui.txtDateTime.setReadOnly(False)
			self.ui.txtGroup.setReadOnly(False)
			self.ui.txtLabel.setReadOnly(False)
			self.ui.txtName.setReadOnly(False)
			self.ui.txtProject.setReadOnly(False)
			self.ui.txtQuality.setReadOnly(False)
			self.ui.txtStatus.setReadOnly(False)
			self.ui.txtSubGroup.setReadOnly(False)
			self.ui.comboBoxSpecies.setEnabled(True)

			self.ui.txtDLocus.setReadOnly(False)
			self.ui.txtDgene.setReadOnly(False)
			self.ui.txtVLocus.setReadOnly(False)
			self.ui.txtVgene.setReadOnly(False)
			self.ui.txtJLocus.setReadOnly(False)
			self.ui.txtJgene.setReadOnly(False)
			self.ui.txtIsotype.setReadOnly(False)
			self.ui.txtID.setReadOnly(False)
			self.ui.txtStop.setReadOnly(False)

			self.ui.textBarcode.setReadOnly(False)
			self.ui.textCluster.setReadOnly(False)
			self.ui.textMutations.setReadOnly(False)
			self.ui.textEdit.setReadOnly(False)
			self.ui.txtPopulation.setReadOnly(False)

			self.ui.txtD.setReadOnly(False)
			self.ui.txtDJ.setReadOnly(False)
			self.ui.txtJend.setReadOnly(False)
			self.ui.txtVD.setReadOnly(False)
			self.ui.txtVend.setReadOnly(False)
			self.ui.txtClonalPool.setReadOnly(False)
			self.ui.txtClonalRank.setReadOnly(False)
			self.ui.txtProductive.setReadOnly(False)
			self.ui.txtReadingFrame.setReadOnly(False)
			self.ui.txtCDR3AA.setReadOnly(False)
			self.ui.txtCDR3DNA.setReadOnly(False)
			self.ui.txtCDR3Length.setReadOnly(False)
			self.ui.txtCDR3MW.setReadOnly(False)
			self.ui.txtCDR3pI.setReadOnly(False)

			self.ui.radioButton_21.setEnabled(True)
			self.ui.radioButton_22.setEnabled(True)
			self.ui.radioButton_23.setEnabled(True)
			self.ui.Autoreactivity.setEnabled(True)
			self.ui.listViewSpecificity.setEnabled(True)
			self.ui.listViewSpecificity_2.setEnabled(True)

			self.ui.txtComments.setReadOnly(False)

			self.enableEdit = True

	@pyqtSlot()
	def on_pushButtonLock_clicked(self):
		if self.ui.btnFieldBulk.isEnabled() == True:
			lock_icon = QIcon()
			lock_icon.addPixmap(QPixmap(":/PNG-Icons/locked.png"), QIcon.Normal, QIcon.Off)
			self.ui.pushButtonLock.setIcon(lock_icon)
			self.ui.btnFieldBulk.setEnabled(False)
			self.ui.pushButtonBatch.setEnabled(False)
		else:
			unlock_icon = QIcon()
			unlock_icon.addPixmap(QPixmap(":/PNG-Icons/unlocked.png"), QIcon.Normal, QIcon.Off)
			self.ui.pushButtonLock.setIcon(unlock_icon)
			self.ui.btnFieldBulk.setEnabled(True)
			self.ui.pushButtonBatch.setEnabled(True)

	@pyqtSlot()
	def on_btnSaveChange_clicked(self):
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return

		global updateMarker
		if self.enableEdit == True:
			SETStatement = 'SET '
			SETStatement += 'Project = "' + self.ui.txtProject.toPlainText() + '",'
			SETStatement += 'Grouping = "' + self.ui.txtGroup.toPlainText() + '",'
			SETStatement += 'SubGroup = "' + self.ui.txtSubGroup.toPlainText() + '",'
			SETStatement += 'SeqName = "' + self.ui.txtName.toPlainText() + '",'
			SETStatement += 'Species = "' + self.ui.comboBoxSpecies.currentText() + '",'
			SETStatement += 'DateEntered = "' + self.ui.txtDateTime.toPlainText() + '",'
			SETStatement += 'Quality = "' + self.ui.txtQuality.toPlainText() + '",'
			SETStatement += 'Blank12 = "' + self.ui.txtLabel.toPlainText() + '",'
			SETStatement += 'Blank12 = "' + self.ui.txtLabel.toPlainText() + '",'
			SETStatement += 'Blank13 = "' + self.ui.txtStatus.toPlainText() + '",'
			SETStatement += 'V1 = "' + self.ui.txtVgene.toPlainText() + '",'
			SETStatement += 'VLocus = "' + self.ui.txtVLocus.toPlainText() + '",'
			SETStatement += 'D1 = "' + self.ui.txtDgene.toPlainText() + '",'
			SETStatement += 'DLocus = "' + self.ui.txtDLocus.toPlainText() + '",'
			SETStatement += 'J1 = "' + self.ui.txtJgene.toPlainText() + '",'
			SETStatement += 'JLocus = "' + self.ui.txtJLocus.toPlainText() + '",'
			SETStatement += 'Isotype = "' + self.ui.txtIsotype.toPlainText() + '",'
			SETStatement += 'IDEvent = "' + self.ui.txtID.toPlainText() + '",'
			SETStatement += 'StopCodon = "' + self.ui.txtStop.toPlainText() + '",'
			SETStatement += 'Blank8 = "' + self.ui.textCluster.toPlainText() + '",'
			SETStatement += 'Blank9 = "' + self.ui.textEdit.toPlainText() + '",'
			SETStatement += 'TotMut = "' + self.ui.textMutations.toPlainText() + '",'
			SETStatement += 'Blank10 = "' + self.ui.textBarcode.toPlainText() + '",'
			SETStatement += 'Blank11 = "' + self.ui.txtPopulation.toPlainText() + '",'
			SETStatement += 'VSeqend = "' + self.ui.txtVend.toPlainText() + '",'
			SETStatement += 'VDJunction = "' + self.ui.txtVD.toPlainText() + '",'
			SETStatement += 'Dregion = "' + self.ui.txtD.toPlainText() + '",'
			SETStatement += 'DJJunction = "' + self.ui.txtDJ.toPlainText() + '",'
			SETStatement += 'begJ = "' + self.ui.txtJend.toPlainText() + '",'
			SETStatement += 'CDR3DNA = "' + self.ui.txtCDR3DNA.toPlainText() + '",'
			SETStatement += 'CDR3AA = "' + self.ui.txtCDR3AA.toPlainText() + '",'
			SETStatement += 'CDR3Length = "' + self.ui.txtCDR3Length.toPlainText() + '",'
			SETStatement += 'CDR3pI = "' + self.ui.txtCDR3pI.toPlainText() + '",'
			SETStatement += 'CDR3MW = "' + self.ui.txtCDR3MW.toPlainText() + '",'
			SETStatement += 'productive = "' + self.ui.txtProductive.toPlainText() + '",'
			SETStatement += 'ReadingFrame = "' + self.ui.txtReadingFrame.toPlainText() + '",'
			SETStatement += 'ClonalPool = "' + self.ui.txtClonalPool.toPlainText() + '",'
			SETStatement += 'ClonalRank = "' + self.ui.txtClonalRank.toPlainText() + '",'
			SETStatement += 'Specificity = "' + self.ui.listViewSpecificity.currentText() + '",'
			SETStatement += 'Subspecificity = "' + self.ui.listViewSpecificity_2.currentText() + '",'
			SETStatement += 'Blank6 = "' + self.ui.Autoreactivity.currentText() + '",'
			SETStatement += 'Comments = "' + self.ui.txtComments.toPlainText() + '" '

			value = self.ui.treeWidget.selectedItems()
			name = value[-1].text(0)
			WHEREStatement = 'WHERE SeqName = "' + name + '"'

			SQLStatement = 'UPDATE vgenesDB ' + SETStatement + WHEREStatement
			foundRecs = VGenesSQL.UpdateMulti(SQLStatement, DBFilename)

			if name == self.ui.txtName.toPlainText():
				pass
			else:
				# update name index
				NameIndex[self.ui.txtName.toPlainText()] = NameIndex[name]
				del NameIndex[name]
				# refresh model
				#model = self.ui.tableView.model()
				#model.refresh()
				# update database marker
				updateMarker = True
				# update tree
				#SQLFields = (
				#	self.ui.cboTreeOp1.currentText(), self.ui.cboTreeOp2.currentText(), self.ui.cboTreeOp3.currentText())
				'''
				index1 = RealNameList.index(self.ui.cboTreeOp1.currentText())
				index2 = RealNameList.index(self.ui.cboTreeOp2.currentText())
				index3 = RealNameList.index(self.ui.cboTreeOp3.currentText())

				SQLFields = (FieldList[index1], FieldList[index2], FieldList[index3])
				'''
				SQLFields = (
					re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText()),
					re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText()),
					re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
				)

				self.initializeTreeView(SQLFields)
				self.ui.treeWidget.expandAll()

			Msg = 'Change saved!'
			QMessageBox.information(self, 'Information', Msg, QMessageBox.Ok, QMessageBox.Ok)
		else:
			Msg = 'Please switch to edit mode first!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
			return

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

		elif option == 'CSV format Entire VDB':
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
		Backfilename = os.path.join(working_prefix, 'BackUP.vdb')
		global DBFilename

		try:
			if DBFilename != None:
				shutil.copy(DBFilename, Backfilename)
		except:
			return

	@pyqtSlot()
	def on_actionRevert_to_previous_triggered(self):
		import shutil

		buttons = 'OKC'
		answer = informationMessage(self, 'Any work done since opening this instance of VGenes will be lost.', buttons)
		if answer == 'Cancel':
			return

		Backfilename = os.path.join(working_prefix, 'BackUP.vdb')

		self.GOOpen(False)
		shutil.move(Backfilename, DBFilename)
		self.GOOpen(False)

		self.refreshDB()
		Msg = 'VGene Database resumed!'
		QMessageBox.information(self, 'information', Msg, QMessageBox.Ok,
		                        QMessageBox.Ok)

	@pyqtSlot()
	def on_actionExport_triggered(self):
		buttons = 'OKC'
		answer = informationMessage(self, 'This function will only export FASTA format sequences\nFor more formats, please use Generate report function.\nContinue?', buttons)
		if answer == 'Cancel':
			return

		selected = self.getTreeCheckedChild()
		selected = selected[3]

		if len(selected) == 0:
			Msg = 'Please select at least one sequence!'
			QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok,
			                        QMessageBox.Ok)
			return

		SQLStatement = 'SELECT SeqName,Sequence FROM vgenesdb WHERE SeqName IN ("' + '","'.join(selected) + '")'
		
		filename = saveFile(self, 'Nucleotide')
		if filename == '' or filename == 'none':
			return

		f = open(filename, 'w')
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)
		for record in DataIn:
			f.write('>' + record[0] + '\n' + record[1] + '\n')
		f.close()

		Msg = 'Sequences exported!'
		QMessageBox.information(self, 'information', Msg, QMessageBox.Ok,
		                        QMessageBox.Ok)

	@pyqtSlot()
	def on_btnExtractRecords_clicked(self):
		New = True
		self.MoveRecords(New)

	@pyqtSlot()
	def on_action_Help_triggered(self):
		Msg = 'VGene was developed and supported by Wilson Lab. Please refer to \nhttp://Wilsonlab.uchicago.edu\n' \
		      'for more information'
		QMessageBox.information(self, 'information', Msg, QMessageBox.Ok,
		                        QMessageBox.Ok)

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

	@pyqtSlot()
	def on_btnUpdateTree_clicked(self):
		if DBFilename == '' or DBFilename == 'none' or DBFilename == None:
			return
		'''
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
		'''

		# fields = self.ui.cboTreeOp1.currentText()
		# field1 = self.TransLateFieldtoReal(fields, True)
		field1 = re.sub(r'\(.+', '', self.ui.cboTreeOp1.currentText())
		i = 0
		for item in FieldList:
			if field1 == item:
				field1Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp2.currentText()
		# field2 = self.TransLateFieldtoReal(fields, True)
		field2 = re.sub(r'\(.+', '', self.ui.cboTreeOp2.currentText())
		i = 0
		for item in FieldList:
			if field2 == item:
				field2Value = data[i]
			i += 1

		# fields = self.ui.cboTreeOp3.currentText()
		# field3 = self.TransLateFieldtoReal(fields, True)
		field3 = re.sub(r'\(.+', '', self.ui.cboTreeOp3.currentText())
		i = 0
		for item in FieldList:
			if field3 == item:
				field3Value = data[i]
			i += 1

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

		#model = self.ui.tableView.model()

		# global RefreshSQL
		if field1 == 'None' or field1 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY SeqName'
		elif field2 == 'None' or field2 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', SeqName'
		elif field3 == 'None' or field3 is None:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', SeqName'
		else:
			RefreshSQL = 'select * from vgenesdb ORDER BY ' + field1 + ', ' + field2 + ', ' + field3 + ', SeqName'

		#model.refresh()

		#self.ui.tableView.sortByColumn(field1Index, Qt.AscendingOrder)
		#self.ui.tableView.sortByColumn(0, Qt.AscendingOrder)
		#self.ui.tableView.sortByColumn(field3Index, Qt.AscendingOrder)
		#self.ui.tableView.sortByColumn(field2Index, Qt.AscendingOrder)

		self.initializeTreeView(SQLFields)
		global DontFindTwice
		DontFindTwice = True
		self.findTreeItem(currentitemIs)
		# self.ui.treeWidget.
		DontFindTwice = False

		self.ui.checkBoxAll.setChecked(False)
		self.ui.checkBoxAll1.setChecked(False)

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

		# clear fields\
		global JustMovedIt
		global FirstupdateF
		JustMovedIt = True
		# JustMoved = True
		if PreVID != ID:
			if ID != -1:
				data.clear()

				self.ui.txtVbeg.clear()
				self.ui.txtVExp.clear()
				self.ui.txtJend_2.clear()
				self.ui.txtJExp.clear()

				#if FirstupdateF == False and ID > -1:
				# MatchingIndex = NameIndex[name]
				if ID == -2: ID = 0
				try:
					newID = list(NameIndex.keys())[list(NameIndex.values()).index(ID)]
				except:
					return
				SQLStatement = 'SELECT * FROM vgenesDB WHERE SeqName = "' + str(newID) + '"'
				DataIs = VGenesSQL.RunSQL(DBFilename, SQLStatement)
				for record in DataIs:
					for item in record:
						data.append(str(item))
				#else:
				#	model = self.ui.tableView.model()
				#	if ID == -2: ID = 0
				#	for i in range(0, 120):
				#		index = model.index(ID, i)
				#		data.append(str(model.data(index)))

				#	FirstupdateF = False

				PreVID = ID

				self.ui.txtName.setText(data[0])
				self.ui.txtName_2.setText(data[0])
				self.ui.label_Name.setText(data[0])

				self.ui.txtGeneType.setText(data[2])
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
					Dbeg = int(data[69])
					self.ui.sbDbeg.setValue(Dbeg)
					Dend = int(data[70])
					self.ui.sbDend.setValue(Dend)
				except:
					print('D error')

				try:
					self.ui.txtVbeg.setText(VbegDisplay)
					self.ui.txtVExp.setText(GVbegDisplay)

					Jbeg = int(data[73])
					self.ui.sbJbeg.setValue(Jbeg)
					Jend = int(data[74])
					GJend = int(data[66])
					self.ui.sbJend.setValue(Jend)

					JendSeq = VSeq[Jend-11:Jend]    # why last 11bp?
					JendAASeq, ErMessage = VGenesSeq.Translator(JendSeq, 0)
					JendDisplay = ' ' + JendAASeq[0] + '   ' + JendAASeq[1] + '   ' + JendAASeq[2] + ' \n' + JendSeq[0:3] + ' ' + JendSeq[3:6] + ' ' + JendSeq[6:9]

					self.ui.txtJend_2.setText(JendDisplay)
				except:
					print('J error')
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

				#self.ui.cboTreeOp1.setCurrentText(data[75])
				#self.ui.cboTreeOp2.setCurrentText(data[76])
				#self.ui.cboTreeOp3.setCurrentText(data[77])
				self.ui.cboTreeOp1.setCurrentText(FieldList[75] + '(' + RealNameList[75] + ')')
				self.ui.cboTreeOp2.setCurrentText(FieldList[76] + '(' + RealNameList[76] + ')')
				self.ui.cboTreeOp3.setCurrentText(FieldList[77] + '(' + RealNameList[77] + ')')

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
				#currentRecord = self.ui.tableView.currentIndex().row()
				#maxRecords = self.ui.tableView.model().rowCount()
				#self.ui.horizontalScrollBar.setMaximum(maxRecords)
				#self.ui.dial.setMaximum(maxRecords)

				currentRecord = ID
				maxRecords = len(NameIndex)
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

	@pyqtSlot()
	def on_actionCreateVDJdb_triggered(self):

		Pathname = openFile(self, 'Nucleotide')
		if Pathname == None:
			return
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
			# clear check first
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

		AAStartSel = math.floor(StartSel / 3)
		AAEndSel = math.floor(EndSel / 3)

		print(str(StartSel) + ',' + str(EndSel) + ',' + str(AAStartSel) + ',' + str(AAEndSel))

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
	def SeqButtonold(self, button):
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

					#model = self.ui.tableView.model()
					#model.refresh()

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
		ErlogFile = os.path.join(temp_folder,'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
		ErLog = 'VGenes input beginning at: ' + time.strftime('%c') + '\n'
		with open(ErlogFile, 'w') as currentFile:  # using with for this automatically closes the file even if you crash
			currentFile.write(ErLog)

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

			self.ui.txtFieldSearch.setText(HSeqName)
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

			self.ui.txtFieldSearch.setText(LSeqName)
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
		#currentRow = self.ui.tableView.currentIndex().row()
		# todo change to app folder
		try:
			filename = os.path.join(working_prefix, 'UpdateRecord.nt')
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
		#model = self.ui.tableView.model()
		DataRow = int(data[119])

		for record in IgBLASTAnalysis:
			for item in record:
				FieldName = FieldList[i]
				ItemValue = str(item) #+ ' '
				if i != 119:
					VGenesSQL.UpdateField(DataRow, ItemValue, FieldName, DBFilename)
				i += 1


		#model.refresh()

		global PreVID
		# NewID = int(data[119])+1
		# PreVID = NewID

		self.updateF(data[119])

	def UpdateMutationAnalysis(self):
		datalist = []
		#currentRow = self.ui.tableView.currentIndex().row()
		# todo change to app folder
		try:
			filename = os.path.join(working_prefix, 'UpdateRecord.nt')
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
			#model = self.ui.tableView.model()
			# DataRow = int(data[119])

			for record in IgBLASTAnalysis:
				for item in record:
					FieldName = FieldList[i]
					ItemValue = str(item)  # + ' '
					if i == 57 or i == 96 or i == 97 or i == 98:  #only SHM fields
						VGenesSQL.UpdateField(DataRow, ItemValue, FieldName, DBFilename)

					i += 1

			#model.refresh()

		self.updateF(data[119])

def IgBlastParserFast(FASTAFile, datalist):
	import os
	# todo change to app folder
	progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
	ErlogFile = os.path.join(temp_folder,  'ErLog.txt')
	ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
	ErLog = 'VGenes input beginning at: ' + time.strftime('%c') + '\n'
	with open(ErlogFile, 'w') as currentFile:  # using with for this automatically closes the file even if you crash
		currentFile.write(ErLog)

	try:
		DBpathname = os.path.join(working_prefix, 'Data', 'VDJGenes.db')
		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)
		GetProductive = False
		conn = db.connect(DBpathname)
	except:
		Msg = 'VDJ DB connect Error!'
		return Msg

	#  then need to create a cursor that lets you traverse the database
	cursor = conn.cursor()
	DATA = []
	Sequences = {}

	project = datalist[0]
	grouping = datalist[1]
	subgroup = datalist[2]
	species = datalist[3]
	GetProductive = datalist[4]
	MaxNum = int(datalist[5])

	ErLog, TotSeqs = IgBLASTer.ProcessFASTA(FASTAFile, MaxNum)
	workingdir = os.path.join(working_prefix, 'IgBlast')
	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt')

	os.chdir(workingdir)
	Sequences.clear()
	now = time.strftime('%c')
	# read input sequence file
	with open(workingfilename, 'r') as currentFile:  #make dictionary of seqs keyed by name
		for FASTAline in currentFile:
			FASTAline = FASTAline.replace('\n', '').replace('\r', '')
			if FASTAline[0] == '>':
				#print(FASTAline)
				SeqNamed  = FASTAline[1:]
				SeqNamed = SeqNamed.strip()

			else:
				Sequence  = FASTAline
				if Sequence != '':
					Sequences[SeqNamed] = Sequence
				SeqNamed = ''
				Sequence = ''

	totel_seq = len(Sequences)

	# run IgBlast
	time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())
	igblast_out_fmt19 = os.path.join(temp_folder, time_stamp + '_igblast_out_fmt19.csv')
	try:
		start = time.time()
		if species == 'Human':
			BLASTCommandLine = workingdir + "/igblastn -germline_db_V HumanVGenes.nt -germline_db_J HumanJGenes.nt -germline_db_D HumanDGenes.nt -organism human -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt 3"
			IgBlastOut_fmt3 = os.popen(BLASTCommandLine)
			print(BLASTCommandLine)
			BLASTCommandLine = workingdir + "/igblastn -germline_db_V HumanVGenes.nt -germline_db_J HumanJGenes.nt -germline_db_D HumanDGenes.nt -organism human -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt 19 > " + igblast_out_fmt19
			IgBlastOut_fmt19 = os.system(BLASTCommandLine)
			print(BLASTCommandLine)
		elif species == 'Mouse':
			BLASTCommandLine = workingdir + "/igblastn -germline_db_V MouseVGenes.nt -germline_db_J MouseJGenes.nt -germline_db_D MouseDGenes.nt -organism mouse -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3"
			IgBlastOut_fmt3 = os.popen(BLASTCommandLine)
			print(BLASTCommandLine)
			BLASTCommandLine = workingdir + "/igblastn -germline_db_V MouseVGenes.nt -germline_db_J MouseJGenes.nt -germline_db_D MouseDGenes.nt -organism mouse -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 19 > " + igblast_out_fmt19
			IgBlastOut_fmt19 = os.system(BLASTCommandLine)
			print(BLASTCommandLine)
		end = time.time()
		print('Run time for IgBlast: ' + str(end - start))
	except:
		ErLog = 'VGenes running Error!\nCurrent CMD: ' + BLASTCommandLine + '\n'
		with open(ErlogFile, 'a') as currentFile:  # using with for this automatically closes the file even if you crash
			currentFile.write(ErLog)

	# parse IgBlast out fmt 19
	with open(igblast_out_fmt19, 'r') as IGBLAST_fmt19:
		result = csv.reader(IGBLAST_fmt19, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				file_handle = open(progressBarFile, 'w')
				progress = str(int((line_id + 1) * 50 / totel_seq))
				file_handle.write(progress)
				file_handle.write(',Step1:' + str(line_id + 1) + '/' + str(totel_seq))
				file_handle.close()

				this_data = [''] * 119

				this_data[0] = record[0]
				this_data[1] = str(len(record[1]))
				if record[2] == 'IGH':
					this_data[2] = 'Heavy'
				elif record[2] == 'IGK':
					this_data[2] = 'Kappa'
				elif record[2] == 'IGL':
					this_data[2] = 'Lambda'
				else:
					this_data[2] = 'Unknown'

				if record[3] == 'F':
					this_data[12] = 'No'
				else:
					this_data[12] = 'Yes'

				if record[4] == '':
					this_data[13] = 'N/A'
				else:
					this_data[13] = record[4]
				if record[5] == '':
					this_data[14] = 'N/A'
				else:
					this_data[14] = record[5]

				if record[6] == 'F':
					this_data[15] = '+'
				else:
					this_data[15] = '-'

				this_data[59] = record[64]
				this_data[60] = record[65]
				this_data[61] = record[68]
				this_data[62] = record[69]
				this_data[63] = ''
				this_data[64] = ''
				this_data[65] = record[72]
				this_data[66] = record[73]

				this_data[67] = record[14]
				this_data[68] = record[15]
				this_data[69] = record[16]
				this_data[70] = record[17]
				this_data[71] = ''
				this_data[72] = ''
				this_data[73] = record[18]
				this_data[74] = record[19]

				# identify grouping
				if project == 'ByFunction':
					this_data[75] = this_data[2]
					if this_data[14] == "Yes":
						this_data[76] = 'Functional'
					else:
						this_data[76] = 'Nonfunctional'
					this_data[77] = this_data[13]
				else:
					this_data[75] = project
					this_data[76] = grouping
					this_data[77] = subgroup

				this_data[78] = species
				this_data[79] = record[10]
				this_data[80] = record[11]

				this_data[86] = "Specificity"
				this_data[87] = "Subspecificity"
				this_data[88] = "0"
				this_data[89] = "0"

				this_data[93] = now

				# identify isotype
				if this_data[2] == 'Heavy':
					IsoSeq = record[1]
					IsoSeq = IsoSeq[int(record[71]):]
					try:
						IsoSeq = IsoSeq.strip('N')
						AGCTs = IsoSeq.count('A') + IsoSeq.count('G') + IsoSeq.count('C') + IsoSeq.count('T')
						if AGCTs > 5:  # todo decide if can determine isotype from < 5 or need more then
							Isotype = VGenesSeq.CallIsotype(IsoSeq)
						else:
							if len(IsoSeq) > 2:
								if IsoSeq[:3] == 'CCT' or IsoSeq == 'CTT':
									Isotype = 'IgG'
								elif IsoSeq[:3] == 'CAT':
									Isotype = 'IgA'
								elif IsoSeq[:3] == 'GGA':
									Isotype = 'IgM'
								elif IsoSeq[:3] == 'CAC':
									Isotype = 'IgD'
								else:
									Isotype = IsoSeq
							else:
								Isotype = 'Unknown'
					except:
						Isotype = 'Unknown'
				else:
					if this_data[2] == 'Kappa':
						Isotype = 'Kappa'
					elif this_data[2] == 'Lambda':
						Isotype = 'Lambda'
				this_data[101] = Isotype

				# import CDR3
				JGeneName = record[9].split(',')[0]
				CDR3_start = int(record[83]) + 1
				if species == 'Human':
					CDR3_end = int(record[70]) + JHuman[JGeneName]
				else:
					CDR3_end = int(record[70]) + JMouse[JGeneName]
				CDR3_NT = record[1]
				CDR3_NT = CDR3_NT[CDR3_start-1:CDR3_end]
				CDR3_AA, msg = Translator(CDR3_NT,0)
				CDR3_len = len(CDR3_AA)
				CDR3_MW = VGenesSeq.OtherParam(CDR3_AA, 'AAMW', 0, True)
				CDR3_pI = VGenesSeq.OtherParam(CDR3_AA, 'AApI', 0, True)

				this_data[81] = CDR3_NT
				this_data[82] = CDR3_AA
				this_data[83] = CDR3_len
				this_data[84] = str(CDR3_start - int(record[62]))
				this_data[85] = str(CDR3_end - int(record[62]))
				this_data[99] = str(CDR3_MW)
				this_data[100] = str(CDR3_pI)

				# import mutation
				mAb_seq = record[10]
				germline_seq = record[11]
				mut, num_mut = IdentifyMutation(mAb_seq, germline_seq)
				this_data[57] = str(num_mut)
				this_data[96] = str(num_mut)
				this_data[97] = mut

				DATA.append(this_data)
			line_id += 1

	# parse IgBlast out fmt 3
	cur_block = ''
	block_id = 0
	read_tag = False
	for IgLine in IgBlastOut_fmt3:
		line_match = re.findall(r'^Query', IgLine)
		if len(line_match) > 0:
			read_tag = True

			if cur_block != '':
				file_handle = open(progressBarFile, 'w')
				progress = str(50 + int((block_id + 1) * 50 / totel_seq))
				file_handle.write(progress)
				file_handle.write(',Step2:' + str(block_id + 1) + '/' + str(totel_seq))
				file_handle.close()

				# import V1,V2,V3,D1,D2,DD3,J1,J2,J3 and V,D,J locus
				ig_match = re.findall('\nIG[^\n]+', cur_block)
				## import V1,V2,V3,D1,D2,DD3,J1,J2,J3
				v_cur_index = 3
				d_cur_index = 6
				j_cur_index = 9
				for line in ig_match[:len(ig_match) - 1]:
					m = re.search('IG\S+', line)
					cur_gene = m.group(0)
					if cur_gene[3] == 'V':
						DATA[block_id][v_cur_index] = cur_gene
						v_cur_index += 1
					elif cur_gene[3] == 'D':
						DATA[block_id][d_cur_index] = cur_gene
						d_cur_index += 1
					else:
						DATA[block_id][j_cur_index] = cur_gene
						j_cur_index += 1
				# V,D,J locus
				m = re.search('IG\S+', DATA[block_id][3])
				try:
					v_locus = m.group(0)
					v_locus = re.sub('^IG', '', v_locus)
					v_locus = re.sub('\*.+', '', v_locus)
					v_locus = v_locus[1] + v_locus[0] + v_locus[2:]
					DATA[block_id][90] = v_locus
				except:
					pass

				m = re.search('IG\S+', DATA[block_id][6])
				try:
					d_locus = m.group(0)
					d_locus = re.sub('^IG', '', d_locus)
					d_locus = re.sub('\*.+', '', d_locus)
					d_locus = d_locus[1] + d_locus[0] + d_locus[2:]
					DATA[block_id][92] = d_locus
				except:
					pass

				m = re.search('IG\S+', DATA[block_id][9])
				try:
					j_locus = m.group(0)
					j_locus = re.sub('^IG', '', j_locus)
					j_locus = re.sub('\*.+', '', j_locus)
					j_locus = j_locus[1] + j_locus[0] + j_locus[2:]
					DATA[block_id][91] = j_locus
				except:
					pass
				
				# import V(D)J junction info
				ig_match = re.findall(r'V-\(D\)-J junction.+\n.+', cur_block)
				junction = ig_match[0]
				junction = junction.split('\n')[1]
				junction_list = junction.split('\t')
				if DATA[block_id][2] == 'Heavy':
					DATA[block_id][16] = junction_list[0]
					DATA[block_id][17] = junction_list[1]
					DATA[block_id][18] = junction_list[2]
					DATA[block_id][19] = junction_list[3]
					DATA[block_id][20] = junction_list[4]
				else:
					DATA[block_id][16] = junction_list[0]
					DATA[block_id][21] = junction_list[1]
					DATA[block_id][20] = junction_list[2]
				
				# import Alignment summary
				ig_match = re.findall(r'\nFR1[^\n]+', cur_block)
				fr1_match = ig_match[0]
				fr1_match = fr1_match[1:]
				fr1_match = fr1_match.split('\t')
				DATA[block_id][22] = 1
				DATA[block_id][23] = str(int(fr1_match[2]) - int(fr1_match[1]) + 1)
				DATA[block_id][24] = fr1_match[3]
				DATA[block_id][25] = fr1_match[4]
				DATA[block_id][26] = fr1_match[5]
				DATA[block_id][27] = fr1_match[6]
				DATA[block_id][28] = fr1_match[7]

				ig_match = re.findall(r'\nCDR1[^\n]+', cur_block)
				cdr1_match = ig_match[0]
				cdr1_match = cdr1_match[1:]
				cdr1_match = cdr1_match.split('\t')
				DATA[block_id][29] = str(int(cdr1_match[1]) - int(fr1_match[1]) + 1)
				DATA[block_id][30] = str(int(cdr1_match[2]) - int(fr1_match[1]) + 1)
				DATA[block_id][31] = cdr1_match[3]
				DATA[block_id][32] = cdr1_match[4]
				DATA[block_id][33] = cdr1_match[5]
				DATA[block_id][34] = cdr1_match[6]
				DATA[block_id][35] = cdr1_match[7]

				ig_match = re.findall(r'\nFR2[^\n]+', cur_block)
				fr2_match = ig_match[0]
				fr2_match = fr2_match[1:]
				fr2_match = fr2_match.split('\t')
				DATA[block_id][36] = str(int(fr2_match[1]) - int(fr1_match[1]) + 1)
				DATA[block_id][37] = str(int(fr2_match[2]) - int(fr1_match[1]) + 1)
				DATA[block_id][38] = fr2_match[3]
				DATA[block_id][39] = fr2_match[4]
				DATA[block_id][40] = fr2_match[5]
				DATA[block_id][41] = fr2_match[6]
				DATA[block_id][42] = fr2_match[7]

				ig_match = re.findall(r'\nCDR2[^\n]+', cur_block)
				cdr2_match = ig_match[0]
				cdr2_match = cdr2_match[1:]
				cdr2_match = cdr2_match.split('\t')
				DATA[block_id][43] = str(int(cdr2_match[1]) - int(fr1_match[1]) + 1)
				DATA[block_id][44] = str(int(cdr2_match[2]) - int(fr1_match[1]) + 1)
				DATA[block_id][45] = cdr2_match[3]
				DATA[block_id][46] = cdr2_match[4]
				DATA[block_id][47] = cdr2_match[5]
				DATA[block_id][48] = cdr2_match[6]
				DATA[block_id][49] = cdr2_match[7]

				ig_match = re.findall(r'\nFR3[^\n]+', cur_block)
				fr3_match = ig_match[0]
				fr3_match = fr3_match[1:]
				fr3_match = fr3_match.split('\t')
				DATA[block_id][50] = str(int(fr3_match[1]) - int(fr1_match[1]) + 1)
				DATA[block_id][51] = str(int(fr3_match[2]) - int(fr1_match[1]) + 1)
				DATA[block_id][52] = fr3_match[3]
				DATA[block_id][53] = fr3_match[4]
				DATA[block_id][54] = fr3_match[5]
				DATA[block_id][55] = fr3_match[6]
				DATA[block_id][56] = fr3_match[7]

				# import Alignment summary
				ig_match = re.findall(r'\nAlignments[\n\S\s]+', cur_block)
				alignment = ig_match[0]
				alignment = alignment[1:]
				alignment = re.sub(r'\n+Lambda[\n\S\s]+','',alignment)
				DATA[block_id][58] = alignment

				block_id += 1

			cur_block = ''

		if read_tag:
			cur_block += IgLine

	if cur_block != '':
		file_handle = open(progressBarFile, 'w')
		progress = str(50 + int((block_id + 1) * 50 / totel_seq))
		file_handle.write(progress)
		file_handle.write(',Step2:' + str(block_id + 1) + '/' + str(totel_seq))
		file_handle.close()

		# import V1,V2,V3,D1,D2,DD3,J1,J2,J3 and V,D,J locus
		ig_match = re.findall('\nIG[^\n]+', cur_block)
		## import V1,V2,V3,D1,D2,DD3,J1,J2,J3
		v_cur_index = 3
		d_cur_index = 6
		j_cur_index = 9
		for line in ig_match[:len(ig_match) - 1]:
			m = re.search('IG\S+', line)
			cur_gene = m.group(0)
			if cur_gene[3] == 'V':
				DATA[block_id][v_cur_index] = cur_gene
				v_cur_index += 1
			elif cur_gene[3] == 'D':
				DATA[block_id][d_cur_index] = cur_gene
				d_cur_index += 1
			else:
				DATA[block_id][j_cur_index] = cur_gene
				j_cur_index += 1
		# V,D,J locus
		m = re.search('IG\S+', DATA[block_id][3])
		try:
			v_locus = m.group(0)
			v_locus = re.sub('^IG', '', v_locus)
			v_locus = re.sub('\*.+', '', v_locus)
			v_locus = v_locus[1] + v_locus[0] + v_locus[2:]
			DATA[block_id][90] = v_locus
		except:
			pass

		m = re.search('IG\S+', DATA[block_id][6])
		try:
			d_locus = m.group(0)
			d_locus = re.sub('^IG', '', d_locus)
			d_locus = re.sub('\*.+', '', d_locus)
			d_locus = d_locus[1] + d_locus[0] + d_locus[2:]
			DATA[block_id][92] = d_locus
		except:
			pass

		m = re.search('IG\S+', DATA[block_id][9])
		try:
			j_locus = m.group(0)
			j_locus = re.sub('^IG', '', j_locus)
			j_locus = re.sub('\*.+', '', j_locus)
			j_locus = j_locus[1] + j_locus[0] + j_locus[2:]
			DATA[block_id][91] = j_locus
		except:
			pass

		# import V(D)J junction info
		ig_match = re.findall(r'V-\(D\)-J junction.+\n.+', cur_block)
		junction = ig_match[0]
		junction = junction.split('\n')[1]
		junction_list = junction.split('\t')
		if DATA[block_id][2] == 'Heavy':
			DATA[block_id][16] = junction_list[0]
			DATA[block_id][17] = junction_list[1]
			DATA[block_id][18] = junction_list[2]
			DATA[block_id][19] = junction_list[3]
			DATA[block_id][20] = junction_list[4]
		else:
			DATA[block_id][16] = junction_list[0]
			DATA[block_id][21] = junction_list[1]
			DATA[block_id][20] = junction_list[2]

		# import Alignment summary
		ig_match = re.findall(r'\nFR1[^\n]+', cur_block)
		fr1_match = ig_match[0]
		fr1_match = fr1_match[1:]
		fr1_match = fr1_match.split('\t')
		DATA[block_id][22] = 1
		DATA[block_id][23] = str(int(fr1_match[2]) - int(fr1_match[1]) + 1)
		DATA[block_id][24] = fr1_match[3]
		DATA[block_id][25] = fr1_match[4]
		DATA[block_id][26] = fr1_match[5]
		DATA[block_id][27] = fr1_match[6]
		DATA[block_id][28] = fr1_match[7]

		ig_match = re.findall(r'\nCDR1[^\n]+', cur_block)
		cdr1_match = ig_match[0]
		cdr1_match = cdr1_match[1:]
		cdr1_match = cdr1_match.split('\t')
		DATA[block_id][29] = str(int(cdr1_match[1]) - int(fr1_match[1]) + 1)
		DATA[block_id][30] = str(int(cdr1_match[2]) - int(fr1_match[1]) + 1)
		DATA[block_id][31] = cdr1_match[3]
		DATA[block_id][32] = cdr1_match[4]
		DATA[block_id][33] = cdr1_match[5]
		DATA[block_id][34] = cdr1_match[6]
		DATA[block_id][35] = cdr1_match[7]

		ig_match = re.findall(r'\nFR2[^\n]+', cur_block)
		fr2_match = ig_match[0]
		fr2_match = fr2_match[1:]
		fr2_match = fr2_match.split('\t')
		DATA[block_id][36] = str(int(fr2_match[1]) - int(fr1_match[1]) + 1)
		DATA[block_id][37] = str(int(fr2_match[2]) - int(fr1_match[1]) + 1)
		DATA[block_id][38] = fr2_match[3]
		DATA[block_id][39] = fr2_match[4]
		DATA[block_id][40] = fr2_match[5]
		DATA[block_id][41] = fr2_match[6]
		DATA[block_id][42] = fr2_match[7]

		ig_match = re.findall(r'\nCDR2[^\n]+', cur_block)
		cdr2_match = ig_match[0]
		cdr2_match = cdr2_match[1:]
		cdr2_match = cdr2_match.split('\t')
		DATA[block_id][43] = str(int(cdr2_match[1]) - int(fr1_match[1]) + 1)
		DATA[block_id][44] = str(int(cdr2_match[2]) - int(fr1_match[1]) + 1)
		DATA[block_id][45] = cdr2_match[3]
		DATA[block_id][46] = cdr2_match[4]
		DATA[block_id][47] = cdr2_match[5]
		DATA[block_id][48] = cdr2_match[6]
		DATA[block_id][49] = cdr2_match[7]

		ig_match = re.findall(r'\nFR3[^\n]+', cur_block)
		fr3_match = ig_match[0]
		fr3_match = fr3_match[1:]
		fr3_match = fr3_match.split('\t')
		DATA[block_id][50] = str(int(fr3_match[1]) - int(fr1_match[1]) + 1)
		DATA[block_id][51] = str(int(fr3_match[2]) - int(fr1_match[1]) + 1)
		DATA[block_id][52] = fr3_match[3]
		DATA[block_id][53] = fr3_match[4]
		DATA[block_id][54] = fr3_match[5]
		DATA[block_id][55] = fr3_match[6]
		DATA[block_id][56] = fr3_match[7]

		# import Alignment summary
		ig_match = re.findall(r'\nAlignments[\n\S\s]+', cur_block)
		alignment = ig_match[0]
		alignment = alignment[1:]
		alignment = re.sub(r'\n+Lambda[\n\S\s]+', '', alignment)
		DATA[block_id][58] = alignment

		block_id += 1

	ErLog = '\nVGenes input ended at: ' + time.strftime('%c')
	with open(ErlogFile2, 'a') as currentFile:  # using with for this automatically closes the file even if you crash
		currentFile.write(ErLog)

	# if productive = TRUE, only keep record with V and J
	start = time.time()

	if GetProductive == 0:  # only keep productive
		# find record
		del_index = []
		for index in range(len(DATA)):
			if DATA[index][14] == 'Yes':
				del_index.append(index)

		for i in del_index:
			ErLog = DATA[i][0] + ' was not a productive rearrangement\n'
			with open(ErlogFile, 'a') as currentfile:
				currentfile.write(ErLog)

		# delete record
		cnt = 0
		for index in del_index:
			del DATA[index - cnt]
			cnt += 1
	elif GetProductive == 1:    #only keep V and J
		# find record
		del_index = []
		for index in range(len(DATA)):
			if DATA[index][90] == '' or DATA[index][91] == '':
				del_index.append(index)

		for i in del_index:
			ErLog = DATA[i][0] + ' missed V or J\n'
			with open(ErlogFile, 'a') as currentfile:
				currentfile.write(ErLog)
		# delete record
		cnt = 0
		for index in del_index:
			del DATA[index - cnt]
			cnt += 1

	end = time.time()
	print('Run time for productive: ' + str(end - start))

	return DATA

def IMGTparser(IMGT_out, data_list):
	progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')

	DBpathname = os.path.join(working_prefix, 'Data', 'VDJGenes.db')
	conn = db.connect(DBpathname)
	cursor = conn.cursor()

	# Extract IMGT files from txz package
	temp_dir, files = extractIMGT(IMGT_out)

	project = data_list[0]
	grouping = data_list[1]
	subgroup = data_list[2]

	now = time.strftime('%c')
	DATA = []
	raw_seq = []
	# Parse
	with open(files['summary'], 'r') as summary, \
			open(files['gapped-nt'], 'r') as gapped, \
			open(files['ntseq'], 'r') as ntseq, \
			open(files['junction'], 'r') as junction, \
			open(files['para'], 'r') as para, \
			open(files['8V-mut'], 'r') as vmut, \
			open(files['aaseq'], 'r') as aaseq, \
			open(files['junction'], 'r') as junction:

		# read species
		result = para.readlines()
		spe = result[3]
		spe = re.sub('Species\t','',spe)
		spe = re.sub('\t', '', spe)
		spe = re.sub('\n', '', spe)
		
		if spe == "Homo sapiens":
			spe = 'Human'
		elif spe == '':
			spe = 'Mouse'
		else:
			Msg = 'Your species\n' + spe + 'is not supported! Only support Human and Mouse!'
			return Msg

		# read records from summary file
		result = csv.reader(summary, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				this_data = [''] * 119

				this_data[0] = record[1]
				this_data[1] = record[27]
				this_data[13] = record[19]
				if record[2] == 'productive':
					this_data[14] = 'Yes'
				else:
					this_data[14] = 'No'
				this_data[15] = record[20]
				#this_data[79] = record[24]
				this_data[78] = spe
				this_data[93] = now

				v_str = record[3]
				v_genes = v_str.split('or ')
				index = 3
				for ele in v_genes:
					m = re.search('IG\S+', ele)
					try:
						this_data[index] = m.group(0)
					except:
						this_data[index] = ele
					index += 1

				d_str = record[11]
				d_genes = d_str.split('or ')
				index = 6
				for ele in d_genes:
					m = re.search('IG\S+', ele)
					try:
						this_data[index] = m.group(0)
					except:
						this_data[index] = ele
					index += 1

				j_str = record[7]
				j_genes = j_str.split('or ')
				index = 9
				for ele in j_genes:
					m = re.search('IG\S+', ele)
					try:
						this_data[index] = m.group(0)
					except:
						this_data[index] = ele
					index += 1

				m = re.search('IG\S+', this_data[3])
				try:
					v_locus = m.group(0)
					v_locus = re.sub('^IG','',v_locus)
					v_locus = re.sub('\*.+', '', v_locus)
					v_locus = v_locus[1] + v_locus[0] + v_locus[2:]
					this_data[90] = v_locus
				except:
					pass

				m = re.search('IG\S+', this_data[6])
				try:
					d_locus = m.group(0)
					d_locus = re.sub('^IG', '', d_locus)
					d_locus = re.sub('\*.+', '', d_locus)
					d_locus = d_locus[1] + d_locus[0] + d_locus[2:]
					this_data[92] = d_locus
				except:
					pass

				m = re.search('IG\S+', this_data[9])
				try:
					j_locus = m.group(0)
					j_locus = re.sub('^IG', '', j_locus)
					j_locus = re.sub('\*.+', '', j_locus)
					j_locus = j_locus[1] + j_locus[0] + j_locus[2:]
					this_data[91] = j_locus
				except:
					pass

				# identify gene type
				GeneType = 'Lambda'
				if v_locus[:2] == 'VH':
					GeneType = 'Heavy'
				elif v_locus[:2] == 'VK':
					GeneType = 'Kappa'
				this_data[2] = GeneType

				# identify grouping
				if project == 'ByFunction':
					this_data[75] = GeneType
					if this_data[14] == "Yes":
						this_data[76] = 'Functional'
					else:
						this_data[76] = 'Nonfunctional'
					this_data[77] = this_data[13]
				else:
					this_data[75] = project
					this_data[76] = grouping
					this_data[77] = subgroup

				# clone and rank
				this_data[86] = 'Specificity'
				this_data[87] = 'Subspecificity'
				this_data[88] = '0'
				this_data[89] = '0'

				DATA.append(this_data)
				raw_seq.append(record[24].upper())
			line_id += 1

		file_handle = open(progressBarFile, 'w')
		progress = str(int(20))
		file_handle.write(progress)
		file_handle.write(',Summary file processed')
		file_handle.close()

		# read records from junction file
		result = csv.reader(junction, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				if len(record) > 73:
					DATA[line_id - 1][99] = record[71]
					DATA[line_id - 1][100] = record[72]

					DATA[line_id - 1][16] = record[9].upper()
					DATA[line_id - 1][17] = record[12].upper()
					DATA[line_id - 1][18] = record[14].upper()
					DATA[line_id - 1][19] = record[19].upper()
					DATA[line_id - 1][20] = record[29].upper()
					DATA[line_id - 1][21] = record[11].upper()

					# make Germline sequence
					if record[6] == 'in-frame':
						if spe == 'Human':
							SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + DATA[line_id - 1][3] + '"'
							SqlStatementD = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + DATA[line_id - 1][6] + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + DATA[line_id - 1][9] + '"'
						elif spe == 'Mouse':
							SqlStatementV = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + DATA[line_id - 1][3] + '"'
							SqlStatementD = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + DATA[line_id - 1][6] + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + DATA[line_id - 1][9] + '"'

						if DATA[line_id - 1][2] == 'Heavy':
							# V gene
							Germline_Vseq = []
							cursor.execute(SqlStatementV)
							for row in cursor:
								for column in row:
									Germline_Vseq.append(column)

							# D gene
							Germline_Dseq = []
							cursor.execute(SqlStatementD)
							for row in cursor:
								for column in row:
									Germline_Dseq.append(column)

							# J gene
							Germline_Jseq = []
							cursor.execute(SqlStatementJ)
							for row in cursor:
								for column in row:
									Germline_Jseq.append(column)

							Germline_seq = Germline_Vseq[4] + record[12].upper() + Germline_Dseq[4] + record[19].upper() + Germline_Jseq[4]  # V + D1 + D + D2 + J

							Germline_Vbeg = 1
							Germline_Vend = len(Germline_Vseq[4])
							offset = Germline_Vend
							if record[12] != '':
								Germline_D1beg = offset + 1
								Germline_D1end = offset + len(record[12])
								offset = Germline_D1end
							else:
								Germline_D1beg = 0
								Germline_D1end = 0

							if Germline_Dseq[4] != '':
								offset = offset + len(Germline_Dseq[4])

							if record[19] != '':
								Germline_D2beg = offset + 1
								Germline_D2end = offset + len(record[19])
								offset = Germline_D2end
							else:
								Germline_D2beg = ''
								Germline_D2end = ''

							Germline_Jbeg = offset + 1
							Germline_Jend = offset + len(Germline_Jseq[4])

							DATA[line_id - 1][80] = Germline_seq
							DATA[line_id - 1][59] = str(Germline_Vbeg)
							DATA[line_id - 1][60] = str(Germline_Vend)
							DATA[line_id - 1][61] = str(Germline_D1beg)
							DATA[line_id - 1][62] = str(Germline_D1end)
							DATA[line_id - 1][63] = str(Germline_D2beg)
							DATA[line_id - 1][64] = str(Germline_D2end)
							DATA[line_id - 1][65] = str(Germline_Jbeg)
							DATA[line_id - 1][66] = str(Germline_Jend)
						else:
							# V gene
							Germline_Vseq = []
							cursor.execute(SqlStatementV)
							for row in cursor:
								for column in row:
									Germline_Vseq.append(column)

							# J gene
							Germline_Jseq = []
							cursor.execute(SqlStatementJ)
							for row in cursor:
								for column in row:
									Germline_Jseq.append(column)

							Germline_seq = Germline_Vseq[4] + record[11].upper() + Germline_Jseq[4]  # V + N + J

							Germline_Vbeg = 1
							Germline_Vend = len(Germline_Vseq[4])
							if record[11] != '':
								offset = Germline_Vend + len(record[11])
							else:
								offset = Germline_Vend
							Germline_D1beg = 0
							Germline_D1end = 0
							Germline_D2beg = ''
							Germline_D2end = ''

							Germline_Jbeg = offset + 1
							Germline_Jend = offset + len(Germline_Jseq[4])

							DATA[line_id - 1][80] = Germline_seq
							DATA[line_id - 1][59] = str(Germline_Vbeg)
							DATA[line_id - 1][60] = str(Germline_Vend)
							DATA[line_id - 1][61] = str(Germline_D1beg)
							DATA[line_id - 1][62] = str(Germline_D1end)
							DATA[line_id - 1][63] = str(Germline_D2beg)
							DATA[line_id - 1][64] = str(Germline_D2end)
							DATA[line_id - 1][65] = str(Germline_Jbeg)
							DATA[line_id - 1][66] = str(Germline_Jend)
			line_id += 1

		file_handle = open(progressBarFile, 'w')
		progress = str(int(40))
		file_handle.write(progress)
		file_handle.write(',Junction file processed')
		file_handle.close()

		# read records from ntseq file
		result = csv.reader(ntseq, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				DATA[line_id-1][22] = record[48]
				DATA[line_id-1][23] = record[49]
				DATA[line_id-1][24] = str(int(record[49]) - int(record[48]) + 1)
				DATA[line_id-1][29] = record[50]
				DATA[line_id-1][30] = record[51]
				DATA[line_id-1][31] = str(int(record[51]) - int(record[50]) + 1)
				DATA[line_id-1][36] = record[52]
				DATA[line_id-1][37] = record[53]
				DATA[line_id-1][38] = str(int(record[53]) - int(record[52]) + 1)
				DATA[line_id-1][43] = record[54]
				DATA[line_id-1][44] = record[55]
				DATA[line_id-1][45] = str(int(record[55]) - int(record[54]) + 1)
				DATA[line_id-1][50] = record[56]
				DATA[line_id-1][51] = record[57]
				DATA[line_id-1][52] = str(int(record[57]) - int(record[56]) + 1)

				DATA[line_id-1][67] = record[46]
				DATA[line_id-1][68] = record[47]
				DATA[line_id-1][71] = record[90]
				DATA[line_id-1][72] = record[91]
				DATA[line_id-1][73] = record[110]
				DATA[line_id-1][74] = record[111]
				if DATA[line_id-1][2] == 'Heavy':
					DATA[line_id-1][79] = record[6].upper()
					DATA[line_id - 1][69] = record[76]
					DATA[line_id - 1][70] = record[77]
				else:
					DATA[line_id - 1][79] = record[7].upper()
					DATA[line_id - 1][69] = record[82]
					DATA[line_id - 1][70] = record[83]
				DATA[line_id-1][81] = record[14].upper()
				DATA[line_id-1][84] = record[58]
				DATA[line_id-1][85] = record[59]

				GeneType = DATA[line_id-1][2]
				Sequence = DATA[line_id-1][79]

				aa_seq,_msg = Translator(DATA[line_id - 1][79], 0)
				if '*' in aa_seq:
					DATA[line_id - 1][12] = 'Yes'
				else:
					DATA[line_id - 1][12] = 'No'

				# identify isotype
				if GeneType == 'Heavy':
					IsoSeq = raw_seq[line_id-1].split(Sequence)
					try:
						IsoSeq = IsoSeq[1]
						# print(SeqName)
						IsoSeq = IsoSeq.strip('N')
						AGCTs = IsoSeq.count('A') + IsoSeq.count('G') + IsoSeq.count('C') + IsoSeq.count('T')
						if AGCTs > 5:  # todo decide if can determine isotype from < 5 or need more then
							Isotype = VGenesSeq.CallIsotype(IsoSeq)
						else:
							if len(IsoSeq) > 2:
								if IsoSeq[:3] == 'CCT' or IsoSeq == 'CTT':
									Isotype = 'IgG'
								elif IsoSeq[:3] == 'CAT':
									Isotype = 'IgA'
								elif IsoSeq[:3] == 'GGA':
									Isotype = 'IgM'
								elif IsoSeq[:3] == 'CAC':
									Isotype = 'IgD'
								else:
									Isotype = IsoSeq
							else:
								Isotype = 'Unknown'
					except:
						Isotype = 'Unknown'
				else:
					if GeneType == 'Kappa':
						Isotype = 'Kappa'
					elif GeneType == 'Lambda':
						Isotype = 'Lambda'
				DATA[line_id - 1][101] = Isotype

			line_id += 1

		file_handle = open(progressBarFile, 'w')
		progress = str(int(60))
		file_handle.write(progress)
		file_handle.write(',NT-sequence file processed')
		file_handle.close()

		# read records from 8-v-mutattion file
		result = csv.reader(vmut, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				record[22] = re.sub('\s.+', '', record[22])
				record[23] = re.sub('\s.+', '', record[23])
				record[24] = re.sub('\s.+', '', record[24])
				record[40] = re.sub('\s.+', '', record[40])
				record[41] = re.sub('\s.+', '', record[41])
				record[42] = re.sub('\s.+', '', record[42])
				record[58] = re.sub('\s.+', '', record[58])
				record[59] = re.sub('\s.+', '', record[59])
				record[60] = re.sub('\s.+', '', record[60])
				record[76] = re.sub('\s.+', '', record[76])
				record[77] = re.sub('\s.+', '', record[77])
				record[78] = re.sub('\s.+', '', record[78])
				record[94] = re.sub('\s.+', '', record[94])
				record[95] = re.sub('\s.+', '', record[95])
				record[96] = re.sub('\s.+', '', record[96])

				DATA[line_id-1][25] = record[24]
				DATA[line_id-1][26] = record[25]
				DATA[line_id-1][27] = str(int(record[22]) - int(record[23]))
				DATA[line_id-1][28] = str(round(int(record[24])/int(record[22])*100,2))
				DATA[line_id-1][32] = record[42]
				DATA[line_id-1][33] = record[43]
				DATA[line_id-1][34] = str(int(record[40]) - int(record[41]))
				DATA[line_id-1][35] = str(round(int(record[42]) / int(record[40]) * 100, 2))
				DATA[line_id-1][39] = record[60]
				DATA[line_id-1][40] = record[61]
				DATA[line_id-1][41] = str(int(record[58]) - int(record[59]))
				DATA[line_id-1][42] = str(round(int(record[60]) / int(record[58]) * 100, 2))
				DATA[line_id-1][46] = record[78]
				DATA[line_id-1][47] = record[79]
				DATA[line_id-1][48] = str(int(record[76]) - int(record[77]))
				DATA[line_id-1][49] = str(round(int(record[78]) / int(record[76]) * 100, 2))
				DATA[line_id-1][53] = record[96]
				DATA[line_id-1][54] = record[97]
				DATA[line_id-1][55] = str(int(record[94]) - int(record[95]))
				DATA[line_id-1][56] = str(round(int(record[96]) / int(record[94]) * 100, 2))
				DATA[line_id-1][57] = re.sub('\s.+','',record[7])
				DATA[line_id-1][96] = re.sub('\s.+','',record[7])
			line_id += 1

		file_handle = open(progressBarFile, 'w')
		progress = str(int(80))
		file_handle.write(progress)
		file_handle.write(',Mutation file processed')
		file_handle.close()

		# read records from aa-seq file
		result = csv.reader(aaseq, delimiter="\t")
		line_id = 0
		for record in result:
			if line_id == 0:
				pass
			else:
				DATA[line_id-1][82] = record[14]
				DATA[line_id-1][83] = str(len(record[14]))
			line_id += 1

		file_handle = open(progressBarFile, 'w')
		progress = str(int(100))
		file_handle.write(progress)
		file_handle.write(',AA-sequence file processed')
		file_handle.close()
	#result = IMGTReader(summary, gapped, ntseq, junction, receptor=False)
		#for record in result:
		#	#print(record)
		#	DATA.append(record)
	# Remove IMGT temporary directory
	temp_dir.cleanup()

	return DATA

def IdentifyMutation(mAb_seq, germline_seq):
	mut_list = []
	if len(mAb_seq) == len(germline_seq):
		for i in range(len(mAb_seq)):
			pos = i + 1
			if mAb_seq[i] != germline_seq[i]:
				if mAb_seq[i] in 'ATCG' and germline_seq[i] in 'ATCG':
					cur_mut = germline_seq[i] + '-' + str(pos) + '-' + mAb_seq[i]
					mut_list.append(cur_mut)
	num_mut = len(mut_list)
	mut = ','.join(mut_list)
	return mut, num_mut

def extractIMGT(imgt_output):
    """
    Extract necessary files from IMGT/HighV-QUEST results. This function is imported from change-o.

    Arguments:
      imgt_output : zipped file or unzipped folder output by IMGT/HighV-QUEST.

    Returns:
      tuple : (temporary directory handle, dictionary with names of extracted IMGT files).
    """
    # Map of IMGT file names
    imgt_names = ('11_Parameters', '1_Summary', '2_IMGT-gapped', '3_Nt-sequences', '4_IMGT-gapped', '5_AA', '6_Junction',
                  '7_V-REGION-mutation', '8_V-REGION-nt-mutation', '9_V-REGION-AA-change')
    imgt_keys = ('para', 'summary', 'gapped-nt', 'ntseq', 'gapped-aa', 'aaseq', 'junction',
                 '7-v-mutation-aa-table', '8V-mut','9-v-mutation-aa')

    # Open temporary directory and intialize return dictionary
    temp_dir = TemporaryDirectory()

    # Zip input
    if zipfile.is_zipfile(imgt_output):
        imgt_zip = zipfile.ZipFile(imgt_output, 'r')
        # Extract required files
        imgt_files = sorted([n for n in imgt_zip.namelist() \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_zip.extractall(temp_dir.name, imgt_files)
        # Define file dictionary
        imgt_dict = {k: os.path.join(temp_dir.name, f) for k, f in zip_longest(imgt_keys, imgt_files)}
    # Folder input
    elif os.path.isdir(imgt_output):
        folder_files = []
        for root, dirs, files in os.walk(imgt_output):
            folder_files.extend([os.path.join(os.path.abspath(root), f) for f in files])
        # Define file dictionary
        imgt_files = sorted([n for n in folder_files \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_dict = {k: f for k, f in zip_longest(imgt_keys, imgt_files)}
    # Tarball input
    elif tarfile.is_tarfile(imgt_output):
        imgt_tar = tarfile.open(imgt_output, 'r')
        # Extract required files
        imgt_files = sorted([n for n in imgt_tar.getnames() \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_tar.extractall(temp_dir.name, [imgt_tar.getmember(n) for n in imgt_files])
        # Define file dictionary
        imgt_dict = {k: os.path.join(temp_dir.name, f) for k, f in zip_longest(imgt_keys, imgt_files)}
    else:
        printError('Unsupported IGMT output file. Must be either a zipped file (.zip), LZMA compressed tarfile (.txz) or a folder.')

    # Check extraction for errors
    if len(imgt_dict) != len(imgt_names):
        printError('Extra files or missing necessary file IMGT output %s.' % imgt_output)

    return temp_dir, imgt_dict

def AlignSequencesHTML(DataSet, template):
	# import tempfile
	import os
	global GLMsg
	global working_prefix
	global clustal_path
	global temp_folder
	global VGenesTextWindows
	global muscle_path

	# align selected sequences (AA) using muscle
	all = dict()
	time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())
	outfilename = os.path.join(temp_folder, "out-" + time_stamp + ".fas")
	aafilename = os.path.join(temp_folder, "in-" + time_stamp + ".fas")
	if len(DataSet) == 1:
		SeqName = DataSet[0][0].replace('\n', '').replace('\r', '')
		SeqName = SeqName.strip()
		NTseq = DataSet[0][1]
		AAseq, ErMessage = VGenesSeq.Translator(NTseq, 0)
		all[SeqName] = [NTseq, AAseq]

		out_handle = open(outfilename,'w')
		out_handle.write('>' + SeqName + '\n')
		out_handle.write(AAseq)
		out_handle.close()
	else:
		aa_handle = open(aafilename,'w')
		for record in DataSet:
			SeqName = record[0].replace('\n', '').replace('\r', '')
			SeqName = SeqName.strip()
			NTseq = record[1]
			# sequence check for NT seq
			pattern = re.compile(r'[^ATCGUatcgu]')
			cur_strange = pattern.findall(NTseq)
			cur_strange = list(set(cur_strange))
			if len(cur_strange) > 0:
				ErrMsg = "We find Unlawful nucleotide: " + ','.join(cur_strange) + '\nfrom \n' + SeqName + \
				         '\nPlease remove those Unlawful nucleotide!'
				return ErrMsg

			AAseq, ErMessage = VGenesSeq.Translator(NTseq, 0)
			AAseq = AAseq.replace('*','X').replace('~','Z').replace('.','J')
			all[SeqName] = [NTseq, AAseq]
			aa_handle.write('>' + SeqName + '\n')
			aa_handle.write(AAseq + '\n')
		aa_handle.close()

		cmd = muscle_path
		cmd += " -in " + aafilename + " -out " + outfilename
		try:
			os.system(cmd)
		except:
			QMessageBox.warning(self, 'Warning', 'Fail to run muscle! Check your muscle path!', QMessageBox.Ok,
			                    QMessageBox.Ok)
			return

	# read alignment file, make alignment NT and AA sequences
	SeqName = ''
	AAseq = ''
	if os.path.isfile(outfilename):
		currentfile = open(outfilename, 'r')
		lines = currentfile.readlines()
		for line in lines:
			Readline = line.replace('\n', '').replace('\r', '')
			Readline = Readline.strip()
			if Readline[0] == '>':
				if SeqName != '':
					AAseq, NTseq = BuildNTalignment(AAseq, all[SeqName][0])
					all[SeqName] = [NTseq, AAseq]
				SeqName = Readline[1:]
				AAseq = ''
			else:
				AAseq += Readline
		AAseq, NTseq = BuildNTalignment(AAseq, all[SeqName][0])
		all[SeqName] = [NTseq, AAseq]
	else:
		return

	#if os.path.exists(outfilename):
	#	os.remove(outfilename)
	#if os.path.exists(aafilename):
	#	os.remove(aafilename)

	# generate consnesus sequences (AA and NT)
	if len(all) == 1:
		for key in all:
			consensusDNA = all[key][0]
			consensusAA = all[key][1]
	else:
		firstOne = all[SeqName]
		seqlen = len(firstOne[0])

		consensusDNA = ''
		tester = ''

		for i in range(seqlen):
			tester = ''
			Cnuc = ''
			for key in all:
				try:
					seq = all[key][0]
					tester += seq[i]
				except:
					Msg = 'Find sequence error in ' + key + ', please check your sequence!'
					QMessageBox.warning(self, 'Warning', Msg, QMessageBox.Ok, QMessageBox.Ok)
					return

			frequencies = [(c, tester.count(c)) for c in set(tester)]
			Cnuc = max(frequencies, key=lambda x: x[1])[0]
			consensusDNA += Cnuc


		consensusAA = ''
		firstOne = all[SeqName]
		seqlen = len(firstOne[1])
		for i in range(seqlen):
			tester = ''
			Caa = ''
			for key in all:
				seq = all[key][1]
				tester += seq[i]

			frequencies = [(c, tester.count(c)) for c in set(tester)]
			Caa = max(frequencies, key=lambda x: x[1])[0]
			consensusAA += Caa

	# align consensus AA sequence with template to generate H1 and H3 numbering
	compact_consensusAA = consensusAA.replace(' ', '')

	# make header HTML
	pos_aa_data = [list(range(1,len(compact_consensusAA)+1)),list(range(1,len(compact_consensusAA)+1))]
	div_pos_aa = MakeDivPosAA('line line_pos_aa', 'Position AA:', 'Original AA position: ', pos_aa_data)
	div_con_aa = MakeDivAA('line con_aa', 'Template AA:', compact_consensusAA)
	pos_nt_data = [list(range(1, len(consensusDNA) + 1)), list(range(1, len(consensusDNA) + 1))]
	div_pos_nt = MakeDivPosNT('line line_pos_nt', 'Position NT:', 'Original NT position: ', pos_nt_data)
	div_con_nt = MakeDivNT('line con_nt', 'Template NT:', consensusDNA)

	# initial and open HTML file
	time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())
	out_html_file = os.path.join(temp_folder, time_stamp + '.html')
	if template == '':
		header_file = os.path.join(working_prefix, 'Data', 'template.html')
	else:
		header_file = os.path.join(working_prefix, 'Data', template + '.html')
	shutil.copyfile(header_file, out_html_file)
	out_file_handle = open(out_html_file, 'a')

	JSdata = '<script type="text/javascript">\n'
	JSdata += 'var data = {\n'
	JSarray = []
	JStext = '"Seq0":["Consensus","' + compact_consensusAA + '","' + consensusDNA + '"]'
	JSarray.append(JStext)

	Optiondata = '<script type="text/javascript">\n'
	Optiondata += '$("#option").append("<option value =\'Seq0\'>Consensus Sequence</option>");\n'

	name_div = '<div class="name_div">\n'
	seq_div = '<div class = "seq_div">\n'
	# write header section
	name_div += div_pos_aa[0] + '\n'
	seq_div += div_pos_aa[1] + '\n'
	name_div += div_con_aa[0] + '\n'
	seq_div += div_con_aa[1] + '\n'
	name_div += div_pos_nt[0] + '\n'
	seq_div += div_pos_nt[1] + '\n'
	name_div += div_con_nt[0] + '\n'
	seq_div += div_con_nt[1] + '\n'
	# make sequence section HTML
	i = 1
	for key in all:
		seq_nick_name = 'Seq' + str(i)

		seq_nt = all[key][0]
		seq_aa = all[key][1]
		con_nt = MakeConSeq(seq_nt, consensusDNA)
		con_aa = MakeConSeq(seq_aa, compact_consensusAA)

		div_aa = MakeDivAA('line line_aa ' + seq_nick_name, key, seq_aa)
		div_aa_mut = MakeDivAA('line line_con_aa ' + seq_nick_name, key, con_aa)
		div_nt = MakeDivNT('line line_nt ' + seq_nick_name, key, seq_nt)
		div_nt_mut = MakeDivNT('line line_con_nt ' + seq_nick_name, key, con_nt)
		# write sequence section
		name_div += div_aa[0] + '\n'
		seq_div += div_aa[1] + '\n'
		name_div += div_aa_mut[0] + '\n'
		seq_div += div_aa_mut[1] + '\n'
		name_div += div_nt[0] + '\n'
		seq_div += div_nt[1] + '\n'
		name_div += div_nt_mut[0] + '\n'
		seq_div += div_nt_mut[1] + '\n'

		JStext = '"' + seq_nick_name + '":["' + key + '","' + seq_aa + '","' + seq_nt + '"]'
		JSarray.append(JStext)
		Optiondata += '$("#option").append("<option value =\'' + seq_nick_name + '\'>' + key + '</option>");\n'
		i += 1

	JSdata += ',\n'.join(JSarray)
	JSdata += '\n}\n</script>\n'
	Optiondata += '</script>\n'

	name_div += '</div>\n'
	seq_div += '</div>\n'

	out_file_handle.write(JSdata)
	out_file_handle.write(Optiondata)
	out_file_handle.write('<div class="box">')
	out_file_handle.write(name_div)
	out_file_handle.write(seq_div)
	out_file_handle.write('\n</div>\n</body>\n</html>')
	out_file_handle.close()
	return out_html_file

def MakeConSeq(seq, con):
	for i in range(len(con)):
		if seq[i] == con[i]:
			seq = seq[:i] + '.' + seq[i+1:]
	return seq

def SparseSeq(seq):
	tmp = list(seq)
	seq = ' ' + '  '.join(tmp) + ' '
	return seq

def BuildNTalignment(aa, nt):
	pos = 0
	new_nt = ''
	for i in range(len(aa)):
		cur_aa = aa[i]
		if cur_aa == '-':
			new_nt += '---'
		elif cur_aa == 'X':
			new_nt += nt[pos:pos + 3]
			aa = aa[:i] + '*' + aa[i+1:]
			pos = pos + 3
		elif cur_aa == 'Z':
			new_nt += nt[pos:] + '-'*(3 - len(nt[pos:]))
			aa = aa[:i] + '~' + aa[i + 1:]
			pos = pos + 3
		elif cur_aa == 'J':
			new_nt += nt[pos:pos + 3]
			pos = pos + 3
		else:
			new_nt += nt[pos:pos + 3]
			pos = pos + 3
	a = 1
	return aa, new_nt

def SequenceCheck(sequence, type):
	Msg = 'none'
	if type == 'aa':
		pattern = re.compile(r'[^ILVFMCAGPTSYWQNHEDKR]')
	else:
		pattern = re.compile(r'[^ATCUG]')

	strange_residues = re.findall(pattern, sequence)

	if len(strange_residues) > 0:
		Msg = ','.join(strange_residues)

	return Msg

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

def AA2NT(sequence, dic):
	nt_seq = ''
	for i in range(len(sequence)):
		nt_seq += dic[sequence[i]]

	return nt_seq

def checkOverlap(x1,y1,x2,y2):
	a = range(x1,y1+1)
	b = range(x2,y2+1)
	inter = set(a).intersection(set(b))

	if len(inter) > 0:
		return True
	else:
		return False

def MakeDivNT(class_name, line_name, data):
	div_name = 	'<div class="' + class_name + ' 1">'
	div_name += '<span class="name">' + line_name + '<span class ="name_tip">' +  line_name + '</span></span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + ' 2">'
	count = 0
	for i in range(len(data)):
		if count == 0:
			div_seq += '<span class="unit_pack">'
		elif count%3 == 0:
			div_seq += '</span><span class="unit_pack">'
		div_seq += '<span class="unit">' + data[i] + '</span>'
		count += 1
	div_seq += '</span>'
	div_seq += '</div>'

	return div_name, div_seq

def MakeDivNTDonor(class_name, line_name, data, ori_seq, donor_region):
	div_name = 	'<div class="' + class_name + ' 1">'
	div_name += '<span class="name">' + line_name + '<span class ="name_tip">' +  line_name + '</span></span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + ' 2">'
	cur_pos = 1
	for i in range(len(data)):
		if i == 0:
			count = int(i/3)
			if ori_seq[count] == "-":
				div_seq += '<span class="unit_pack">'
			else:
				if cur_pos in donor_region:
					div_seq += '<span class="unit_pack donor">'
				else:
					div_seq += '<span class="unit_pack">'
				cur_pos += 1
		elif i%3 == 0:
			count = int(i / 3)
			if ori_seq[count] == "-":
				div_seq += '</span><span class="unit_pack">'
			else:
				if cur_pos in donor_region:
					div_seq += '</span><span class="unit_pack donor">'
				else:
					div_seq += '</span><span class="unit_pack">'
				cur_pos += 1

		div_seq += '<span class="unit">' + data[i] + '</span>'
	div_seq += '</span>'
	div_seq += '</div>'

	return div_name, div_seq

def MakeDivAA(class_name, line_name, data):
	div_name = '<div class="' + class_name + ' 1">'
	div_name += '<span class="name">' + line_name + '<span class ="name_tip">' +  line_name + '</span></span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + ' 2">'
	for i in range(len(data)):
		div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + data[i] + '</span><span class="insert">&nbsp;</span></span>'
	div_seq += '</div>'

	return div_name, div_seq

def MakeDivAADonor(class_name, line_name, data, ori_seq, donor_region):
	div_name = '<div class="' + class_name + ' 1">'
	div_name += '<span class="name">' + line_name + '<span class ="name_tip">' +  line_name + '</span></span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + ' 2">'
	cur_pos = 1
	for i in range(len(data)):
		if ori_seq[i] == "-":
			div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + \
			           data[i] + '</span><span class="insert">&nbsp;</span></span>'
		else:
			if cur_pos in donor_region:
				div_seq += '<span class="unit_pack donor"><span class="insert">&nbsp;</span><span class="unit">' + \
			           data[i] + '</span><span class="insert">&nbsp;</span></span>'
			else:
				div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + \
				           data[i] + '</span><span class="insert">&nbsp;</span></span>'
			cur_pos += 1

	div_seq += '</div>'

	return div_name, div_seq

def MakeDivPosAA(class_name, line_name, tip_text, data):
	div_name = '<div class="' + class_name + '">'
	div_name += '<span class="name">' + line_name + '</span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + '">'
	for i in range(len(data[0])):
		if data[0][i] != '-':
			if int(data[0][i]) % 5 == 0:
				div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + str(data[0][i]) + \
				               '<span class ="unit_tip">' + tip_text + str(data[1][i]) + \
				               '</span></span><span class="insert">&nbsp;</span></span>'
			else:
				div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + '.' + \
				               '<span class ="unit_tip">' + tip_text + str(data[1][i]) + \
				               '</span></span><span class="insert">&nbsp;</span></span>'
		else:
			div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + str(data[0][i]) + \
			               '<span class ="unit_tip">' + tip_text + str(data[1][i]) + \
			               '</span></span><span class="insert">&nbsp;</span></span>'
	div_seq += '</div>'

	return div_name, div_seq

def MakeDivH1N3(class_name, line_name, tip_text, data):
	div_name = '<div class="' + class_name + '">'
	div_name += '<span class="name">' + line_name + '</span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + '">'
	for i in range(len(data)):
		if data[i][2] == '':
			if data[i][0] != '-':
				if int(data[i][0]) % 5.0 == 0:
					div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + str(data[i][0]) + \
					               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
					               '</span></span><span class="insert">&nbsp;</span></span>'
				else:
					div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + '.' + \
					               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
					               '</span></span><span class="insert">&nbsp;</span></span>'
			else:
				div_seq += '<span class="unit_pack"><span class="insert">&nbsp;</span><span class="unit">' + str(data[i][0]) + \
				               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
				               '</span></span><span class="insert">&nbsp;</span></span>'
		else:
			if data[i][0] != '-':
				if int(data[i][0]) % 5.0 == 0:
					div_seq += '<span class="unit_pack"><span class="insert ' + data[i][2] + '">&nbsp;</span><span class="unit ' + \
					               data[i][2] + '">' + str(data[i][0]) + \
					               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
					               '</span></span><span class="insert ' + data[i][2] + '">&nbsp;</span></span>'
				else:
					div_seq += '<span class="unit_pack"><span class="insert ' + data[i][2] + '">&nbsp;</span><span class="unit ' + \
					               data[i][2] + '">' + '.' + \
					               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
					               '</span></span><span class="insert ' + data[i][2] + '">&nbsp;</span></span>'
			else:
				div_seq += '<span class="unit_pack"><span class="insert ' + data[i][2] + '">&nbsp;</span><span class="unit ' + \
				               data[i][2] + '">' + str(data[i][0]) + \
				               '<span class ="unit_tip">' + tip_text + str(data[i][1]) + \
				               '</span></span><span class="insert ' + data[i][2] + '">&nbsp;</span></span>'
	div_seq += '</div>'

	return div_name, div_seq

def MakeDivPosNT(class_name, line_name, tip_text, data):
	div_name = '<div class="' + class_name + '">'
	div_name += '<span class="name">' + line_name + '</span>'
	div_name += '</div>'
	div_seq = '<div class="' + class_name + '">'
	count = 0
	for i in range(len(data[0])):
		if count == 0:
			div_seq += '<span class="unit_pack">'
			if data[0][i] % 5.0 == 0:
				div_seq += '<span class="unit">' + str(data[0][i]) + '<span class ="unit_tip">' + tip_text + \
				               str(data[1][i]) + ' - ' + str(int(data[1][i]) + 2) +  '</span></span>'
			else:
				div_seq += '<span class="unit">' + '.' + '<span class ="unit_tip">' + tip_text + str(data[1][i]) + \
				           ' - ' + str(int(data[1][i]) + 2) +  '</span></span>'
		elif count % 3 == 0:
			div_seq += '</span><span class="unit_pack">'
			if data[0][i] % 5.0 == 0:
				div_seq += '<span class="unit">' + str(data[0][i]) + '<span class ="unit_tip">' + tip_text + \
				               str(data[1][i])  + ' - ' + str(int(data[1][i]) + 2) +  '</span></span>'
			else:
				div_seq += '<span class="unit">' + '.' + '<span class ="unit_tip">' + tip_text + str(data[1][i]) + \
				           ' - ' + str(int(data[1][i]) + 2) +  '</span></span>'
		else:
			if data[0][i] % 5.0 == 0:
				div_seq += '<span class="unit">' + str(data[0][i]) + '</span>'
			else:
				div_seq += '<span class="unit">' + '.' + '</span>'
		count += 1
	div_seq += '</span>'
	div_seq += '</div>'
	return div_name, div_seq

def MakeSeqWithInseetion(class_name,id,AAseq,info):
	start_dict = {}
	end_dict = {}
	if len(info) > 0:
		for ele in info:
			start_dict[info[ele][0]] = ele
			end_dict[info[ele][1]] = ele

	div_seq = '<div class="' + class_name + '" id="' + id + '">'
	i = 0
	for aa in AAseq:
		pos = i + 1
		if end_dict.__contains__(pos):
			div_seq += '<span class="unit">' + aa + '</span>'
			div_seq += '<span class="insertion" style="margin-top: 10px;">' + info[end_dict[pos]][2] + '</span>'
			div_seq += '</span>'
			i += 1
			continue
		if start_dict.__contains__(pos):
			div_seq += '<span class="replace" id="' + str(start_dict[pos]) + '" title="' + info[start_dict[pos]][4] + '">'
			div_seq += '<span class="unit">' + aa + '</span>'
			i += 1
			continue
		div_seq += '<span class="unit">' + aa + '</span>'
		i += 1
	div_seq += '</div>'

	return div_seq

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

AACodonDict={'I':'ATT','L':'CTT','V':'GTT','F':'TTT','M':'ATG','C':'TGT',
			 'A':'GCT','G':'GGT','P':'CCT','T':'ACT','S':'TCT','Y':'TAT',
			 'W':'TGG','Q':'CAA','N':'AAT','H':'CAT','E':'GAA','D':'GAT',
			 'K':'AAA','R':'CGT'}

# Dictionaries indicating 1st nucleotide after the CDR3 ends in each J gene
JHuman = {'IGKJ1*01':8, 'IGKJ2*01':9, 'IGKJ2*02':8, 'IGKJ2*03':8,
		  'IGKJ2*04':9, 'IGKJ3*01':8, 'IGKJ4*01':8, 'IGKJ4*02':8,
		  'IGKJ5*01':8, 'IGHJ1*01':18, 'IGHJ2*01':19, 'IGHJ3*01':16,
		  'IGHJ3*02':16, 'IGHJ4*01':14, 'IGHJ4*02':14, 'IGHJ4*03':14,
		  'IGHJ5*01':17, 'IGHJ5*02':17, 'IGHJ6*01':29, 'IGHJ6*02':29,
		  'IGHJ6*03':29,    'IGHJ6*04':29, 'IGLJ1*01':8, 'IGLJ2*01':8,
		  'IGLJ3*01':8, 'IGLJ3*02':8, 'IGLJ4*01':8, 'IGLJ5*01':8,
		  'IGLJ5*02':8, 'IGLJ6*01':8, 'IGLJ7*01':8, 'IGLJ7*02':8}

JMouse = {'IGKJ1*01':8, 'IGKJ1*02':7, 'IGKJ2*01':9, 'IGKJ2*02':9,
		  'IGKJ2*03':9, 'IGKJ3*01':8, 'IGKJ3*02':8, 'IGKJ4*01':8,
		  'IGKJ4*02':8, 'IGKJ5*01':8, 'IGHJ1*01':19, 'IGHJ1*02':19,
		  'IGHJ1*03':19, 'IGHJ2*01':14, 'IGHJ2*02':14, 'IGHJ2*03':12,
		  'IGHJ3*01':14, 'IGHJ3*02':14, 'IGHJ4*01':20}

if __name__ == '__main__':
	import sys

	app = QtWidgets.QApplication(sys.argv)
	Vgenes = VGenesForm()

	Vgenes.ApplicationStarted()
	sys.exit(app.exec_())
