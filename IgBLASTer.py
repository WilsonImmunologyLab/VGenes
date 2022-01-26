__author__ = 'wilsonp'
# Command line for IgBLAST is:
# igblastn -germline_db_V mouse_gl_V -germline_db_J mouse_gl_J -germline_db_D mouse_gl_D -organism mouse -domain_system kabat -query MVtests1.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3
#
# where MVtests1 is the FASTA file name

# must first switch to a directory containing both the database and the FASTA file
# then use the subprocess command to run program and grab the output as a string
import os, time, sys, re
import sqlite3 as db
from PyQt5 import QtWidgets
from multiprocessing.dummy import Pool as ThreadPool
import VGenesDialogues #import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, setText
import VGenesSeq
from platform import system

# IgBlastThreadTest = ''

global working_prefix
global temp_folder
global igblast_path

working_prefix = os.path.dirname(os.path.realpath(sys.argv[0]))
temp_folder = os.path.join(working_prefix, 'Temp')

if system() == 'Windows':
	igblast_path = os.path.join(working_prefix, 'IgBlast', 'igblastn.exe')
else:
	igblast_path = os.path.join(working_prefix, 'IgBlast', 'igblastn')

def ProcessFASTAold(FASTAfile, MaxNum):
	ErLog = ''
	ErlogFile = ''

	TotSeqs = 0

	with open(FASTAfile, 'r') as currentFile:  #using with for this automatically closes the file even if you crash

		Newline = ''
		Titleline = ''
		ErLog = ''
		CleanSeq = ''
		for FASTAline in currentFile:

			FASTAline = FASTAline.replace('\n', '').replace('\r', '')
			if FASTAline == '':
				FASTAline = ' '

			if FASTAline[0] == '>':
				if len(Newline) > 1:
					Titleline += '\n'
					Newline += '\n'
					CleanSeq += Titleline + Newline
					TotSeqs += 1
				else:
					if Newline != ' ' or  Newline != '':
						if Titleline != '':
							ErLog += Titleline + ': Sequence error\n'

				Titleline = FASTAline
				Newline = ''
			else:
				FASTAline = FASTAline.upper()
				for nuc in FASTAline:
					if nuc == 'N' or nuc == 'A' or nuc == 'T' or nuc == 'G' or nuc == 'C'or nuc == '-':
						Newline += nuc
			if MaxNum > 0:
				if TotSeqs > MaxNum-1:
					break
		# need to write the last sequence into the FASTA file
		if len(Newline) > 1:
			Titleline += '\n'
			Newline += '\n'
			CleanSeq += Titleline + Newline
			TotSeqs += 1
		else:
			if Newline != ' ' or  Newline != '':
				ErLog += Titleline + '\n'

	# if ErLog != '':
	#     ErLog = 'The following sequences appear to have errors and were not processed: \n' + ErLog
	#     dirname = os.getcwd()
	#     now = 'Error log_' + time.strftime('%c')
	#     ErlogFile = os.path.join(dirname, now)
	#
	#     with open(ErlogFile, 'w') as currentFile:  #using with for this automatically closes the file even if you crash
	#         currentFile.write(ErLog)



	WorkingDir  = os.path.join(working_prefix, 'IgBlast')
	os.chdir(WorkingDir)

	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt')

	if CleanSeq == '':
		msg = 'There were no good variable gene seqeunces in this set'
		buttons = 'OK'
		# answer = VGenesDialogues.informationMessage(self, msg, buttons)
		return

	with open(workingfilename, 'w') as currentFile:
		currentFile.write(CleanSeq)

	return ErLog, TotSeqs


def ProcessFASTA(FASTAfile, MaxNum):
	ErLog = ''
	ErlogFile = ''

	TotSeqs = 0

	with open(FASTAfile, 'r') as currentFile:  #using with for this automatically closes the file even if you crash

		SeqName = ''
		Sequence = ''
		ErLog = ''
		CleanSeq = ''
		for FASTAline in currentFile:

			FASTAline = FASTAline.replace('\n', '').replace('\r', '')
			if FASTAline == '':
				FASTAline = ' '

			if FASTAline[0] == '>':
				if SeqName == "":
					pass
				else:
					SeqName = SeqNameParse(SeqName,35)
					CleanSeq += SeqName + '\n' + Sequence + '\n'
					TotSeqs += 1

				SeqName = FASTAline
				Sequence = ''
			else:
				#Sequence += re.sub(r'[^NATCG\-]','',FASTAline.upper())
				Sequence += FASTAline.upper()
			if MaxNum > 0:
				if TotSeqs > MaxNum-1:
					break
		# need to write the last sequence into the FASTA file
		SeqName = SeqNameParse(SeqName, 35)
		CleanSeq += SeqName + '\n' + Sequence + '\n'
		TotSeqs += 1

	WorkingDir = os.path.join(working_prefix, 'IgBlast')
	os.chdir(WorkingDir)

	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt')

	if CleanSeq == '':
		msg = 'There were no good variable gene seqeunces in this set'
		buttons = 'OK'
		# answer = VGenesDialogues.informationMessage(self, msg, buttons)
		return

	with open(workingfilename, 'w') as currentFile:
		currentFile.write(CleanSeq)

	return ErLog, TotSeqs

def SeqNameParse(SeqName, maxLen):
	SeqName = re.sub(r'[^\w\d\-\_\>]', '_', SeqName)
	if len(SeqName) > maxLen:
		SeqName = SeqName[0:maxLen]

	return SeqName

def TimeCheck(listname):
	countim = 0
	listim = 0
	print('Start ' + listname + ' counts at ' + time.strftime('%c'))
	for IgLine in IgBlastThreadTest:
		# testList += IgLine
		if IgLine[0:7] == 'Query= ':
			countim+= 1
			listim += 1
			if countim == 1000:
				print(listname + str(listim)+ ' at ' + time.strftime('%c'))
				countim = 0
	return countim

def IgBLASTit(FASTAFile, datalist, signal):
	import os
	#todo change to app folder
	progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
	ErlogFile = os.path.join(temp_folder, 'ErLog.txt')
	ErlogFile2 = os.path.join(temp_folder, 'ErLog2.txt')
	ErLog = 'VGenes input beginning at: '+ time.strftime('%c') + '\n'
	with open(ErlogFile, 'w') as currentFile:  #using with for this automatically closes the file even if you crash
		currentFile.write(ErLog)

	try:
		DBpathname = os.path.join(working_prefix, 'Data','VDJGenes.db')
		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)

		GetProductive = 2
		conn = db.connect(DBpathname)

	except:
		DBpathname = '/Volumes/Promise Pegasus/Dropbox/VGenes/VDJGenes.db'
		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)

		GetProductive = 2
		conn = db.connect(DBpathname)

	#  then need to create a cursor that lets you traverse the database
	cursor = conn.cursor()

	IgBLASTAnalysis = []
	IgBLASTset = []
	Sequences = {}

	project = datalist[0]
	grouping  = datalist[1]
	subgroup = datalist[2]
	species = datalist[3]
	GetProductive = datalist[4]
	MaxNum = int(datalist[5])

	try:
		multiProject = datalist[6]
	except:
		print('ops')

	# IgBlastOut = ''
	ErLog, TotSeqs = ProcessFASTA(FASTAFile, MaxNum)


	workingdir = os.path.join(working_prefix, 'IgBlast')
	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt')
	# add code in ProcessFASTA to also return number of seqs then split into 10,000 seq
	# sets to send through independently with different final names...might need split
	# IgBlastit out from this code to get parallel. Actually, it's after IgBLASTn that
	# needs to be split because BLASTngoes faster the more it has

	os.chdir(workingdir)

	Sequences.clear()
	IgBLASTset.clear()

	with open(workingfilename, 'r') as currentFile:  #make dictionary of seqs keyed by name
		# Sequences = fasta_dictionary(currentFile)

		# # Then use it to retrieve sequences with names as key
		# if len(currentFile) == 0:
		#     print('no')

		for FASTAline in currentFile:
			FASTAline = FASTAline.replace('\n', '').replace('\r', '')
			if FASTAline[0] == '>':
				#print(FASTAline)
				SeqNamed  = FASTAline[1:]
				SeqNamed = SeqNamed.strip()

			else:
				Sequence  = FASTAline
				#print(Sequence)
				if Sequence != '':
					Sequences[SeqNamed] = Sequence
				SeqNamed = ''
				Sequence = ''

	# TODO on deployment need to ensure base user /Applications/IgBlast/database contains igblastn, and databases

	try:
		start = time.time()
		if species == 'Human':
			BLASTCommandLine = igblast_path + " -germline_db_V IG/Human/HumanVGenes.nt -germline_db_J IG/Human/HumanJGenes.nt -germline_db_D IG/Human/HumanDGenes.nt -organism human -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt 3"
			IgBlastOut = os.popen(BLASTCommandLine)
		elif species == 'Mouse':
			BLASTCommandLine = igblast_path + " -germline_db_V IG/Mouse/MouseVGenes.nt -germline_db_J IG/Mouse/MouseJGenes.nt -germline_db_D IG/Mouse/MouseDGenes.nt -organism mouse -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3"
			IgBlastOut = os.popen(BLASTCommandLine)
		end = time.time()
		print('Run time for IgBlast: ' + str(end - start))
	except:
		ErLog = 'VGenes running Error!\nCurrent CMD: ' + BLASTCommandLine + '\n'
		with open(ErlogFile, 'a') as currentFile:  # using with for this automatically closes the file even if you crash
			currentFile.write(ErLog)

	# IgBlastOut += '\nBLASTend'
	# igblastn -germline_db_V human_gl_V_IMGT -germline_db_J human_gl_J_IMGT -germline_db_D human_gl_D_IMGT -organism human -domain_system kabat -query SFV-005H.nt -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt 3
	# IgBlastOut = os.popen("igblastn -germline_db_V mouse_gl_V -germline_db_J mouse_gl_J -germline_db_D mouse_gl_D -organism mouse -domain_system kabat -query MVtests1.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3")
	#print('stop')

	# todo Below is code to test IgBlast out put saving to a file...need comment out before use
	# testfilename = '/Applications/IgBlast/database/testfilename.txt'
	# with open(testfilename, 'w') as currentFile:
	#     for IgLine in IgBlastOut:
	#         currentFile.write(IgLine)
	# i = 0
	# SeqCount = 0
	# IgBlastList =  ''
	# for IgLine in IgBlastOut:
	#     IgBlastList += IgLine
	#     if IgLine[0:7] == 'Query= ':
	#         SeqCount += 1
	#         if SeqCount > 10000:  # if not the first in a new list
	#             print('test')
	#
	#

	testList = ''

	# ErlogFile2 = '/Applications/IgBlast/database/ErLog2.txt'  # NoErrors  NoGoodSeqs
	# ErLog2  = 'VGenes input beginning at: '+ time.strftime('%c') + '\n'


	# listread1 ='List1'
	# listread2 ='list2'
	# listread3 ='list3'
	# TimeCheck(IgBlastOut, listread1)
	# TimeCheck(IgBlastOut, listread2)
	# TimeCheck(IgBlastOut, listread3)
	#
	global IgBlastThreadTest
	# IgBlastThreadTest = IgBlastOut
	# listnames = ['list1', 'list2']#, 'list3','list4', 'list5', 'list6']
	# pool = ThreadPool(2)
	# results = pool.map(TimeCheck, listnames)
	# pool.close()
	# pool.join()
	#


	# ErLog2  = 'VGenes input ended at: '+ time.strftime('%c') + '\n'
	# return

	start = time.time()


	Importing  = False
	NotValid = False
	GrabAlignment = False
	GrabRecomb = False
	GrabJunction = False
	LineCount = 0
	TotMut = 0
	x = ''
	SeqAlignment = ''
	SeqNumber = 0
	BadSeqs = 0
	ErReport = ''
	GermlineSeq = ''
	GVgene = ''
	GDgene = ''
	GJgene = ''
	GVseq = []
	GDseq = []
	GJseq = []
	Aligned = []
	MutList = []
	DgeneLocus = ''

	Isotype = ''
	GVbeg = 0
	GVend = 0
	GD1beg = 0
	GD1end = 0
	GD2beg = 0
	GD2end = 0
	GJbeg = 0
	GJend = 0
	Vbeg = 0
	Vend = 0
	D1beg = 0
	D1end = 0
	D2beg = 0
	D2end = 0
	Jbeg = 0
	Jend = 0
	ORF = ""

	current_seq = 0
	totel_seq = len(Sequences)
	for IgLine in IgBlastOut:

		# Alignment is actually last thing that comes up but have to code first to retain /n and /r just for the alignment field
		if GrabAlignment == False:
			IgLine = IgLine.replace('\n', '').replace('\r', '')

		if IgLine[0:7] == 'Matrix:':
		# Analysis is completed but last seq analysis needs saved as no new seqname to trigger
			if NotValid == False:

				IgBLASTAnalysis.append(project)
				IgBLASTAnalysis.append(grouping)
				IgBLASTAnalysis.append(subgroup)
				IgBLASTAnalysis.append(species)
				try:

					Sequence = Sequences[IgBLASTAnalysis[0]]
				except:
					if SeqName:
						print(SeqName)
					else:
						return

				if Sequence != '':
					IgBLASTAnalysis.append(Sequence)
				else:
					NotValid = True
					ErLog = SeqName + 'has no sequence (line 210)\n'
					with open(ErlogFile, 'a') as currentfile:
						currentfile.write(ErLog)
					BadSeqs += 1
				#
				# TODO need code for following for now just added as blanks for database fields: can delete each as completed
				# CDR3DNA text, CDR3AA text, CDR3Length text, CDR3beg text, CDR3end text,


				IgBLASTAnalysis.append(GermlineSeq)
				IgBLASTAnalysis.append(CDR3DNA)
				IgBLASTAnalysis.append(CDR3AA)
				IgBLASTAnalysis.append(CDR3AALength)
				IgBLASTAnalysis.append(CDR3beg)
				IgBLASTAnalysis.append(CDR3end)
				IgBLASTAnalysis.append('Specificity')
				IgBLASTAnalysis.append('Subspecificity')
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(VGeneLocus)
				IgBLASTAnalysis.append(JgeneLocus)
				IgBLASTAnalysis.append(DgeneLocus)
				IgBLASTAnalysis.append(now)
				IgBLASTAnalysis.append(' ')
				IgBLASTAnalysis.append(' ')
				IgBLASTAnalysis.append(TotalMuts)
				IgBLASTAnalysis.append(MutListDB)
				IgBLASTAnalysis.append(IDevent)
				IgBLASTAnalysis.append(CDR3MW)
				IgBLASTAnalysis.append(CDR3pI)
				if Isotype == '':
					Isotype = "Unknown"
				IgBLASTAnalysis.append(Isotype)
				IgBLASTAnalysis.append(GCDR3beg)
				IgBLASTAnalysis.append(GCDR3end)
				IgBLASTAnalysis.append('Blank6')
				IgBLASTAnalysis.append(ORF)
				ORF = ''
				IgBLASTAnalysis.append('Blank8')
				IgBLASTAnalysis.append('Blank9')
				IgBLASTAnalysis.append('Blank10')
				IgBLASTAnalysis.append('Blank11')
				IgBLASTAnalysis.append('Blank12')
				IgBLASTAnalysis.append('Blank13')
				IgBLASTAnalysis.append('Blank14')
				IgBLASTAnalysis.append('Blank15')
				IgBLASTAnalysis.append('Blank16')
				IgBLASTAnalysis.append('Blank17')
				IgBLASTAnalysis.append('Blank18')
				IgBLASTAnalysis.append('Blank19')
				IgBLASTAnalysis.append('Blank20')

				# IgBLASTAnalysis.append('TheEnd')  #  entry to signal SQL generator properly

				# IgBLASTset[SeqNumber-1] = IgBLASTAnalysis  #  if not first record it saves

				if project == 'ByFunction':
					if Importing == False:
						IgBLASTAnalysis[75] = GeneType
						if multiProject != '':
							IgBLASTAnalysis[75] = multiProject
						if Productive == "Yes":
							IgBLASTAnalysis[76] = 'Functional'
						else:
							IgBLASTAnalysis[76] = 'Nonfunctional'

						if multiProject != '':
							IgBLASTAnalysis[76] = GeneType

					# grouping = Productive
					IgBLASTAnalysis[77] = ReadingFrame
					if multiProject != '':
						if Productive == "Yes":
							IgBLASTAnalysis[77] = 'Functional'
						else:
							IgBLASTAnalysis[77] = 'Nonfunctional'


				IgBLASTset.append((IgBLASTAnalysis))
			NotValid = False


			# if ErLog != '':
			NumberAnalyzed = SeqNumber - BadSeqs
			# ErLog = ErLog.replace('>', '')
			# if NumberAnalyzed > 0:
			#     ErLog = '\n' + str(NumberAnalyzed) + ' sequences were analyzed by IgBLAST \nThe following sequences appear to have errors and were not processed: \n' + ErLog
			# else:
			#     ErReport = 'NoGoodSeqs'
			if GetProductive == 1:  # only keep V and J
				# find record
				del_index = []
				for index in range(len(IgBLASTset)):
					if IgBLASTset[index][90] == 'N/A' or IgBLASTset[index][91] == 'N/A':
						del_index.append(index)

				for i in del_index:
					ErLog = IgBLASTset[i][0] + ' missed V or J\n'
					with open(ErlogFile, 'a') as currentfile:
						currentfile.write(ErLog)
				# delete record
				cnt = 0
				for index in del_index:
					del IgBLASTset[index - cnt]
					cnt += 1

			ErLog  = ErLog + '\nVGenes input ended at: '+ time.strftime('%c')
			with open(ErlogFile, 'a') as currentFile:  #using with for this automatically closes the file even if you crash
				currentFile.write(ErLog)

			# TODO need to use error log contents to determine if seqs should be added to database
			end = time.time()
			print('Run time for original mode: ' + str(end - start))

			return IgBLASTset

		if IgLine[0:10] == 'Alignments':
			GrabAlignment = True
		elif GrabAlignment == True:
			if NotValid == False:
				if IgLine[0:6] != "Lambda":
					SeqAlignment += IgLine
				else:
					GrabAlignment = False
					IgBLASTAnalysis.append(SeqAlignment)
					AlignParts = []

					# get ORF info from alignment
					if ORF == "":
						lines = SeqAlignment.split('\n')
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
								break
							line_num += 1

					AlignParts = ParseAlignment(SeqAlignment)
					if AlignParts == 'None':
						NotValid = True
						ErLog = SeqName + 'was problematic ParseAlignment of IgBLAST, line 323\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

					if AlignParts == 'short':
						NotValid = True
						# ErLog = SeqName + 'was problematic ParseAlignment of IgBLAST\n'
						# with open(ErlogFile, 'a') as currentfile:
						#     currentfile.write(ErLog)
						# BadSeqs += 1


					if NotValid == False:
						Dend = 0
						# for part in AlignParts:
						GVbeg = AlignParts[0]
						GVend = AlignParts[1]
						GD1beg = AlignParts[2]
						GD1end = AlignParts[3]
						GJbeg = AlignParts[4]
						GJend = AlignParts[5]
						AdjustEnd = AlignParts[6]
						SeqOrient = AlignParts[7]
						SeqBegin = AlignParts[8]
						SeqAdjust = SeqBegin - 1
						Vbeg = GVbeg + AdjustEnd
						Vend = GVend + AdjustEnd

						try:
							if SeqBegin > Vbeg and int(FR1From) == SeqBegin:
								IgBLASTAnalysis[22] = int(FR1From) - SeqAdjust
								IgBLASTAnalysis[23] = int(FR1To) - SeqAdjust
								IgBLASTAnalysis[29] = int(CDR1From) - SeqAdjust
								IgBLASTAnalysis[30] = int(CDR1To) - SeqAdjust
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust

								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust

									# print('no')
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(CDR1From) == SeqBegin:
								IgBLASTAnalysis[29] = int(CDR1From) - SeqAdjust
								IgBLASTAnalysis[30] = int(CDR1To) - SeqAdjust
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust
								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(FR2From) == SeqBegin:
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust
								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(CDR2From) == SeqBegin:
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(FR3From) == SeqBegin:
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust


							if SeqOrient == 'reversed':
								Sequence = Sequences[IgBLASTAnalysis[0]]
								Sequence2 = VGenesSeq.ReverseComp(Sequence)
								Sequences[IgBLASTAnalysis[0]]=Sequence2
								Sequence = Sequence2

							if SeqBegin > Vbeg:
								Sequence = Sequences[IgBLASTAnalysis[0]]
								Sequence = Sequence[SeqBegin-1:]
								Sequences[IgBLASTAnalysis[0]] = Sequence

							if GD1beg != 0:
								VDJun = IgBLASTAnalysis[17]
								if VDJun[:2] == 'N/A' or VDJun == '':  #  if no VD junction

									D1beg = Vend #+ 1
									D1end = Vend + len(IgBLASTAnalysis[18])

								elif VDJun[0] != '(':  # not overlapping with end of V
									D1beg = Vend + len(VDJun)+1
									D1end = D1beg + len(IgBLASTAnalysis[18])-1
								else:

									D1beg = Vend + 1
									D1end = Vend + len(IgBLASTAnalysis[18])-1
							else:
								D1beg = 0
								D1end = 0
						except:

							NotValid = True
							ErLog = SeqName + 'was problematic with sequence segments, line 348\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)

						GD2beg = 'GD2beg'
						GD2end = 'GD2end'
						D2beg = 'D2beg'
						D2end = 'D2end'


						if GeneType == 'Heavy':
							if IgBLASTAnalysis[17] != 'N/A':
								VD = IgBLASTAnalysis[17]
							else:
								VD = ''

							if IgBLASTAnalysis[19] != 'N/A':
								DJ = IgBLASTAnalysis[19]
							else:
								DJ = ''
							if DJ != '':

								try:
									if DJ[0] == '(':
										DJ = ''
								except:
									DJ = ''

							if VD != '':
								try:
									if VD[0] == '(':
										VD = ''
								except:
									VD = ''

							if GJbeg != 0:   # if there is a J
								Jbeg = (Vend+1) + len(VD) + (D1end-D1beg+1) +len(DJ)
								Jend = (Jbeg) + (GJend-GJbeg + 1)

							else:
								Jbeg = 0
								Jend = 0
						else:  # light chain
							try:
								Jjunc = IgBLASTAnalysis[21]
								if Jjunc == 'N/A':
									Jjunc = ''
									IgBLASTAnalysis[21] = ''

								Jbeg = Vend + len(Jjunc) + 1
								Jend = Jbeg + (GJend-GJbeg)
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1

						# if GJbeg != 0:
						#
						#     # Jend = int(IgBLASTAnalysis[1])
						#     Jend = GVend +
						#




						SeqDelin = [GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend]




						for segment in SeqDelin:
							IgBLASTAnalysis.append(segment)
						SeqAlignment = ''

						# Vgene =
						if species == 'Human':
							SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Vgene1 + '"'
							# SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'
							if GeneType == 'Heavy':
								SqlStatementD = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Dgene1 + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Jgene1 + '"'
						elif species == 'Mouse':
							SqlStatementV = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Vgene1 + '"'
							# SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'
							if GeneType == 'Heavy':
								SqlStatementD = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Dgene1 + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Jgene1 + '"'

						GVseq.clear()
						GDseq.clear()
						GJseq.clear()
						cursor.execute(SqlStatementV)
						for row in cursor:
							for column in row:
								GVseq.append(column)     # 0 SeqName, 1 Allele, 2 Species, 3 CodingStartNucleotide, 4 Sequence, 5 IMGTSequence
						if len(GVseq) < 6:
							VGeneLocus = 'NA'
							GVgene = 'NA'
							GStart = ''
							IMGTGVgene = 'NA'
							IMGTend = ''
							IMGTVDJ = 'NA'
						else:
							VGeneLocus = GVseq[0]
							GVgene = GVseq[4]
							GStart = int(GVseq[3])-1
							GVgene = GVgene[GStart:GVend]
							IMGTVgene = ''
							IMGTGVgene = GVseq[5]
							dotsIn = IMGTGVgene.count('.')
							# CharCount = 0

							# if GStart-1 != 0:
							#      #make sure no '.' in beginning (like if seq starts after CDR1)
							#     while CharCount != GStart-1:
							#         if IMGTVgene[CharCount] == '.':
							#             dotsIn -= 1

							IMGTend = GVend + dotsIn

							IMGTVDJ = IMGTGVgene[int(GVseq[3])-1:IMGTend]  #full length germline V with IMGT spacers



						if GeneType == 'Heavy':

							cursor.execute(SqlStatementD)
							for row in cursor:
								for column in row:
									GDseq.append(column)
							try:
								GDgene = GDseq[4]
								DgeneLocus = GDseq[0]
								GDgene = GDgene[(GD1beg-1):(GD1end)]
							except:
								GDgene = ''
								GDgeneLocus = 'N/A'


						else:
							GDgene = ''
							GDgeneLocus = 'N/A'
							DgeneLocus = 'N/A'
							for thing in range(0, 4):
								GDseq.append('N/A')

						cursor.execute(SqlStatementJ)
						for row in cursor:
							for column in row:
								GJseq.append(column)

							try:
								GJgene = GJseq[4]
								JgeneLocus = GJseq[0]
								GJgene = GJgene[(GJbeg-1):(GJend)]
							except:
								GJgene = ''
								GJgeneLocus = 'none'

						# GJgene = GJseq[4]
						# GJgene = GJgene[(GJbeg-1):(GJend)]


						if GeneType == 'Heavy':
							if IgBLASTAnalysis[17] != 'N/A':
								VD = IgBLASTAnalysis[17]
							else:
								VD = ''

							if IgBLASTAnalysis[19] != 'N/A':
								DJ = IgBLASTAnalysis[19]
							else:
								DJ = ''
							if DJ != '':
								try:
									if DJ[0] == '(':
										Overlap = len(DJ)-2
										GDgene = GDgene[Overlap:]
										DJ = ''
								except:
									DJ = ''
							if VD != '':
								try:
									if VD[0] == '(':
										Overlap = len(VD)-2
										GDgene = GDgene[Overlap:]
										VD = ''
								except:
									VD = ''


							GermlineSeq = GVgene + VD + GDgene + DJ + GJgene
							IMGTVDJ = IMGTVDJ + VD + GDgene + DJ + GJgene
							JunctionLength = len(VD) + len(GDgene) + len(DJ)

						else:

							try:
								GermlineSeq = GVgene + IgBLASTAnalysis[21] + GJgene
								IMGTVDJ = IMGTVDJ + IgBLASTAnalysis[21] + GJgene
								JunctionLength = len(IgBLASTAnalysis[21])
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1

						if int(GVend) < 200:
							NotValid = True
							ErLog = SeqName + 'V gene insufficiently long\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)
							BadSeqs += 1


						try:
							SeqIs = Sequences[IgBLASTAnalysis[0]]
						except:
							print("shit")
						# SeqIs = SeqIs[:Jend]
						# todo showing use of ClustalO
						ToClustalO = []
						Seq = (SeqName, SeqIs)
						ToClustalO.append(Seq)
						Seq = ('Germline', GermlineSeq)


						ToClustalO.append(Seq)
						if len(SeqIs) >= len(GermlineSeq):
							wraplength = len(SeqIs)+1
						else:
							wraplength = len(GermlineSeq)+1
						wraplength  = 80
						Aligned.clear()


						if NotValid == False:

							try:
								outfilename = VGenesSeq.ClustalO(ToClustalO,wraplength,False)
								if os.path.isfile(outfilename):
									# outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-out-seqs.fa')
									# while Aligned is None:
									Aligned = VGenesSeq.readClustalOutput(outfilename)
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1
								print(SeqName + ' goofed at clustal call in IgBlast line 777')


							finally:
								os.remove(outfilename)

							MutList.clear()

						try:
							MutList, SeqFixed, IDLength, AlignedSeq = VGenesSeq.Mutations(Aligned, GVend)
						except:
							print(SeqName)

						try:
							Sequences[IgBLASTAnalysis[0]] = SeqFixed
							if IDLength > 0:  #have to adjust for insertion...seq already adjusted for del
								IgBLASTAnalysis[69] = IgBLASTAnalysis[69]+ IDLength #  D1beg
								IgBLASTAnalysis[70] = IgBLASTAnalysis[70]+ IDLength #  D1end
								# IgBLASTAnalysis[71] = IgBLASTAnalysis[71]+ IDLength #  todo D2beg for when D2 found
								# IgBLASTAnalysis[72] = IgBLASTAnalysis[72]+ IDLength #  D2end
								IgBLASTAnalysis[73] = IgBLASTAnalysis[73]+ IDLength #  Jbeg
								IgBLASTAnalysis[74] = IgBLASTAnalysis[74]+ IDLength #  Jend
						except:
							NotValid = True
							ErLog = SeqName + 'was problematic (line 621)\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)
							BadSeqs += 1

						MutListDB = ''
						TotalMuts = 0
						IDevent = 'None'

						if SeqNumber >0:
							if NotValid == False:
								try:
									fw1end = int(IgBLASTAnalysis[23])
									cw1beg = int(IgBLASTAnalysis[29])
									cw1end = int(IgBLASTAnalysis[30])
									fw2beg = int(IgBLASTAnalysis[36])
									fw2end = int(IgBLASTAnalysis[37])
									cw2beg = int(IgBLASTAnalysis[43])
									cw2end = int(IgBLASTAnalysis[44])
									fw3beg = int(IgBLASTAnalysis[50])
									fw3end = int(IgBLASTAnalysis[51])
								except:
									NotValid = True
									ErLog = SeqName + 'was problematic (line 644)\n'
									with open(ErlogFile, 'a') as currentfile:
										currentfile.write(ErLog)
									BadSeqs += 1



					change = 0
					if NotValid == False:
						if len(MutList) > 0:
							for Mutation in MutList:
								listnum = 0
								for data in Mutation:
									MutListDB += str(data)
									if listnum <2: MutListDB+= '-'
									listnum +=1
									if data == 'Insertion' or data == 'Deletion':
										if IDevent == 'None':
											IDevent = data
										elif IDevent != data:
											IDevent = 'Both'
								#             (ins , position, sequence)
										if data == 'Insertion':
											change = len(Mutation[2])
										# elif data == 'Deletion':
										#     change = 0-len(Mutation[2])

										if SeqNumber >0:

											if int(Mutation[1]) < fw1end:
												IgBLASTAnalysis[23] = str(fw1end+change)
												IgBLASTAnalysis[29] = str(cw1beg+change)
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw1beg:
												IgBLASTAnalysis[29] = str(cw1beg+change)
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw1end:
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw2beg:
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw2end:
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw2beg:
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw2end:
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw3beg:
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw3end:
												IgBLASTAnalysis[51] = str(fw3end+change)


								MutListDB += ','
								try:
									TotalMuts += 1
								except:
									#print(SeqName)
									NotValid = True
									ErLog = SeqName + 'was problematic (line 735)\n'
									with open(ErlogFile, 'a') as currentfile:
										currentfile.write(ErLog)
									BadSeqs += 1

							try:
								if MutListDB[len(MutListDB)-1] == ',':
									MutListDB = MutListDB[:len(MutListDB)-1]
							except:
								#print(SeqName)
								NotValid = True
								ErLog = SeqName + 'was problematic (line 746)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1


					# FindCDR3(SeqFixed, IMGTVDJ, Vend, JGeneName, GJbeg, species, IDSeqLen, JunctionLength)
					# todo seq SFV-4C02 has CDR3 off because

					try:
						GCDR3beg, GCDR3end, CDR3DNA, CDR3AA, CDR3AALength, CDR3MW, CDR3pI = FindCDR3(AlignedSeq, IMGTVDJ, IgBLASTAnalysis[60], IgBLASTAnalysis[9], GJbeg, species, IDLength, JunctionLength) # DermVDJ, IMGTGermVDJ, Jgene1Name, Jbeg
						# IgBLASTAnalysis[51] = CDR3beg-1  # because IgBLAST made FW3 to end of V...
						CDR3beg = SeqFixed.find(CDR3DNA)+1
						CDR3end = CDR3beg + len(CDR3DNA)-1
						IgBLASTAnalysis[51] = CDR3beg-1

						now = time.strftime('%c')
					except:
						#print(SeqName)
						NotValid = True
						ErLog = SeqName + 'was problematic (line 766)\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

		if IgLine[0:7] == 'Query= ':
			# write current progress to file
			file_handle = open(progressBarFile,'w')
			progress = str(int(current_seq*100/totel_seq))
			file_handle.write(progress)
			file_handle.write(',' + str(current_seq) + '/' + str(totel_seq))
			file_handle.close()

			if signal == False:
				pass
			else:
				pct = int(current_seq*100/totel_seq)
				label = 'Processing: ' + str(current_seq) + '/' + str(totel_seq)
				signal.emit(pct, label)
			#print('Current Progress: ' + progress)
			current_seq += 1


			if SeqNumber >0:
				if NotValid == False:


					IgBLASTAnalysis.append(project)
					IgBLASTAnalysis.append(grouping)
					IgBLASTAnalysis.append(subgroup)
					IgBLASTAnalysis.append(species)

					Sequence = Sequences[IgBLASTAnalysis[0]]
					if Sequence != '':
						IgBLASTAnalysis.append(Sequence)
					else:
						NotValid = True
						ErLog = SeqName + 'has no sequence (line 210)\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

					if GeneType == 'Heavy':
						IsoSeq = (Sequence[(Jend):])
						#print(SeqName)
						IsoSeq = IsoSeq.strip('N')
						AGCTs = IsoSeq.count('A') + IsoSeq.count('G') +IsoSeq.count('C') +IsoSeq.count('T')
						if AGCTs > 5:  # todo decide if can determine isotype from < 5 or need more then
							#print('Start ' + SeqName + ' counts at ' + time.strftime('%c'))
							if species == 'Human':
								Isotype = VGenesSeq.CallIsotype(IsoSeq)
							elif species == 'Mouse':
								Isotype = VGenesSeq.CallIsotypeMouse(IsoSeq)
							else:
								Msg = 'Your current species is: ' + species + \
								      '\nWe do not support this species!'
								return Msg
							#print('end ' + SeqName + ' counts at ' + time.strftime('%c'))

						else:
							if len(IsoSeq) > 2:
								if IsoSeq[:3] == 'CCT' or IsoSeq == 'CTT':
									Isotype = 'IgG'
								elif IsoSeq[:3] == 'CAT':
									Isotype  = 'IgA'
								elif IsoSeq[:3] == 'GGA':
									Isotype  = 'IgM'
								elif IsoSeq[:3] == 'CAC':
									Isotype  = 'IgD'
								else:
									Isotype = IsoSeq
							else:
								Isotype = 'Unknown'
					else:
						if GeneType == 'Kappa':
							Isotype = 'Kappa'
						elif GeneType == 'Lambda':
							Isotype = 'Lambda'

					# TODO need code for following for now just added as blanks for database fields: can delete each as completed

					IgBLASTAnalysis.append(GermlineSeq)
					IgBLASTAnalysis.append(CDR3DNA)
					IgBLASTAnalysis.append(CDR3AA)
					IgBLASTAnalysis.append(CDR3AALength)
					IgBLASTAnalysis.append(CDR3beg)
					IgBLASTAnalysis.append(CDR3end)
					IgBLASTAnalysis.append('Specificity')
					IgBLASTAnalysis.append('Subspecificity')
					IgBLASTAnalysis.append(0)
					IgBLASTAnalysis.append(0)
					IgBLASTAnalysis.append(VGeneLocus)
					IgBLASTAnalysis.append(JgeneLocus)
					# if DgeneLocus == '':
					#     DgeneLocus = 'none'
					IgBLASTAnalysis.append(DgeneLocus)
					IgBLASTAnalysis.append(now)
					IgBLASTAnalysis.append(' ')
					IgBLASTAnalysis.append(' ')
					IgBLASTAnalysis.append(TotalMuts)
					IgBLASTAnalysis.append(MutListDB)
					IgBLASTAnalysis.append(IDevent)
					IgBLASTAnalysis.append(CDR3MW)
					IgBLASTAnalysis.append(CDR3pI)
					IgBLASTAnalysis.append(Isotype)
					IgBLASTAnalysis.append(GCDR3beg)
					IgBLASTAnalysis.append(GCDR3end)
					IgBLASTAnalysis.append('Blank6')
					IgBLASTAnalysis.append(ORF)
					ORF = ''
					IgBLASTAnalysis.append('Blank8')
					IgBLASTAnalysis.append('Blank9')
					IgBLASTAnalysis.append('Blank10')
					IgBLASTAnalysis.append('Blank11')
					IgBLASTAnalysis.append('Blank12')
					IgBLASTAnalysis.append('Blank13')
					IgBLASTAnalysis.append('Blank14')
					IgBLASTAnalysis.append('Blank15')
					IgBLASTAnalysis.append('Blank16')
					IgBLASTAnalysis.append('Blank17')
					IgBLASTAnalysis.append('Blank18')
					IgBLASTAnalysis.append('Blank19')
					IgBLASTAnalysis.append('Blank20')
					# IgBLASTAnalysis.append('TheEnd')

					# IgBLASTset[SeqNumber-1] = IgBLASTAnalysis  #  if not first record it saves
					if project == 'ByFunction':
						if Importing == False:
							IgBLASTAnalysis[75] = GeneType
							if multiProject != '':
								IgBLASTAnalysis[75] = multiProject
							if Productive == "Yes":
								IgBLASTAnalysis[76] = 'Functional'
							else:
								IgBLASTAnalysis[76]  = 'Nonfunctional'

							if multiProject != '':
								IgBLASTAnalysis[76] = GeneType

						# grouping = Productive
						IgBLASTAnalysis[77] = ReadingFrame
						if multiProject != '':
							if Productive == "Yes":
								IgBLASTAnalysis[77] = 'Functional'
							else:
								IgBLASTAnalysis[77]  = 'Nonfunctional'


					IgBLASTset.append((IgBLASTAnalysis))


				IgBLASTAnalysis = []  # can now empty IgBLASTANalysis for next seq
				NotValid = False




			SeqNumber += 1
			LineCount = 0
			TotMut = 0
			GrabAlignment = False
			GrabRecomb = False
			GrabJunction = False

			# NotValid = False
			SeqAlignment = ''
			FR1hit = False
			CDR1hit = False
			FW2hit = False
			CDR2hit = False




			SeqName = IgLine[7:]
			SeqName = SeqName.strip()  # removes any spaces before or after
			if SeqName[(len(SeqName)-6):] == 'Import':
				Importing = True
				SeqNameSplit = Segments.split("_")
				SeqName = SeqNameSplit[0]
				project = SeqNameSplit[1]
				grouping = SeqNameSplit[2]



			IgBLASTAnalysis.append(SeqName)
		elif IgLine[0:7] == 'Length=':
			LenSeq = int(IgLine[7:])

			IgBLASTAnalysis.append(str(LenSeq))
		elif IgLine[0:25] == '***** No hits found *****':
			NotValid = True

			ErLog = SeqName + ': No hits found\n'
			with open(ErlogFile, 'a') as currentfile:
				currentfile.write(ErLog)
			BadSeqs += 1
		elif IgLine[0:21] == 'V-(D)-J rearrangement':  # trick to get next line based on line count
			GrabRecomb = True  #  trick to get next line based on line count
		elif GrabRecomb == True: #  spring this only if previous line was V-(D)-J...
			# Should make this a different function to clean up
			GrabRecomb = False
			Segments = IgLine

			RecombParts = Segments.split("\t")
			#  Generates a list of each item in the string
			# seperated by tabs or any character '/n' is para mark

			GeneType = ""
			for segment in RecombParts:
				if segment == 'VH':
					GeneType = 'Heavy'
					break
				elif segment == 'VK':
					GeneType = 'Kappa'
					break
				else:
					GeneType = 'Lambda'

			IgBLASTAnalysis.append(GeneType)


			Vgene1 = ''
			Vgene2 = ''
			Vgene3 = ''
			Dgene1 = ''
			Dgene2 = ''
			Dgene3 = ''
			Jgene1 = ''
			Jgene2 = ''
			Jgene3 = ''

			Vgene = RecombParts[0]
			VChoice = Vgene.split(',')
			i = 0
			for i in range(0, len(VChoice)):
				Vgene1 = VChoice[0]

				if len(VChoice) == 2:
					Vgene2 = VChoice[1]
				elif len(VChoice) == 3:
					Vgene3 = VChoice[2]
			if Vgene1 != 'N/A':
				IgBLASTAnalysis.append(Vgene1)
				IgBLASTAnalysis.append(Vgene2)
				IgBLASTAnalysis.append(Vgene3)
			else:
				NotValid = True

				ErLog = SeqName + ': No V gene found\n'
				with open(ErlogFile, 'a') as currentfile:
					currentfile.write(ErLog)
				BadSeqs += 1

			IndexN = 0

			if GeneType == 'Heavy':
				IndexN = 1
			Dgene = RecombParts[IndexN]
			DChoice = Dgene.split(',')
			i = 0
			for i in range(0, len(DChoice)):
				Dgene1 = DChoice[0]
				if len(DChoice) == 2:
					Dgene2 = DChoice[1]
				elif len(DChoice) == 3:
					Dgene3 = DChoice[2]

			if GeneType == 'Heavy':  #only ned insert D for H chains
				IgBLASTAnalysis.append(Dgene1)
				IgBLASTAnalysis.append(Dgene2)
				IgBLASTAnalysis.append(Dgene3)
			else:
				IgBLASTAnalysis.append('None')
				IgBLASTAnalysis.append('None')
				IgBLASTAnalysis.append('None')


			IndexN += 1

			Jgene = RecombParts[IndexN]
			JChoice = Jgene.split(',')
			i = 0
			for i in range(0, len(JChoice)):
				Jgene1 = JChoice[0]
				if len(JChoice) == 2:
					Jgene2 = JChoice[1]
				elif len(JChoice) == 3:
					Jgene3 = JChoice[2]

			if Jgene1 != 'N/A':
				IgBLASTAnalysis.append(Jgene1)
				IgBLASTAnalysis.append(Jgene2)
				IgBLASTAnalysis.append(Jgene3)
			else:
				if GetProductive == 0:
					NotValid = True

				ErLog = SeqName + ': No J gene found\n'
				with open(ErlogFile, 'a') as currentfile:
					currentfile.write(ErLog)
				BadSeqs += 1

			IndexN += 2
			StopCodon = RecombParts[IndexN]
			IgBLASTAnalysis.append(StopCodon)

			IndexN += 1
			ReadingFrame = RecombParts[IndexN]
			IgBLASTAnalysis.append(ReadingFrame)

			IndexN += 1
			Productive = RecombParts[IndexN]
			#print(Productive)
			if GetProductive == 0:
				if Productive == 'Yes':
					IgBLASTAnalysis.append(Productive)
				else:
					NotValid = True
					ErLog = SeqName + ' was not a productive rearrangement\n'
					with open(ErlogFile, 'a') as currentfile:
						currentfile.write(ErLog)
					BadSeqs += 1
			else:
				IgBLASTAnalysis.append(Productive)

			IndexN += 1
			Strand = RecombParts[IndexN]
			IgBLASTAnalysis.append(Strand)
		elif IgLine[0:16] == 'V-(D)-J junction':  # trick to get next line based on line count
			GrabJunction = True  #  trick to get next line based on line count
			GrabIt = LineCount
		elif GrabJunction == True:
			GrabJunction = False
			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark
			VSeqend = RecombParts[0]
			IgBLASTAnalysis.append(VSeqend)

			if GeneType == "Heavy":
				VDJunction = RecombParts[1]
				Dregion = RecombParts[2]
				DJJunction = RecombParts[3]
				begJ = RecombParts[4]
				VJunction = ''  # only for Light chains but need to add field for database integrity

				IgBLASTAnalysis.append(VDJunction)
				IgBLASTAnalysis.append(Dregion)
				IgBLASTAnalysis.append(DJJunction)
				IgBLASTAnalysis.append(begJ)
				IgBLASTAnalysis.append(VJunction) # only for Light chains but need to add field for database integrity

			else: #  TODO need to set up for light chain so not out of whack with fewer fields

				VJunction = RecombParts[1]
				begJ = RecombParts[2]
				VDJunction = ''  # only for heavy chains but need to add field for database integrity
				Dregion = ''  # only for heavy chains but need to add field for database integrity
				DJJunction = ''  # only for heavy chains but need to add field for database integrity

				IgBLASTAnalysis.append(VDJunction)  # only for heavy chains but need to add field for database integrity
				IgBLASTAnalysis.append(Dregion)  # only for heavy chains but need to add field for database integrity
				IgBLASTAnalysis.append(DJJunction)  # only for heavy chains but need to add field for database integrity

				IgBLASTAnalysis.append(begJ)
				IgBLASTAnalysis.append(VJunction)
		elif IgLine[0:3] == 'FR1':
			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR1From = RecombParts[1]
			FR1To = RecombParts[2]
			FR1length = RecombParts[3]
			FR1matches = RecombParts[4]
			FR1mis = RecombParts[5]
			FR1gaps = RecombParts[6]
			FR1PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR1From)
			IgBLASTAnalysis.append(FR1To)
			IgBLASTAnalysis.append(FR1length)
			IgBLASTAnalysis.append(FR1matches)
			IgBLASTAnalysis.append(FR1mis)
			IgBLASTAnalysis.append(FR1gaps)
			IgBLASTAnalysis.append(FR1PercentIdentity)

			FR1hit = True



			TotMut += int(FR1mis)
			TotMut += int(FR1gaps)
		elif IgLine[0:4] == 'CDR1':
			if FR1hit == False:
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine

			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR1From = RecombParts[1]
			CDR1To = RecombParts[2]
			CDR1length = RecombParts[3]
			CDR1matches = RecombParts[4]
			CDR1mis = RecombParts[5]
			CDR1gaps = RecombParts[6]
			CDR1PercentIdentity = RecombParts[7]


			IgBLASTAnalysis.append(CDR1From)
			IgBLASTAnalysis.append(CDR1To)
			IgBLASTAnalysis.append(CDR1length)
			IgBLASTAnalysis.append(CDR1matches)
			IgBLASTAnalysis.append(CDR1mis)
			IgBLASTAnalysis.append(CDR1gaps)
			IgBLASTAnalysis.append(CDR1PercentIdentity)

			CDR1hit = True

			try:
				TotMut += int(CDR1mis)
			except:
				TotMut += 0

			try:
				TotMut += int(CDR1gaps)
			except:
				TotMut += 0
		elif IgLine[0:3] == 'FR2':
			if CDR1hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR2From = RecombParts[1]
			FR2To = RecombParts[2]
			FR2length = RecombParts[3]
			FR2matches = RecombParts[4]
			FR2mis = RecombParts[5]
			FR2gaps = RecombParts[6]
			FR2PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR2From)
			IgBLASTAnalysis.append(FR2To)
			IgBLASTAnalysis.append(FR2length)
			IgBLASTAnalysis.append(FR2matches)
			IgBLASTAnalysis.append(FR2mis)
			IgBLASTAnalysis.append(FR2gaps)
			IgBLASTAnalysis.append(FR2PercentIdentity)

			FW2hit = True

			TotMut += int(FR2mis)
			TotMut += int(FR2gaps)
		elif IgLine[0:4] == 'CDR2':
			if FW2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR2From = RecombParts[1]
			CDR2To = RecombParts[2]
			CDR2length = RecombParts[3]
			CDR2matches = RecombParts[4]
			CDR2mis = RecombParts[5]
			CDR2gaps = RecombParts[6]
			CDR2PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(CDR2From)
			IgBLASTAnalysis.append(CDR2To)
			IgBLASTAnalysis.append(CDR2length)
			IgBLASTAnalysis.append(CDR2matches)
			IgBLASTAnalysis.append(CDR2mis)
			IgBLASTAnalysis.append(CDR2gaps)
			IgBLASTAnalysis.append(CDR2PercentIdentity)

			CDR2hit = True

			if CDR2hit != 'N/A' and CDR2mis != 'N/A':
				try:
					TotMut += int(CDR2mis)
					TotMut += int(CDR2gaps)
				except:
					print('oops')
		elif IgLine[0:3] == 'FR3':
			if CDR2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #CDR2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)


			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR3From = RecombParts[1]
			FR3To = RecombParts[2]
			FR3length = RecombParts[3]
			FR3matches = RecombParts[4]
			FR3mis = RecombParts[5]
			FR3gaps = RecombParts[6]
			FR3PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR3From)
			IgBLASTAnalysis.append(FR3To)
			IgBLASTAnalysis.append(FR3length)
			IgBLASTAnalysis.append(FR3matches)
			IgBLASTAnalysis.append(FR3mis)
			IgBLASTAnalysis.append(FR3gaps)
			IgBLASTAnalysis.append(FR3PercentIdentity)
			try:
				if FR3mis == 'N/A':
					FR3mis = 0
				if FR3gaps == 'N/A':
					FR3gaps = 0
				TotMut += int(FR3mis)
				TotMut += int(FR3gaps)
			except:
				print('mistake')

			IgBLASTAnalysis.append(str(TotMut))
		# elif IgLine[0:4] == 'CDR3':

		# SeqNumber += 1

		# TODO maybe need to analyze IgBLAST output for SHM

	conn.close()

	end = time.time()
	print('Run time for original mode: ' + str(end - start))

	ErLog = '\nVGenes input ended at: ' + time.strftime('%c')
	with open(ErlogFile2, 'a') as currentFile:  # using with for this automatically closes the file even if you crash
		currentFile.write(ErLog)

def IgBLASTitResults(FASTAFile, IgBlastOutFile, datalist, signal):
	import os
	#todo change to app folder
	#progressBarFile = os.path.join(temp_folder, 'progressBarFile.txt')
	ErlogFile = os.path.join(working_prefix, 'IgBlast', 'ErLog.txt')  # '/Applications/IgBlast/database/ErLog.txt'  # NoErrors  NoGoodSeqs
	ErLog  = 'VGenes input beginning at: '+ time.strftime('%c') + '\n'
	with open(ErlogFile, 'w') as currentFile:  #using with for this automatically closes the file even if you crash
		currentFile.write(ErLog)

	try:
		DBpathname = os.path.join(working_prefix, 'Data','VDJGenes.db')
		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)

		GetProductive = 2
		conn = db.connect(DBpathname)

	except:
		DBpathname = '/Volumes/Promise Pegasus/Dropbox/VGenes/VDJGenes.db'
		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)

		GetProductive = 2
		conn = db.connect(DBpathname)

	#  then need to create a cursor that lets you traverse the database
	cursor = conn.cursor()

	IgBLASTAnalysis = []
	IgBLASTset = []
	Sequences = {}

	project = datalist[0]
	grouping  = datalist[1]
	subgroup = datalist[2]
	species = datalist[3]
	GetProductive = datalist[4]
	MaxNum = int(datalist[5])

	try:
		multiProject = datalist[6]
	except:
		print('ops')

	# IgBlastOut = ''
	ErLog, TotSeqs = ProcessFASTA(FASTAFile, MaxNum)

	workingdir = os.path.join(working_prefix, 'IgBlast')
	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt')
	# add code in ProcessFASTA to also return number of seqs then split into 10,000 seq
	# sets to send through independently with different final names...might need split
	# IgBlastit out from this code to get parallel. Actually, it's after IgBLASTn that
	# needs to be split because BLASTngoes faster the more it has

	os.chdir(workingdir)

	Sequences.clear()
	IgBLASTset.clear()

	with open(workingfilename, 'r') as currentFile:  # make dictionary of seqs keyed by name
		# Sequences = fasta_dictionary(currentFile)

		# # Then use it to retrieve sequences with names as key
		# if len(currentFile) == 0:
		#     print('no')

		for FASTAline in currentFile:
			FASTAline = FASTAline.replace('\n', '').replace('\r', '')
			if FASTAline[0] == '>':
				# print(FASTAline)
				SeqNamed = FASTAline[1:]
				SeqNamed = SeqNamed.strip()

			else:
				Sequence = FASTAline
				# print(Sequence)
				if Sequence != '':
					Sequences[SeqNamed] = Sequence
				SeqNamed = ''
				Sequence = ''

	IgBlastOut = ''
	fh = open(IgBlastOutFile,'r')
	IgBlastOut = fh.readlines()

	#
	global IgBlastThreadTest
	# IgBlastThreadTest = IgBlastOut
	# listnames = ['list1', 'list2']#, 'list3','list4', 'list5', 'list6']
	# pool = ThreadPool(2)
	# results = pool.map(TimeCheck, listnames)
	# pool.close()
	# pool.join()
	#


	# ErLog2  = 'VGenes input ended at: '+ time.strftime('%c') + '\n'
	# return



	Importing  = False
	NotValid = False
	GrabAlignment = False
	GrabRecomb = False
	GrabJunction = False
	LineCount = 0
	TotMut = 0
	x = ''
	SeqAlignment = ''
	SeqNumber = 0
	BadSeqs = 0
	ErReport = ''
	GermlineSeq = ''
	GVgene = ''
	GDgene = ''
	GJgene = ''
	GVseq = []
	GDseq = []
	GJseq = []
	Aligned = []
	MutList = []
	DgeneLocus = ''

	Isotype = ''
	GVbeg = 0
	GVend = 0
	GD1beg = 0
	GD1end = 0
	GD2beg = 0
	GD2end = 0
	GJbeg = 0
	GJend = 0
	Vbeg = 0
	Vend = 0
	D1beg = 0
	D1end = 0
	D2beg = 0
	D2end = 0
	Jbeg = 0
	Jend = 0

	current_seq = 0
	totel_seq = len(Sequences)
	for IgLine in IgBlastOut:

		# Alignment is actually last thing that comes up but have to code first to retain /n and /r just for the alignment field
		if GrabAlignment == False:
			IgLine = IgLine.replace('\n', '').replace('\r', '')

		if IgLine[0:7] == 'Matrix:':
		# Analysis is completed but last seq analysis needs saved as no new seqname to trigger
			if NotValid == False:

				IgBLASTAnalysis.append(project)
				IgBLASTAnalysis.append(grouping)
				IgBLASTAnalysis.append(subgroup)
				IgBLASTAnalysis.append(species)
				try:

					Sequence = Sequences[IgBLASTAnalysis[0]]
				except:
					if SeqName:
						print(SeqName)
					else:
						return

				if Sequence != '':
					IgBLASTAnalysis.append(Sequence)
				else:
					NotValid = True
					ErLog = SeqName + 'has no sequence (line 210)\n'
					with open(ErlogFile, 'a') as currentfile:
						currentfile.write(ErLog)
					BadSeqs += 1
				#
				# TODO need code for following for now just added as blanks for database fields: can delete each as completed
				# CDR3DNA text, CDR3AA text, CDR3Length text, CDR3beg text, CDR3end text,


				IgBLASTAnalysis.append(GermlineSeq)
				IgBLASTAnalysis.append(CDR3DNA)
				IgBLASTAnalysis.append(CDR3AA)
				IgBLASTAnalysis.append(CDR3AALength)
				IgBLASTAnalysis.append(CDR3beg)
				IgBLASTAnalysis.append(CDR3end)
				IgBLASTAnalysis.append('Specificity')
				IgBLASTAnalysis.append('Subspecificity')
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(VGeneLocus)
				IgBLASTAnalysis.append(JgeneLocus)
				IgBLASTAnalysis.append(DgeneLocus)
				IgBLASTAnalysis.append(now)
				IgBLASTAnalysis.append(' ')
				IgBLASTAnalysis.append(' ')
				IgBLASTAnalysis.append(TotalMuts)
				IgBLASTAnalysis.append(MutListDB)
				IgBLASTAnalysis.append(IDevent)
				IgBLASTAnalysis.append(CDR3MW)
				IgBLASTAnalysis.append(CDR3pI)
				if Isotype == '':
					Isotype = "Unknown"
				IgBLASTAnalysis.append(Isotype)
				IgBLASTAnalysis.append(GCDR3beg)
				IgBLASTAnalysis.append(GCDR3end)
				IgBLASTAnalysis.append('Blank6')
				IgBLASTAnalysis.append('Blank7')
				IgBLASTAnalysis.append('Blank8')
				IgBLASTAnalysis.append('Blank9')
				IgBLASTAnalysis.append('Blank10')
				IgBLASTAnalysis.append('Blank11')
				IgBLASTAnalysis.append('Blank12')
				IgBLASTAnalysis.append('Blank13')
				IgBLASTAnalysis.append('Blank14')
				IgBLASTAnalysis.append('Blank15')
				IgBLASTAnalysis.append('Blank16')
				IgBLASTAnalysis.append('Blank17')
				IgBLASTAnalysis.append('Blank18')
				IgBLASTAnalysis.append('Blank19')
				IgBLASTAnalysis.append('Blank20')

				# IgBLASTAnalysis.append('TheEnd')  #  entry to signal SQL generator properly

				# IgBLASTset[SeqNumber-1] = IgBLASTAnalysis  #  if not first record it saves

				if project == 'ByFunction':
					if Importing == False:
						IgBLASTAnalysis[75] = GeneType
						if multiProject != '':
							IgBLASTAnalysis[75] = multiProject
						if Productive == "Yes":
							IgBLASTAnalysis[76] = 'Functional'
						else:
							IgBLASTAnalysis[76] = 'Nonfunctional'

						if multiProject != '':
							IgBLASTAnalysis[76] = GeneType

					# grouping = Productive
					IgBLASTAnalysis[77] = ReadingFrame
					if multiProject != '':
						if Productive == "Yes":
							IgBLASTAnalysis[77] = 'Functional'
						else:
							IgBLASTAnalysis[77] = 'Nonfunctional'


				IgBLASTset.append((IgBLASTAnalysis))
			NotValid = False


			# if ErLog != '':
			NumberAnalyzed = SeqNumber - BadSeqs
			# ErLog = ErLog.replace('>', '')
			# if NumberAnalyzed > 0:
			#     ErLog = '\n' + str(NumberAnalyzed) + ' sequences were analyzed by IgBLAST \nThe following sequences appear to have errors and were not processed: \n' + ErLog
			# else:
			#     ErReport = 'NoGoodSeqs'

				# TODO make ErLog one permanent file deleted every time unless user wants to save on Main query

			ErLog  = ErLog + '\nVGenes input ended at: '+ time.strftime('%c')
			with open(ErlogFile, 'a') as currentFile:  #using with for this automatically closes the file even if you crash
				currentFile.write(ErLog)

			# TODO need to use error log contents to determine if seqs should be added to database

			return IgBLASTset

		if IgLine[0:10] == 'Alignments':
			GrabAlignment = True
		elif GrabAlignment == True:
			if NotValid == False:
				if IgLine[0:6] != "Lambda":
					SeqAlignment += IgLine
				else:
					GrabAlignment = False
					IgBLASTAnalysis.append(SeqAlignment)
					AlignParts = []
					AlignParts = ParseAlignment(SeqAlignment)
					if AlignParts == 'None':
						NotValid = True
						ErLog = SeqName + 'was problematic ParseAlignment of IgBLAST, line 323\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

					if AlignParts == 'short':
						NotValid = True
						# ErLog = SeqName + 'was problematic ParseAlignment of IgBLAST\n'
						# with open(ErlogFile, 'a') as currentfile:
						#     currentfile.write(ErLog)
						# BadSeqs += 1


					if NotValid == False:
						Dend = 0
						# for part in AlignParts:
						GVbeg = AlignParts[0]
						GVend = AlignParts[1]
						GD1beg = AlignParts[2]
						GD1end = AlignParts[3]
						GJbeg = AlignParts[4]
						GJend = AlignParts[5]
						AdjustEnd = AlignParts[6]
						SeqOrient = AlignParts[7]
						SeqBegin = AlignParts[8]
						SeqAdjust = SeqBegin - 1
						#Vbeg = GVbeg + AdjustEnd
						#Vend = GVend + AdjustEnd
						Vbeg = GVbeg
						Vend = GVend

						try:
							if SeqBegin > Vbeg and int(FR1From) == SeqBegin:
								IgBLASTAnalysis[22] = int(FR1From) - SeqAdjust
								IgBLASTAnalysis[23] = int(FR1To) - SeqAdjust
								IgBLASTAnalysis[29] = int(CDR1From) - SeqAdjust
								IgBLASTAnalysis[30] = int(CDR1To) - SeqAdjust
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust

								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust

									# print('no')
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(CDR1From) == SeqBegin:
								IgBLASTAnalysis[29] = int(CDR1From) - SeqAdjust
								IgBLASTAnalysis[30] = int(CDR1To) - SeqAdjust
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust
								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(FR2From) == SeqBegin:
								IgBLASTAnalysis[36] = int(FR2From) - SeqAdjust
								IgBLASTAnalysis[37] = int(FR2To) - SeqAdjust
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(CDR2From) == SeqBegin:
								IgBLASTAnalysis[43] = int(CDR2From) - SeqAdjust
								IgBLASTAnalysis[44] = int(CDR2To) - SeqAdjust
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust

							elif SeqBegin > Vbeg and int(FR3From) == SeqBegin:
								IgBLASTAnalysis[50] = int(FR3From) - SeqAdjust
								IgBLASTAnalysis[51] = int(FR3To) - SeqAdjust


							if SeqOrient == 'reversed':
								Sequence = Sequences[IgBLASTAnalysis[0]]
								Sequence2 = VGenesSeq.ReverseComp(Sequence)
								Sequences[IgBLASTAnalysis[0]]=Sequence2
								Sequence = Sequence2

							if SeqBegin > Vbeg:
								Sequence = Sequences[IgBLASTAnalysis[0]]
								Sequence = Sequence[SeqBegin-1:]
								Sequences[IgBLASTAnalysis[0]] = Sequence

							if GD1beg != 0:
								VDJun = IgBLASTAnalysis[17]
								if VDJun[:2] == 'N/A' or VDJun == '':  #  if no VD junction

									D1beg = Vend #+ 1
									D1end = Vend + len(IgBLASTAnalysis[18])

								elif VDJun[0] != '(':  # not overlapping with end of V
									D1beg = Vend + len(VDJun)+1
									D1end = D1beg + len(IgBLASTAnalysis[18])-1
								else:

									D1beg = Vend + 1
									D1end = Vend + len(IgBLASTAnalysis[18])-1
							else:
								D1beg = 0
								D1end = 0
						except:

							NotValid = True
							ErLog = SeqName + 'was problematic with sequence segments, line 348\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)

						GD2beg = 'GD2beg'
						GD2end = 'GD2end'
						D2beg = 'D2beg'
						D2end = 'D2end'


						if GeneType == 'Heavy':
							if IgBLASTAnalysis[17] != 'N/A':
								VD = IgBLASTAnalysis[17]
							else:
								VD = ''

							if IgBLASTAnalysis[19] != 'N/A':
								DJ = IgBLASTAnalysis[19]
							else:
								DJ = ''
							if DJ != '':

								try:
									if DJ[0] == '(':
										DJ = ''
								except:
									DJ = ''

							if VD != '':
								try:
									if VD[0] == '(':
										VD = ''
								except:
									VD = ''

							if GJbeg != 0:   # if there is a J
								Jbeg = (Vend+1) + len(VD) + (D1end-D1beg+1) +len(DJ)
								Jend = (Jbeg) + (GJend-GJbeg + 1)

							else:
								Jbeg = 0
								Jend = 0
						else:  # light chain
							try:
								Jjunc = IgBLASTAnalysis[21]
								if Jjunc == 'N/A':
									Jjunc = ''
									IgBLASTAnalysis[21] = ''

								Jbeg = Vend + len(Jjunc) + 1
								Jend = Jbeg + (GJend-GJbeg)
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1

						# if GJbeg != 0:
						#
						#     # Jend = int(IgBLASTAnalysis[1])
						#     Jend = GVend +
						#




						SeqDelin = [GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend]




						for segment in SeqDelin:
							IgBLASTAnalysis.append(segment)
						SeqAlignment = ''

						# Vgene =
						if species == 'Human':
							SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Vgene1 + '"'
							# SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'
							if GeneType == 'Heavy':
								SqlStatementD = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Dgene1 + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Human" AND Allele = "' + Jgene1 + '"'
						elif species == 'Mouse':
							SqlStatementV = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Vgene1 + '"'
							# SqlStatementV = 'SELECT SeqName, Allele, Species, CodingStartNucleotide, Sequence FROM GermLineDB'
							if GeneType == 'Heavy':
								SqlStatementD = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Dgene1 + '"'
							SqlStatementJ = 'SELECT SeqName, Allele, Strain, CodingStartNucleotide, Sequence, IMGTSequence FROM GermLineDB WHERE Species = "Mouse" AND Allele = "' + Jgene1 + '"'

						GVseq.clear()
						GDseq.clear()
						GJseq.clear()
						cursor.execute(SqlStatementV)
						for row in cursor:
							for column in row:
								GVseq.append(column)     # 0 SeqName, 1 Allele, 2 Species, 3 CodingStartNucleotide, 4 Sequence, 5 IMGTSequence
						VGeneLocus = GVseq[0]
						GVgene = GVseq[4]
						GStart = int(GVseq[3])-1
						GVgene = GVgene[GStart:GVend]
						IMGTVgene = ''
						IMGTGVgene = GVseq[5]
						dotsIn = IMGTGVgene.count('.')
						# CharCount = 0

						# if GStart-1 != 0:
						#      #make sure no '.' in beginning (like if seq starts after CDR1)
						#     while CharCount != GStart-1:
						#         if IMGTVgene[CharCount] == '.':
						#             dotsIn -= 1

						IMGTend = GVend + dotsIn

						IMGTVDJ = IMGTGVgene[int(GVseq[3])-1:IMGTend]  #full length germline V with IMGT spacers



						if GeneType == 'Heavy':

							cursor.execute(SqlStatementD)
							for row in cursor:
								for column in row:
									GDseq.append(column)
							try:
								GDgene = GDseq[4]
								DgeneLocus = GDseq[0]
								GDgene = GDgene[(GD1beg-1):(GD1end)]
							except:
								GDgene = ''
								GDgeneLocus = 'N/A'


						else:
							GDgene = ''
							GDgeneLocus = 'N/A'
							DgeneLocus = 'N/A'
							for thing in range(0, 4):
								GDseq.append('N/A')

						cursor.execute(SqlStatementJ)
						for row in cursor:
							for column in row:
								GJseq.append(column)

							try:
								GJgene = GJseq[4]
								JgeneLocus = GJseq[0]
								GJgene = GJgene[(GJbeg-1):(GJend)]
							except:
								GJgene = ''
								GJgeneLocus = 'none'

						# GJgene = GJseq[4]
						# GJgene = GJgene[(GJbeg-1):(GJend)]


						if GeneType == 'Heavy':
							if IgBLASTAnalysis[17] != 'N/A':
								VD = IgBLASTAnalysis[17]
							else:
								VD = ''

							if IgBLASTAnalysis[19] != 'N/A':
								DJ = IgBLASTAnalysis[19]
							else:
								DJ = ''
							if DJ != '':
								try:
									if DJ[0] == '(':
										Overlap = len(DJ)-2
										GDgene = GDgene[Overlap:]
										DJ = ''
								except:
									DJ = ''
							if VD != '':
								try:
									if VD[0] == '(':
										Overlap = len(VD)-2
										GDgene = GDgene[Overlap:]
										VD = ''
								except:
									VD = ''


							GermlineSeq = GVgene + VD + GDgene + DJ + GJgene
							IMGTVDJ = IMGTVDJ + VD + GDgene + DJ + GJgene
							JunctionLength = len(VD) + len(GDgene) + len(DJ)

						else:

							try:
								GermlineSeq = GVgene + IgBLASTAnalysis[21] + GJgene
								IMGTVDJ = IMGTVDJ + IgBLASTAnalysis[21] + GJgene
								JunctionLength = len(IgBLASTAnalysis[21])
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1

						if int(GVend) < 200:
							NotValid = True
							ErLog = SeqName + 'V gene insufficiently long\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)
							BadSeqs += 1


						try:
							SeqIs = Sequences[IgBLASTAnalysis[0]]
						except:
							print("shit")
						# SeqIs = SeqIs[:Jend]
						# todo showing use of ClustalO
						ToClustalO = []
						Seq = (SeqName, SeqIs)
						ToClustalO.append(Seq)
						Seq = ('Germline', GermlineSeq)


						ToClustalO.append(Seq)
						if len(SeqIs) >= len(GermlineSeq):
							wraplength = len(SeqIs)+1
						else:
							wraplength = len(GermlineSeq)+1
						wraplength  = 80
						Aligned.clear()


						if NotValid == False:

							try:
								outfilename = VGenesSeq.ClustalO(ToClustalO,wraplength,False)
								if os.path.isfile(outfilename):
									# outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'ClustalOmega', 'my-out-seqs.fa')
									# while Aligned is None:
									Aligned = VGenesSeq.readClustalOutput(outfilename)
							except:
								NotValid = True
								ErLog = SeqName + 'was problematic (line 621)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1
								print(SeqName + ' goofed at clustal call in IgBlast line 777')


							finally:
								os.remove(outfilename)

							MutList.clear()

						try:
							MutList, SeqFixed, IDLength, AlignedSeq = VGenesSeq.Mutations(Aligned, GVend)
						except:
							print(SeqName)

						try:
							Sequences[IgBLASTAnalysis[0]] = SeqFixed
							if IDLength > 0:  #have to adjust for insertion...seq already adjusted for del
								IgBLASTAnalysis[69] = IgBLASTAnalysis[69]+ IDLength #  D1beg
								IgBLASTAnalysis[70] = IgBLASTAnalysis[70]+ IDLength #  D1end
								# IgBLASTAnalysis[71] = IgBLASTAnalysis[71]+ IDLength #  todo D2beg for when D2 found
								# IgBLASTAnalysis[72] = IgBLASTAnalysis[72]+ IDLength #  D2end
								IgBLASTAnalysis[73] = IgBLASTAnalysis[73]+ IDLength #  Jbeg
								IgBLASTAnalysis[74] = IgBLASTAnalysis[74]+ IDLength #  Jend
						except:
							NotValid = True
							ErLog = SeqName + 'was problematic (line 621)\n'
							with open(ErlogFile, 'a') as currentfile:
								currentfile.write(ErLog)
							BadSeqs += 1

						MutListDB = ''
						TotalMuts = 0
						IDevent = 'None'

						if SeqNumber >0:
							if NotValid == False:
								try:
									fw1end = int(IgBLASTAnalysis[23])
									cw1beg = int(IgBLASTAnalysis[29])
									cw1end = int(IgBLASTAnalysis[30])
									fw2beg = int(IgBLASTAnalysis[36])
									fw2end = int(IgBLASTAnalysis[37])
									cw2beg = int(IgBLASTAnalysis[43])
									cw2end = int(IgBLASTAnalysis[44])
									fw3beg = int(IgBLASTAnalysis[50])
									fw3end = int(IgBLASTAnalysis[51])
								except:
									NotValid = True
									ErLog = SeqName + 'was problematic (line 644)\n'
									with open(ErlogFile, 'a') as currentfile:
										currentfile.write(ErLog)
									BadSeqs += 1



					change = 0
					if NotValid == False:
						if len(MutList) > 0:
							for Mutation in MutList:
								listnum = 0
								for data in Mutation:
									MutListDB += str(data)
									if listnum <2: MutListDB+= '-'
									listnum +=1
									if data == 'Insertion' or data == 'Deletion':
										if IDevent == 'None':
											IDevent = data
										elif IDevent != data:
											IDevent = 'Both'
								#             (ins , position, sequence)
										if data == 'Insertion':
											change = len(Mutation[2])
										# elif data == 'Deletion':
										#     change = 0-len(Mutation[2])

										if SeqNumber >0:

											if int(Mutation[1]) < fw1end:
												IgBLASTAnalysis[23] = str(fw1end+change)
												IgBLASTAnalysis[29] = str(cw1beg+change)
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw1beg:
												IgBLASTAnalysis[29] = str(cw1beg+change)
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw1end:
												IgBLASTAnalysis[30] = str(cw1end+change)
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw2beg:
												IgBLASTAnalysis[36] = str(fw2beg+change)
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw2end:
												IgBLASTAnalysis[37] = str(fw2end+change)
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw2beg:
												IgBLASTAnalysis[43] = str(cw2beg+change)
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < cw2end:
												IgBLASTAnalysis[44] = str(cw2end+change)
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw3beg:
												IgBLASTAnalysis[50] = str(fw3beg+change)
												IgBLASTAnalysis[51] = str(fw3end+change)
											elif int(Mutation[1]) < fw3end:
												IgBLASTAnalysis[51] = str(fw3end+change)


								MutListDB += ','
								try:
									TotalMuts += 1
								except:
									#print(SeqName)
									NotValid = True
									ErLog = SeqName + 'was problematic (line 735)\n'
									with open(ErlogFile, 'a') as currentfile:
										currentfile.write(ErLog)
									BadSeqs += 1

							try:
								if MutListDB[len(MutListDB)-1] == ',':
									MutListDB = MutListDB[:len(MutListDB)-1]
							except:
								#print(SeqName)
								NotValid = True
								ErLog = SeqName + 'was problematic (line 746)\n'
								with open(ErlogFile, 'a') as currentfile:
									currentfile.write(ErLog)
								BadSeqs += 1


					# FindCDR3(SeqFixed, IMGTVDJ, Vend, JGeneName, GJbeg, species, IDSeqLen, JunctionLength)
					# todo seq SFV-4C02 has CDR3 off because

					try:
						GCDR3beg, GCDR3end, CDR3DNA, CDR3AA, CDR3AALength, CDR3MW, CDR3pI = FindCDR3(AlignedSeq, IMGTVDJ, IgBLASTAnalysis[60], IgBLASTAnalysis[9], GJbeg, species, IDLength, JunctionLength) # DermVDJ, IMGTGermVDJ, Jgene1Name, Jbeg
						# IgBLASTAnalysis[51] = CDR3beg-1  # because IgBLAST made FW3 to end of V...
						CDR3beg = SeqFixed.find(CDR3DNA)+1
						CDR3end = CDR3beg + len(CDR3DNA)-1
						IgBLASTAnalysis[51] = CDR3beg-1

						now = time.strftime('%c')
					except:
						#print(SeqName)
						NotValid = True
						ErLog = SeqName + 'was problematic (line 766)\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

		if IgLine[0:7] == 'Query= ':
			# write current progress to file
			'''
			file_handle = open(progressBarFile,'w')
			progress = str(int(current_seq*100/totel_seq))
			file_handle.write(progress)
			file_handle.write(',' + str(current_seq) + '/' + str(totel_seq))
			file_handle.close()
			'''

			pct = int(current_seq * 100 / totel_seq)
			label = 'Processing: ' + str(current_seq) + '/' + str(totel_seq)
			signal.emit(pct, label)
			#print('Current Progress: ' + progress)
			current_seq += 1


			if SeqNumber >0:
				if NotValid == False:


					IgBLASTAnalysis.append(project)
					IgBLASTAnalysis.append(grouping)
					IgBLASTAnalysis.append(subgroup)
					IgBLASTAnalysis.append(species)

					Sequence = Sequences[IgBLASTAnalysis[0]]
					if Sequence != '':
						IgBLASTAnalysis.append(Sequence)
					else:
						NotValid = True
						ErLog = SeqName + 'has no sequence (line 210)\n'
						with open(ErlogFile, 'a') as currentfile:
							currentfile.write(ErLog)
						BadSeqs += 1

					if GeneType == 'Heavy':
						IsoSeq = (Sequence[(Jend):])
						#print(SeqName)
						IsoSeq = IsoSeq.strip('N')
						AGCTs = IsoSeq.count('A') + IsoSeq.count('G') +IsoSeq.count('C') +IsoSeq.count('T')
						if AGCTs > 5:  # todo decide if can determine isotype from < 5 or need more then
							#print('Start ' + SeqName + ' counts at ' + time.strftime('%c'))
							if species == 'Human':
								Isotype = VGenesSeq.CallIsotype(IsoSeq)
							elif species == 'Mouse':
								Isotype = VGenesSeq.CallIsotypeMouse(IsoSeq)
							else:
								Msg = 'Your current species is: ' + species + \
								      '\nWe do not support this species!'
								return Msg
							#print('end ' + SeqName + ' counts at ' + time.strftime('%c'))

						else:
							if len(IsoSeq) > 2:
								if IsoSeq[:3] == 'CCT' or IsoSeq == 'CTT':
									Isotype = 'IgG'
								elif IsoSeq[:3] == 'CAT':
									Isotype  = 'IgA'
								elif IsoSeq[:3] == 'GGA':
									Isotype  = 'IgM'
								elif IsoSeq[:3] == 'CAC':
									Isotype  = 'IgD'
								else:
									Isotype = IsoSeq
							else:
								Isotype = 'Unknown'
					else:
						if GeneType == 'Kappa':
							Isotype = 'Kappa'
						elif GeneType == 'Lambda':
							Isotype = 'Lambda'

					# TODO need code for following for now just added as blanks for database fields: can delete each as completed

					IgBLASTAnalysis.append(GermlineSeq)
					IgBLASTAnalysis.append(CDR3DNA)
					IgBLASTAnalysis.append(CDR3AA)
					IgBLASTAnalysis.append(CDR3AALength)
					IgBLASTAnalysis.append(CDR3beg)
					IgBLASTAnalysis.append(CDR3end)
					IgBLASTAnalysis.append('Specificity')
					IgBLASTAnalysis.append('Subspecificity')
					IgBLASTAnalysis.append(0)
					IgBLASTAnalysis.append(0)
					IgBLASTAnalysis.append(VGeneLocus)
					IgBLASTAnalysis.append(JgeneLocus)
					# if DgeneLocus == '':
					#     DgeneLocus = 'none'
					IgBLASTAnalysis.append(DgeneLocus)
					IgBLASTAnalysis.append(now)
					IgBLASTAnalysis.append(' ')
					IgBLASTAnalysis.append(' ')
					IgBLASTAnalysis.append(TotalMuts)
					IgBLASTAnalysis.append(MutListDB)
					IgBLASTAnalysis.append(IDevent)
					IgBLASTAnalysis.append(CDR3MW)
					IgBLASTAnalysis.append(CDR3pI)
					IgBLASTAnalysis.append(Isotype)
					IgBLASTAnalysis.append(GCDR3beg)
					IgBLASTAnalysis.append(GCDR3end)
					IgBLASTAnalysis.append('Blank6')
					IgBLASTAnalysis.append('Blank7')
					IgBLASTAnalysis.append('Blank8')
					IgBLASTAnalysis.append('Blank9')
					IgBLASTAnalysis.append('Blank10')
					IgBLASTAnalysis.append('Blank11')
					IgBLASTAnalysis.append('Blank12')
					IgBLASTAnalysis.append('Blank13')
					IgBLASTAnalysis.append('Blank14')
					IgBLASTAnalysis.append('Blank15')
					IgBLASTAnalysis.append('Blank16')
					IgBLASTAnalysis.append('Blank17')
					IgBLASTAnalysis.append('Blank18')
					IgBLASTAnalysis.append('Blank19')
					IgBLASTAnalysis.append('Blank20')
					# IgBLASTAnalysis.append('TheEnd')

					# IgBLASTset[SeqNumber-1] = IgBLASTAnalysis  #  if not first record it saves
					if project == 'ByFunction':
						if Importing == False:
							IgBLASTAnalysis[75] = GeneType
							if multiProject != '':
								IgBLASTAnalysis[75] = multiProject
							if Productive == "Yes":
								IgBLASTAnalysis[76] = 'Functional'
							else:
								IgBLASTAnalysis[76]  = 'Nonfunctional'

							if multiProject != '':
								IgBLASTAnalysis[76] = GeneType

						# grouping = Productive
						IgBLASTAnalysis[77] = ReadingFrame
						if multiProject != '':
							if Productive == "Yes":
								IgBLASTAnalysis[77] = 'Functional'
							else:
								IgBLASTAnalysis[77]  = 'Nonfunctional'


					IgBLASTset.append((IgBLASTAnalysis))


				IgBLASTAnalysis = []  # can now empty IgBLASTANalysis for next seq
				NotValid = False




			SeqNumber += 1
			LineCount = 0
			TotMut = 0
			GrabAlignment = False
			GrabRecomb = False
			GrabJunction = False

			# NotValid = False
			SeqAlignment = ''
			FR1hit = False
			CDR1hit = False
			FW2hit = False
			CDR2hit = False




			SeqName = IgLine[7:]
			SeqName = SeqName.strip()  # removes any spaces before or after
			if SeqName[(len(SeqName)-6):] == 'Import':
				Importing = True
				SeqNameSplit = Segments.split("_")
				SeqName = SeqNameSplit[0]
				project = SeqNameSplit[1]
				grouping = SeqNameSplit[2]



			IgBLASTAnalysis.append(SeqName)
		elif IgLine[0:7] == 'Length=':
			LenSeq = int(IgLine[7:])

			IgBLASTAnalysis.append(str(LenSeq))
		elif IgLine[0:25] == '***** No hits found *****':
			NotValid = True

			ErLog = SeqName + ': No hits found\n'
			with open(ErlogFile, 'a') as currentfile:
				currentfile.write(ErLog)
			BadSeqs += 1
		elif IgLine[0:21] == 'V-(D)-J rearrangement':  # trick to get next line based on line count
			GrabRecomb = True  #  trick to get next line based on line count
		elif GrabRecomb == True: #  spring this only if previous line was V-(D)-J...
			# Should make this a different function to clean up
			GrabRecomb = False
			Segments = IgLine

			RecombParts = Segments.split("\t")
			#  Generates a list of each item in the string
			# seperated by tabs or any character '/n' is para mark

			GeneType = ""
			for segment in RecombParts:
				if segment == 'VH':
					GeneType = 'Heavy'
					break
				elif segment == 'VK':
					GeneType = 'Kappa'
					break
				else:
					GeneType = 'Lambda'

			IgBLASTAnalysis.append(GeneType)


			Vgene1 = ''
			Vgene2 = ''
			Vgene3 = ''
			Dgene1 = ''
			Dgene2 = ''
			Dgene3 = ''
			Jgene1 = ''
			Jgene2 = ''
			Jgene3 = ''

			Vgene = RecombParts[0]
			VChoice = Vgene.split(',')
			i = 0
			for i in range(0, len(VChoice)):
				Vgene1 = VChoice[0]

				if len(VChoice) == 2:
					Vgene2 = VChoice[1]
				elif len(VChoice) == 3:
					Vgene3 = VChoice[2]
			if Vgene1 != 'N/A':
				IgBLASTAnalysis.append(Vgene1)
				IgBLASTAnalysis.append(Vgene2)
				IgBLASTAnalysis.append(Vgene3)
			else:
				NotValid = True

				ErLog = SeqName + ': No V gene found\n'
				with open(ErlogFile, 'a') as currentfile:
					currentfile.write(ErLog)
				BadSeqs += 1

			IndexN = 0

			if GeneType == 'Heavy':
				IndexN = 1
			Dgene = RecombParts[IndexN]
			DChoice = Dgene.split(',')
			i = 0
			for i in range(0, len(DChoice)):
				Dgene1 = DChoice[0]
				if len(DChoice) == 2:
					Dgene2 = DChoice[1]
				elif len(DChoice) == 3:
					Dgene3 = DChoice[2]

			if GeneType == 'Heavy':  #only ned insert D for H chains
				IgBLASTAnalysis.append(Dgene1)
				IgBLASTAnalysis.append(Dgene2)
				IgBLASTAnalysis.append(Dgene3)
			else:
				IgBLASTAnalysis.append('None')
				IgBLASTAnalysis.append('None')
				IgBLASTAnalysis.append('None')


			IndexN += 1

			Jgene = RecombParts[IndexN]
			JChoice = Jgene.split(',')
			i = 0
			for i in range(0, len(JChoice)):
				Jgene1 = JChoice[0]
				if len(JChoice) == 2:
					Jgene2 = JChoice[1]
				elif len(JChoice) == 3:
					Jgene3 = JChoice[2]

			if Jgene1 != 'N/A':
				IgBLASTAnalysis.append(Jgene1)
				IgBLASTAnalysis.append(Jgene2)
				IgBLASTAnalysis.append(Jgene3)
			else:
				if GetProductive == 0:
					NotValid = True

				ErLog = SeqName + ': No J gene found\n'
				with open(ErlogFile, 'a') as currentfile:
					currentfile.write(ErLog)
				BadSeqs += 1

			IndexN += 2
			StopCodon = RecombParts[IndexN]
			IgBLASTAnalysis.append(StopCodon)

			IndexN += 1
			ReadingFrame = RecombParts[IndexN]
			IgBLASTAnalysis.append(ReadingFrame)

			IndexN += 1
			Productive = RecombParts[IndexN]
			#print(Productive)
			if GetProductive == 0:
				if Productive == 'Yes':
					IgBLASTAnalysis.append(Productive)
				else:
					NotValid = True
					ErLog = SeqName + ' was not a productive rearrangement\n'
					with open(ErlogFile, 'a') as currentfile:
						currentfile.write(ErLog)
					BadSeqs += 1
			else:
				IgBLASTAnalysis.append(Productive)

			IndexN += 1
			Strand = RecombParts[IndexN]
			IgBLASTAnalysis.append(Strand)
		elif IgLine[0:16] == 'V-(D)-J junction':  # trick to get next line based on line count
			GrabJunction = True  #  trick to get next line based on line count
			GrabIt = LineCount
		elif GrabJunction == True:
			GrabJunction = False
			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark
			VSeqend = RecombParts[0]
			IgBLASTAnalysis.append(VSeqend)

			if GeneType == "Heavy":
				VDJunction = RecombParts[1]
				Dregion = RecombParts[2]
				DJJunction = RecombParts[3]
				begJ = RecombParts[4]
				VJunction = ''  # only for Light chains but need to add field for database integrity

				IgBLASTAnalysis.append(VDJunction)
				IgBLASTAnalysis.append(Dregion)
				IgBLASTAnalysis.append(DJJunction)
				IgBLASTAnalysis.append(begJ)
				IgBLASTAnalysis.append(VJunction) # only for Light chains but need to add field for database integrity

			else: #  TODO need to set up for light chain so not out of whack with fewer fields

				VJunction = RecombParts[1]
				begJ = RecombParts[2]
				VDJunction = ''  # only for heavy chains but need to add field for database integrity
				Dregion = ''  # only for heavy chains but need to add field for database integrity
				DJJunction = ''  # only for heavy chains but need to add field for database integrity

				IgBLASTAnalysis.append(VDJunction)  # only for heavy chains but need to add field for database integrity
				IgBLASTAnalysis.append(Dregion)  # only for heavy chains but need to add field for database integrity
				IgBLASTAnalysis.append(DJJunction)  # only for heavy chains but need to add field for database integrity

				IgBLASTAnalysis.append(begJ)
				IgBLASTAnalysis.append(VJunction)
		elif IgLine[0:3] == 'FR1':
			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR1From = RecombParts[1]
			FR1To = RecombParts[2]
			FR1length = RecombParts[3]
			FR1matches = RecombParts[4]
			FR1mis = RecombParts[5]
			FR1gaps = RecombParts[6]
			FR1PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR1From)
			IgBLASTAnalysis.append(FR1To)
			IgBLASTAnalysis.append(FR1length)
			IgBLASTAnalysis.append(FR1matches)
			IgBLASTAnalysis.append(FR1mis)
			IgBLASTAnalysis.append(FR1gaps)
			IgBLASTAnalysis.append(FR1PercentIdentity)

			FR1hit = True



			TotMut += int(FR1mis)
			TotMut += int(FR1gaps)
		elif IgLine[0:4] == 'CDR1':
			if FR1hit == False:
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine

			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR1From = RecombParts[1]
			CDR1To = RecombParts[2]
			CDR1length = RecombParts[3]
			CDR1matches = RecombParts[4]
			CDR1mis = RecombParts[5]
			CDR1gaps = RecombParts[6]
			CDR1PercentIdentity = RecombParts[7]


			IgBLASTAnalysis.append(CDR1From)
			IgBLASTAnalysis.append(CDR1To)
			IgBLASTAnalysis.append(CDR1length)
			IgBLASTAnalysis.append(CDR1matches)
			IgBLASTAnalysis.append(CDR1mis)
			IgBLASTAnalysis.append(CDR1gaps)
			IgBLASTAnalysis.append(CDR1PercentIdentity)

			CDR1hit = True

			try:
				TotMut += int(CDR1mis)
			except:
				TotMut += 0

			try:
				TotMut += int(CDR1gaps)
			except:
				TotMut += 0
		elif IgLine[0:3] == 'FR2':
			if CDR1hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR2From = RecombParts[1]
			FR2To = RecombParts[2]
			FR2length = RecombParts[3]
			FR2matches = RecombParts[4]
			FR2mis = RecombParts[5]
			FR2gaps = RecombParts[6]
			FR2PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR2From)
			IgBLASTAnalysis.append(FR2To)
			IgBLASTAnalysis.append(FR2length)
			IgBLASTAnalysis.append(FR2matches)
			IgBLASTAnalysis.append(FR2mis)
			IgBLASTAnalysis.append(FR2gaps)
			IgBLASTAnalysis.append(FR2PercentIdentity)

			FW2hit = True

			TotMut += int(FR2mis)
			TotMut += int(FR2gaps)
		elif IgLine[0:4] == 'CDR2':
			if FW2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR2From = RecombParts[1]
			CDR2To = RecombParts[2]
			CDR2length = RecombParts[3]
			CDR2matches = RecombParts[4]
			CDR2mis = RecombParts[5]
			CDR2gaps = RecombParts[6]
			CDR2PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(CDR2From)
			IgBLASTAnalysis.append(CDR2To)
			IgBLASTAnalysis.append(CDR2length)
			IgBLASTAnalysis.append(CDR2matches)
			IgBLASTAnalysis.append(CDR2mis)
			IgBLASTAnalysis.append(CDR2gaps)
			IgBLASTAnalysis.append(CDR2PercentIdentity)

			CDR2hit = True

			if CDR2hit != 'N/A' and CDR2mis != 'N/A':
				try:
					TotMut += int(CDR2mis)
					TotMut += int(CDR2gaps)
				except:
					print('oops')
		elif IgLine[0:3] == 'FR3':
			if CDR2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #CDR2
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)


			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR3From = RecombParts[1]
			FR3To = RecombParts[2]
			FR3length = RecombParts[3]
			FR3matches = RecombParts[4]
			FR3mis = RecombParts[5]
			FR3gaps = RecombParts[6]
			FR3PercentIdentity = RecombParts[7]

			IgBLASTAnalysis.append(FR3From)
			IgBLASTAnalysis.append(FR3To)
			IgBLASTAnalysis.append(FR3length)
			IgBLASTAnalysis.append(FR3matches)
			IgBLASTAnalysis.append(FR3mis)
			IgBLASTAnalysis.append(FR3gaps)
			IgBLASTAnalysis.append(FR3PercentIdentity)
			try:
				if FR3mis == 'N/A':
					FR3mis = 0
				if FR3gaps == 'N/A':
					FR3gaps = 0
				TotMut += int(FR3mis)
				TotMut += int(FR3gaps)
			except:
				print('mistake')

			IgBLASTAnalysis.append(str(TotMut))
		# elif IgLine[0:4] == 'CDR3':

		# SeqNumber += 1

		# TODO maybe need to analyze IgBLAST output for SHM

	conn.close()
		# os.remove(workingfilename)

def parseIgBlast(Files, datalist):
	pass

def parseIMGT(Files, datalist):
	pass

def GetGLCDRs(Sequence, species):
	import os


	#todo change to app folder
	try:
		#DBpathname = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'VDJGenes.db')
		DBpathname = os.path.join(working_prefix, 'Data', 'VDJGenes.db')


		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)


		conn = db.connect(DBpathname)


	except:
		DBpathname = '/Volumes/Promise Pegasus/Dropbox/VGenes/VDJGenes.db'

		(dirname, filename) = os.path.split(DBpathname)
		os.chdir(dirname)


		conn = db.connect(DBpathname)

	IgBLASTAnalysis = []


	workingdir = os.path.join(working_prefix, 'IgBlast')  #'/Applications/IgBlast/database'
	workingfilename = os.path.join(working_prefix, 'IgBlast', 'WorkingFile.nt') #'/Applications/IgBlast/database/WorkingFile.nt'
	os.chdir(workingdir)

	Sequence = '>It\n' + Sequence + '\n'



	with open(workingfilename, 'w') as currentfile:
		currentfile.write(Sequence)


	if species == 'Human':
		BLASTCommandLine = igblast_path + " -germline_db_V IG/Human/HumanVGenes.nt -germline_db_J IG/Human/HumanJGenes.nt -germline_db_D IG/Human/HumanDGenes.nt -organism human -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt 3"
		IgBlastOut = os.popen(BLASTCommandLine)
	elif species == 'Mouse':
		BLASTCommandLine = igblast_path + " -germline_db_V IG/Mouse/MouseVGenes.nt -germline_db_J IG/Mouse/MouseJGenes.nt -germline_db_D IG/Mouse/MouseDGenes.nt -organism mouse -domain_system kabat -query WorkingFile.nt -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3"
		IgBlastOut = os.popen(BLASTCommandLine)


	for IgLine in IgBlastOut:


		if IgLine[0:3] == 'FR1':
			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR1From = RecombParts[1]
			FR1To = RecombParts[2]


			IgBLASTAnalysis.append(FR1From)
			IgBLASTAnalysis.append(FR1To)

			FR1hit = True


		elif IgLine[0:4] == 'CDR1':
			if FR1hit == False:
				IgBLASTAnalysis.append(0)
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR1From = RecombParts[1]
			CDR1To = RecombParts[2]

			IgBLASTAnalysis.append(CDR1From)
			IgBLASTAnalysis.append(CDR1To)
			CDR1hit = True

		elif IgLine[0:3] == 'FR2':
			if CDR1hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR2From = RecombParts[1]
			FR2To = RecombParts[2]

			IgBLASTAnalysis.append(FR2From)
			IgBLASTAnalysis.append(FR2To)

			FW2hit = True

		elif IgLine[0:4] == 'CDR2':
			if FW2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			CDR2From = RecombParts[1]
			CDR2To = RecombParts[2]

			IgBLASTAnalysis.append(CDR2From)
			IgBLASTAnalysis.append(CDR2To)

			CDR2hit = True

		elif IgLine[0:3] == 'FR3':
			if CDR2hit == False:
				IgBLASTAnalysis.append(0) #FW1
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0) #CDR1
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #FW2
				IgBLASTAnalysis.append(0)

				IgBLASTAnalysis.append(0)  #CDR2
				IgBLASTAnalysis.append(0)

			Junction = IgLine
			#print(Junction)
			RecombParts = Junction.split("\t")
			#  Generates a list of each item in the string seperated by tabs or any character '/n' is para mark

			FR3From = RecombParts[1]
			FR3To = RecombParts[2]

			IgBLASTAnalysis.append(FR3From)
			IgBLASTAnalysis.append(FR3To)

	conn.close()

	return IgBLASTAnalysis

def ParseAlignment(Alignment):
	global LineLen
	AlignmentContents = []
	LineNum = 0
	LineInc = 0
	FirstOne = True
	FirstQ = True
	PartsS = ''
	line = ''
	Vend = 0
	PartsS = Alignment
	PartsS += '\nZ                             '  #to catch any without a true D
	Parts = PartsS.split('\n')
	SeqOrient = "forward"
	global element, endelement
	try:
		for line in Parts:
			if len(line) > 25:
				Tester = line[0:27]
				Tester = Tester.strip(' ')
				if Tester[0:2] == 'Qu':
					LineNum += 1

					if FirstQ == True:
						LineLen = len(line)
						StartParse = len(line) -100
						EndParse = StartParse+3
						element = line[StartParse:EndParse]
						element = element.strip()
						SeqBegin = int(element)
						FirstQ = False
						#AlignmentContents.append(SeqBegin)


				elif Tester[0:5] == 'lcl|Q':
					LineNum += 1

					if FirstQ == True:
						LineLen = len(line)
						StartParse = len(line) -100
						EndParse = StartParse+3
						try:
							element = line[StartParse:EndParse]
						except:
							return 'short'
						element = element.strip()
						SeqBegin = int(element)
						FirstQ = False
						#AlignmentContents.append(SeqBegin)
					#DidntWork = 'None'
				   # return DidntWork
					SeqOrient = "reversed"

				elif line[0] == 'V':


					# todo damn Ig Blast does different spacing so current doesn't work...maybe can split by space or read an iterate number
					if LineNum != LineInc:   # it's the first V in a section

						if FirstOne == True:
							LineLen = len(line)
							StartParse = len(line) -100
							EndParse = StartParse+3
							element = line[StartParse:EndParse]
							element = element.strip()
							Vbeg = int(element)
							AlignmentContents.append(Vbeg)
							if Vbeg > SeqBegin:
								AdjustEnd = 0 -(Vbeg - SeqBegin)
							elif SeqBegin > Vbeg:
								AdjustEnd  = 0#SeqBegin-Vbeg
							else:
								AdjustEnd = 0

							FirstOne = False

						elif FirstOne == False:
							StartParse = len(line) -3
							EndParse = StartParse+3
							element = line[StartParse:EndParse]
						LineInc = LineNum
				elif line[0] == 'D' or line[0] == 'J' or line[0] == 'Z':
					element = element.strip()
					Vend = int(element) #+ AdjustEnd

					AlignmentContents.append(Vend)
					break

		LineNum = 0
		LineInc = 0
		FirstOne = True
		for line in Parts:
			if len(line) > 25:
				Tester = line[0:27]
				Tester = Tester.strip(' ')
				if Tester[0:2] == 'Qu':
					LineNum += 1
					# FirstOne = True
				elif Tester[0:5] == 'lcl|Q':
					LineNum += 1
					# FirstOne = True
				elif line[0] == 'D':
					NewLine = len(line)
					ShortLine = LineLen - NewLine
					if LineNum != LineInc:   # it's the first D in a section


						if FirstOne == True:
							StartParse = len(line) -(100-ShortLine)
							EndParse = StartParse+4
							element = line[StartParse:EndParse]
							element = element.strip()
							Vbeg = int(element)
							FirstOne = False
							AlignmentContents.append(Vbeg)  #writes the beginning
							# need to grab if to short
							StartParse = len(line) -3
							EndParse = StartParse+3
							endelement = line[StartParse:EndParse]  #grabs it but doesn't write it yet
						elif FirstOne == False:
							# code to ensure still same D not new one start next section when two Ds
							StartParse2 = len(line) -(100-ShortLine)
							EndParse2 = StartParse2+4
							element2 = line[StartParse2:EndParse2]
							element2 = element2.strip()
							test = int(element2)
							test2 = endelement.strip()
							if test == int(test2)+1:
								StartParse = len(line) -3
								EndParse = StartParse+3             #as it continues reincrements end element
								endelement = line[StartParse:EndParse]
						LineInc = LineNum
				elif line[0] == 'J' or line[0] == 'Z':

					if FirstOne == True: # no D was detected
						AlignmentContents.append(0)
						AlignmentContents.append(0)
					else:
						endelement = endelement.strip()
						Vend = int(endelement)
						AlignmentContents.append(Vend)

					break

		LineNum = 0
		LineInc = 0
		FirstOne = True
		for line in Parts:
			if len(line) > 25:
				Tester = line[0:27]
				Tester = Tester.strip(' ')
				if Tester[0:2] == 'Qu':
					LineNum += 1
					# FirstOne = True
				elif Tester[0:5] == 'lcl|Q':
					LineNum += 1
					# FirstOne = True
				elif line[0] == 'J':
					NewLine = len(line)
					ShortLine = LineLen - NewLine
					if LineNum != LineInc:   # it's the first V in a section


						if FirstOne == True:
							StartParse = len(line) -(100-ShortLine)
							EndParse = StartParse+4
							element = line[StartParse:EndParse]
							element = element.strip()
							Vbeg = int(element)
							FirstOne = False
							AlignmentContents.append(Vbeg)
							# need to grab if to short
							StartParse = len(line) -3
							EndParse = StartParse+3
							endelement = line[StartParse:EndParse]
						elif FirstOne == False:
							StartParse = len(line) -3
							EndParse = StartParse+3
							endelement = line[StartParse:EndParse]
						LineInc = LineNum
				elif line[0] == 'Z':
					endelement = endelement.strip()
					Vend = int(endelement)

					if FirstOne == True: # no J was detected
						AlignmentContents.append(0)
						AlignmentContents.append(0)
					else:
						AlignmentContents.append(Vend)

					AlignmentContents.append(AdjustEnd)
					AlignmentContents.append(SeqOrient)
					AlignmentContents.append(SeqBegin)

					return AlignmentContents
	except:
		DidntWork = 'None'
		return DidntWork

def FindCDR3(SeqFixed, IMGTVDJ, Vend, JGeneName, GJbeg, species, IDSeqLen, JunctionLength):


	if JGeneName == 'N/A':
		return 0, 0, 'N/A', 'N/A', 0, 0, 0

	spacers = IMGTVDJ.count('.')
	Vend = Vend + spacers
	for i in range(0,spacers):
		SeqFixed = '.' + SeqFixed

	CDR3beg = 312

	if IDSeqLen>0: CDR3beg += IDSeqLen

	VPart = Vend - CDR3beg  #'codon 108'



	if species == 'Human':
		JCDR3 = JHuman[JGeneName]
	else:
		JCDR3 = JMouse[JGeneName]

	if GJbeg < JCDR3:
		JPart = JCDR3 - GJbeg
	else:
		JPart = 0

	CDR3Length  = VPart + JunctionLength + JPart

	while CDR3Length % 3 != 0:    # it's not a full codon...either J is trimmed or OOF
		CDR3Length +=1

	CDR3end = CDR3beg + CDR3Length + IDSeqLen
	CDR3 = SeqFixed[CDR3beg:CDR3end]
	CDR3AA1, ErMes = VGenesSeq.Translator(CDR3,0)

	CDR3AA = str(CDR3AA1)

	CDR3AALength = len(CDR3AA)
	CDR3beg -= spacers
	CDR3end -= spacers
	CDR3MW = VGenesSeq.OtherParam(CDR3AA, 'AAMW',0, True)
	CDR3pI = VGenesSeq.OtherParam(CDR3AA, 'AApI',0, True)


	return CDR3beg, CDR3end, CDR3, CDR3AA, CDR3AALength, CDR3MW, CDR3pI

# CDR3DNA text, CDR3AA text, CDR3Length text, CDR3beg text, CDR3end text,


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
#
# """
# The CDR3 is delimited by (but does not include) the anchor positions 2nd-CYS 104 and J-PHE or J-TRP 118.
# As for the CDR1 and CDR2 these anchor positions (that belong to the neighbouring FR) are shown as
# squares in the IMGT Collier de Perles.
# The JUNCTION includes 2nd-CYS 104 and J-PHE or J-TRP 118 and is therefore two amino acids longer than the CDR3.
# The CDR3 numbering goes from 105 to 117 and, if necessary, gaps or additional positions are added at the top
# of the loop http://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
# Note that:
# the J-PHE or J-TRP belongs to the characteristic J-REGION motif 'F/W-G-X-G' at positions 118-121
# the CDR3 is delimited by the same anchor positions (2nd-CYS 104 and J-PHE or J-TRP 118), whatever the
# receptor type (IG or TR), the chain type (heavy or light for IG; alpha, beta, gamma or delta for TR) or the species.
#
# Steps:
# 1. build germline seq as GL from the IMGTGermV with spacing, my genes junction, and GL J
# 2. determine length difference when IMGT spacing included by including ...'s or not an getting seq length
#     -note, need to identify if ins/del present at this point as when dots removed from IMGT V should be same length
#     if not same compare by base to find ID event or maybe it's in IgBLAST report.
#     ---maybe store ID events as start position (last nuc before event) and length (+/-). then can use for aligments. (23,3) or (23, -6)
# 3. then calculate which codon is last V Cys based on seq lenght differences and verify really a Cys. CDR3 starts on next full codon.
# 4. because few Js can use if, elif statement to find particular J (can group and maybe even use J locus)...compare
# where it should be to Jbeg. It's the last full codon before the F/W-G-X-G or if that's deleted, before the next J codon.
# 5. Then can grab end of V (from Cys), junction (from IgBLAST) and beginning of J and consrtuct.
# """



# todo need codeto find second d...can use ifo above to figure it out
