__author__ = 'wilsonp'

from VGenesDialogues import openFile, openFiles, newFile, saveFile, questionMessage, informationMessage, setItem, \
	setText
from collections import Counter
import VGenesSeq
import VGenesSQL
from operator import itemgetter
import itertools
import sys

def MapHotspots(self, ReportType, DBFilename, SeqName):


	errrep = '\nThe following errors occured during processing:\n'
	DataIn = []
	#     scores = (Nuc, MutType, CTType, RScore, MutIndex, TENDA, TENDG, TENDC, TENDT)
	if ReportType != 'AnalysisDB Heat Map':
		fields = ['SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To', 'FR2To', 'CDR2To', 'FR3To',
				  'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend']

		SQLStatement = VGenesSQL.MakeSQLStatement(self, fields, SeqName)
		DataIn = VGenesSQL.RunSQL(DBFilename, SQLStatement)  # returns list of tuples where seqname is first
	else:
		# IMGT_ID, SeqName, Species, functionality, IMGTSequence, FWR1, CDR1, FWR2, CDR2, FWR3, ID

		DataIn = VGenesSQL.RunAnalysisSQL(DBFilename)
		NewMap = []
		# for item in AnalysisDB:
		#     record  = []
		#     record.append(item[1])
		#     IMGTSequence = item[4]

	if ReportType == 'Standard':
		Document = 'Name, Sequence, Region, Targeting, Type, Conservation, Mutation Index, Tendency-A, Tendency-G, Tendency-C, Tendency-T, Germline, Sequece, Targeting, Type, Conservation, Mutation Index, Tendency-A, Tendency-G, Tendency-C, Tendency-T\n'
	elif ReportType == 'Hotspots' or ReportType == 'HS-Summary' or ReportType == 'Combined' or ReportType == 'Heat Map'  or ReportType == 'GC-Heat Map' or ReportType == 'AnalysisDB Heat Map':
		Header = []
		FW1 = []
		CDR1 = []
		FW2 = []
		CDR2 = []
		FW3 = []
		CDR3 = []
		FW4 = []
		GFW1 = []
		GCDR1 = []
		GFW2 = []
		GCDR2 = []
		GFW3 = []
		GCDR3 = []
		GFW4 = []

		F1GermSet = []
		F1SeqSet = []
		C1GermSet = []
		C1SeqSet = []
		F2GermSet = []
		F2SeqSet = []
		C2GermSet = []
		C2SeqSet = []
		F3GermSet = []
		F3SeqSet = []
		C3GermSet = []
		C3SeqSet = []
		F4GermSet = []
		F4SeqSet = []
		counts = []
		Item = ()
		FW1len = 0
		CDR1len = 0
		FW2len = 0
		CDR2len = 0
		FW3len = 0
		CDR3len = 0
		FW4len = 0
		mutcount = 0
		Document = ''
		SumCounts = []
		SumTup = (
		'Name', 'Mutations', 'Total WRCY', 'Total WRC', 'Total WA', 'Total TW', 'Total AID-CS', 'Total Neutral'
		, 'FW WRCY', 'FW WRC', 'FW WA', 'FW TW', 'FW AID-CS', 'FW Neutral'
		, 'CDR WRCY', 'CDR WRC', 'CDR WA', 'CDR TW', 'CDR AID-CS', 'CDR Neutral'
		, 'FW1 WRCY', 'FW1 WRC', 'FW1 WA', 'FW1 TW', 'FW1 AID-CS', 'FW1 Neutral'
		, 'CDR1 WRCY', 'CDR1 WRC', 'CDR1 WA', 'CDR1 TW', 'CDR1 AID-CS', 'CDR1 Neutral'
		, 'FW2 WRCY', 'FW2 WRC', 'FW2 WA', 'FW2 TW', 'FW2 AID-CS', 'FW2 Neutral'
		, 'CDR2 WRCY', 'CDR2 WRC', 'CDR2 WA', 'CDR2 TW', 'CDR2 AID-CS', 'CDR2 Neutral'
		, 'FW3 WRCY', 'FW3 WRC', 'FW3 WA', 'FW3 TW', 'FW3 AID-CS', 'FW3 Neutral'
		, 'CDR3 WRCY', 'CDR3 WRC', 'CDR3 WA', 'CDR3 TW', 'CDR3 AID-CS', 'CDR3 Neutral'
		, 'FW4 WRCY', 'FW4 WRC', 'FW4 WA', 'FW4 TW', 'FW4 AID-CS', 'FW4 Neutral'
		, 'GTotal WRCY', 'GTotal WRC', 'GTotal WA', 'GTotal TW', 'GTotal AID-CS', 'GTotal Neutral'
		, 'GFW WRCY', 'GFW WRC', 'GFW WA', 'GFW TW', 'GFW AID-CS', 'GFW Neutral'
		, 'GCDR WRCY', 'GCDR WRC', 'GCDR WA', 'GCDR TW', 'GCDR AID-CS', 'GCDR Neutral'
		, 'GFW1 WRCY', 'GFW1 WRC', 'GFW1 WA', 'GFW1 TW', 'GFW1 AID-CS', 'GFW1 Neutral'
		, 'GCDR1 WRCY', 'GCDR1 WRC', 'GCDR1 WA', 'GCDR1 TW', 'GCDR1 AID-CS', 'GCDR1 Neutral'
		, 'GFW2 WRCY', 'GFW2 WRC', 'GFW2 WA', 'GFW2 TW', 'GFW2 AID-CS', 'GFW2 Neutral'
		, 'GCDR2 WRCY', 'GCDR2 WRC', 'GCDR2 WA', 'GCDR2 TW', 'GCDR2 AID-CS', 'GCDR2 Neutral'
		, 'GFW3 WRCY', 'GFW3 WRC', 'GFW3 WA', 'GFW3 TW', 'GFW3 AID-CS', 'GFW3 Neutral'
		, 'GCDR3 WRCY', 'GCDR3 WRC', 'GCDR3 WA', 'GCDR3 TW', 'GCDR3 AID-CS', 'GCDR3 Neutral'
		, 'GFW4 WRCY', 'GFW4 WRC', 'GFW4 WA', 'GFW4 TW', 'GFW4 AID-CS', 'GFW4 Neutral')

		Summary = []
		Summary.append(SumTup)

	DataIn.sort(key=itemgetter(1))
	if ReportType == 'Mutation Frequencies':
		NumScored = 0
		Smuts = 0
		Rmuts = 0
		SMutations = {'AWA': 0, 'ANeutral': 0, 'TTW': 0, 'TNeutral': 0, 'CWRCY': 0
			, 'CWRC': 0, 'CAID-CS': 0, 'CNeutral': 0, 'GWRCY': 0, 'GWRC': 0, 'GAID-CS': 0, 'GNeutral': 0
			, 'ATWA': 0, 'ATNeutral': 0, 'AGWA': 0, 'AGNeutral': 0, 'ACWA': 0, 'ACNeutral': 0
			, 'TATW': 0, 'TANeutral': 0, 'TGTW': 0, 'TGNeutral': 0, 'TCTW': 0, 'TCNeutral': 0
			, 'CAWRCY': 0, 'CAWRC': 0, 'CAAID-CS': 0, 'CANeutral': 0
			, 'CTWRCY': 0, 'CTWRC': 0, 'CTAID-CS': 0, 'CTNeutral': 0
			, 'CGWRCY': 0, 'CGWRC': 0, 'CGAID-CS': 0, 'CGNeutral': 0
			, 'GAWRCY': 0, 'GAWRC': 0, 'GAAID-CS': 0, 'GANeutral': 0
			, 'GTWRCY': 0, 'GTWRC': 0, 'GTAID-CS': 0, 'GTNeutral': 0
			, 'GCWRCY': 0, 'GCWRC': 0, 'GCAID-CS': 0, 'GCNeutral': 0}

		RMutations = {'AWA': 0, 'ANeutral': 0, 'TTW': 0, 'TNeutral': 0, 'CWRCY': 0
			, 'CWRC': 0, 'CAID-CS': 0, 'CNeutral': 0, 'GWRCY': 0, 'GWRC': 0, 'GAID-CS': 0, 'GNeutral': 0
			, 'ATWA': 0, 'ATNeutral': 0, 'AGWA': 0, 'AGNeutral': 0, 'ACWA': 0, 'ACNeutral': 0
			, 'TATW': 0, 'TANeutral': 0, 'TGTW': 0, 'TGNeutral': 0, 'TCTW': 0, 'TCNeutral': 0
			, 'CAWRCY': 0, 'CAWRC': 0, 'CAAID-CS': 0, 'CANeutral': 0
			, 'CTWRCY': 0, 'CTWRC': 0, 'CTAID-CS': 0, 'CTNeutral': 0
			, 'CGWRCY': 0, 'CGWRC': 0, 'CGAID-CS': 0, 'CGNeutral': 0
			, 'GAWRCY': 0, 'GAWRC': 0, 'GAAID-CS': 0, 'GANeutral': 0
			, 'GTWRCY': 0, 'GTWRC': 0, 'GTAID-CS': 0, 'GTNeutral': 0
			, 'GCWRCY': 0, 'GCWRC': 0, 'GCAID-CS': 0, 'GCNeutral': 0}

		SPositions = {'AWA': 0, 'ANeutral': 0, 'TTW': 0, 'TNeutral': 0, 'CWRCY': 0
			, 'CWRC': 0, 'CAID-CS': 0, 'CNeutral': 0, 'GWRCY': 0, 'GWRC': 0, 'GAID-CS': 0, 'GNeutral': 0
			, 'ATWA': 0, 'ATNeutral': 0, 'AGWA': 0, 'AGNeutral': 0, 'ACWA': 0, 'ACNeutral': 0
			, 'TATW': 0, 'TANeutral': 0, 'TGTW': 0, 'TGNeutral': 0, 'TCTW': 0, 'TCNeutral': 0
			, 'CAWRCY': 0, 'CAWRC': 0, 'CAAID-CS': 0, 'CANeutral': 0
			, 'CTWRCY': 0, 'CTWRC': 0, 'CTAID-CS': 0, 'CTNeutral': 0
			, 'CGWRCY': 0, 'CGWRC': 0, 'CGAID-CS': 0, 'CGNeutral': 0
			, 'GAWRCY': 0, 'GAWRC': 0, 'GAAID-CS': 0, 'GANeutral': 0
			, 'GTWRCY': 0, 'GTWRC': 0, 'GTAID-CS': 0, 'GTNeutral': 0
			, 'GCWRCY': 0, 'GCWRC': 0, 'GCAID-CS': 0, 'GCNeutral': 0}

		RPositions = {'AWA': 0, 'ANeutral': 0, 'TTW': 0, 'TNeutral': 0, 'CWRCY': 0
			, 'CWRC': 0, 'CAID-CS': 0, 'CNeutral': 0, 'GWRCY': 0, 'GWRC': 0, 'GAID-CS': 0, 'GNeutral': 0
			, 'ATWA': 0, 'ATNeutral': 0, 'AGWA': 0, 'AGNeutral': 0, 'ACWA': 0, 'ACNeutral': 0
			, 'TATW': 0, 'TANeutral': 0, 'TGTW': 0, 'TGNeutral': 0, 'TCTW': 0, 'TCNeutral': 0
			, 'CAWRCY': 0, 'CAWRC': 0, 'CAAID-CS': 0, 'CANeutral': 0
			, 'CTWRCY': 0, 'CTWRC': 0, 'CTAID-CS': 0, 'CTNeutral': 0
			, 'CGWRCY': 0, 'CGWRC': 0, 'CGAID-CS': 0, 'CGNeutral': 0
			, 'GAWRCY': 0, 'GAWRC': 0, 'GAAID-CS': 0, 'GANeutral': 0
			, 'GTWRCY': 0, 'GTWRC': 0, 'GTAID-CS': 0, 'GTNeutral': 0
			, 'GCWRCY': 0, 'GCWRC': 0, 'GCAID-CS': 0, 'GCNeutral': 0}

	if ReportType == 'Heat Map' or ReportType == 'GC-Heat Map' or ReportType == 'AnalysisDB Heat Map':
		for sequence in DataIn:
			#     forfields  = ['SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To','FR2To', 'CDR2To','FR3To', 'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend']
			Vbeg = int(sequence[12])
			if FW1len < int(sequence[3]) + Vbeg:
				FW1len = int(sequence[3]) + Vbeg
			if FW2len < int(sequence[5]) - int(sequence[4]):
				FW2len = int(sequence[5]) - int(sequence[4])
			if FW3len < int(sequence[7]) - int(sequence[6]):
				FW3len = int(sequence[7]) - int(sequence[6])
			if FW4len < int(sequence[9]) - int(sequence[8]):
				FW4len = int(sequence[9]) - int(sequence[8])
			if CDR1len < int(sequence[4]) - int(sequence[3]):
				CDR1len = int(sequence[4]) - int(sequence[3])
			if CDR2len < int(sequence[6]) - int(sequence[5]):
				CDR2len = int(sequence[6]) - int(sequence[5])
			if CDR3len < int(sequence[8]) - int(sequence[7]):
				CDR3len = int(sequence[8]) - int(sequence[7])

		FW1.clear()
		GFW1.clear()
		CDR1.clear()
		GCDR1.clear()
		FW2.clear()
		GFW2.clear()
		CDR2.clear()
		GCDR2.clear()
		FW3.clear()
		GFW3.clear()
		CDR3.clear()
		GCDR3.clear()
		FW4.clear()
		GFW4.clear()

	for sequence in DataIn:
		mutcount = 0
		SeqName = sequence[0]
		try:
			Seq = sequence[1]
			if ReportType == 'AnalysisDB Heat Map':
				IMGTSeq = sequence[1]
			Seq = Seq.replace('.', '')
			Seq = Seq[int(sequence[12]) - 1:]
			Seq = Seq.upper()
			GSeq = sequence[2]
			GSeq = GSeq[int(sequence[13]) - 1:]
			GSeq = GSeq.upper()
			SeqSet = list(set(Seq))
			for chr in SeqSet:
				if chr != "A" and chr != "G" and chr != "C" and chr != "T" and chr != "N":
					Seq = Seq.replace(chr, '')
			SeqSet = list(set(GSeq))
			for chr in SeqSet:
				if chr != "A" and chr != "G" and chr != "C" and chr != "T" and chr != "N":
					GSeq = GSeq.replace(chr, '')

			if sequence[10] == 'Insertion' or 'Deletion':
				mutations = sequence[11]
				MutEvents = mutations.split(',')
				for event in MutEvents:
					insert = ''
					Event = event.split('-')
					if Event[0] == "Insertion":
						for i in range(0, len(Event[2])):
							insert += '.'
							mutcount += 1
						GSeq = GSeq[0:int(Event[1])] + insert + GSeq[int(Event[1]):]
					elif Event[0] == 'Deletion':
						for i in range(0, len(Event[2])):
							insert += '.'
							mutcount += 1
						Seq = Seq[0:int(Event[1])] + insert + Seq[int(Event[1]):]

			Seq = Seq[:int(sequence[9])]
			GSeq = GSeq[:int(sequence[9])]

			if len(Seq) > len(GSeq):
				CheckTo = len(GSeq)
			else:
				CheckTo = len(Seq)

			Seq = Seq.upper()
			GSeq = GSeq.upper()
			for i in range(0, CheckTo):
				if Seq[i] != GSeq[i]:

					if (Seq[i] == 'G' or Seq[i] == 'A' or Seq[i] == 'C' or Seq[i] == 'T') and (
									GSeq[i] == 'G' or GSeq[i] == 'A' or GSeq[i] == 'C' or GSeq[i] == 'T'):
						mutcount += 1
					elif Seq[i] == 'N' and (GSeq[i] == 'G' or GSeq[i] == 'A' or GSeq[i] == 'C' or GSeq[i] == 'T'):

						Seq = Seq[:i] + GSeq[i] + Seq[i + 1:]

			# Nuc, MutType, CTType, RScore, MutIndex, TENDA, TENDG, TENDC, TENDT
			if ReportType == 'Heat Map':

				MapThisSeq = VGenesSeq.MutMap(Seq)
				MapThisGSeq = VGenesSeq.MutMap(GSeq)

			elif ReportType == 'GC-Heat Map':
				MapThisSeq = VGenesSeq.GCMutMap(Seq)
				MapThisGSeq = VGenesSeq.GCMutMap(GSeq)


			if ReportType == 'AnalysisDB Heat Map':

				NewMap.clear()
				n = 0
				w = 0
				# blank = ('.','.','.','.','.','.','.','.','.')
				for w in range(0, len(IMGTSeq) - 1):
					nuc = IMGTSeq[w]
					blank = (nuc, '.', '.', '.', '.', '.', '.', '.', '.')
					if n < len(MapThisSeq):
						if MapThisSeq[n][0] == nuc:
							ToAppend = MapThisSeq[n]
							NewMap.append(ToAppend)
							n += 1
						elif IMGTSeq[w] == '.':
							NewMap.append(blank)
					else:

						NewMap.append(blank)

				MapThisSeq.clear
				MapThisGSeq.clear
				MapThisSeq = NewMap
				MapThisGSeq = NewMap

			# (Nuc, MutType, CTType, RScore)

			if ReportType == 'Heat Map'  or ReportType == 'GC-Heat Map' or ReportType == 'AnalysisDB Heat Map':
				#     scores = (Nuc, MutType, CTType, RScore, MutIndex, TENDA, TENDG, TENDC, TENDT)
				#     forfields  = ['SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To','FR2To', 'CDR2To','FR3To', 'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend']
				F1SeqSet.clear()
				F1GermSet.clear()
				C1SeqSet.clear()
				C1GermSet.clear()
				F2SeqSet.clear()
				F2GermSet.clear()
				C2SeqSet.clear()
				C2GermSet.clear()
				F3SeqSet.clear()
				F3GermSet.clear()
				C3SeqSet.clear()
				C3GermSet.clear()
				F4SeqSet.clear()
				F4GermSet.clear()

				Vbeg = int(sequence[13])
				F1SeqSet.append(SeqName)
				F1GermSet.append(('Germline ' + SeqName))
				for i in range(1, Vbeg):  # to get all similar length and regions
					F1SeqSet.append(' ')
					F1GermSet.append(' ')

				for i in range(0, int(sequence[3])):  # i is base#
					F1SeqSet.append(MapThisSeq[i][4])
					F1GermSet.append(MapThisGSeq[i][4])
				FW1.append(tuple(F1SeqSet))
				GFW1.append(tuple(F1GermSet))

				for i in range(int(sequence[3]) + 1, int(sequence[4])):  # CDR1
					C1SeqSet.append(MapThisSeq[i][4])
					C1GermSet.append(MapThisGSeq[i][4])

				CDR1.append(tuple(C1SeqSet))
				GCDR1.append(tuple(C1GermSet))

				for i in range(int(sequence[4]) + 1, int(sequence[5])):  # FW2
					F2SeqSet.append(MapThisSeq[i][4])
					F2GermSet.append(MapThisGSeq[i][4])

				FW2.append(tuple(F2SeqSet))
				GFW2.append(tuple(F2GermSet))

				for i in range(int(sequence[5]) + 1, int(sequence[6])):  # CDR2
					C2SeqSet.append(MapThisSeq[i][4])
					C2GermSet.append(MapThisGSeq[i][4])

				CDR2.append(tuple(C2SeqSet))
				GCDR2.append(tuple(C2GermSet))

				GoTo = int(sequence[7])
				LenMap = len(MapThisSeq)
				if GoTo > LenMap:
					GoTo = LenMap - 1
				for i in range(int(sequence[6]) + 1, GoTo):  # FW3
					F3SeqSet.append(MapThisSeq[i][4])
					F3GermSet.append(MapThisGSeq[i][4])

				FW3.append(tuple(F3SeqSet))
				GFW3.append(tuple(F3GermSet))

				if (int(sequence[7]) + 1) < len(MapThisSeq):
					if int(sequence[8]) > len(MapThisSeq): GoTo = len(MapThisSeq) - 1
					for i in range(int(sequence[7]) + 1, GoTo):  # CDR3
						C3SeqSet.append(MapThisSeq[i][4])
						C3GermSet.append(MapThisGSeq[i][4])

				CDR3.append(tuple(C3SeqSet))
				GCDR3.append(tuple(C3GermSet))

				if len(MapThisSeq) > len(MapThisGSeq):
					F4end = len(MapThisGSeq)
				else:
					F4end = len(MapThisSeq)

				if (int(sequence[8]) + 1) < len(MapThisSeq):
					# if F4end> len(MapThisSeq): GoTo = len(MapThisSeq)-1
					for i in range(int(sequence[8]) + 1, F4end):  # FW4
						F4SeqSet.append(MapThisSeq[i][4])
						F4GermSet.append(MapThisGSeq[i][4])

				FW4.append(tuple(F4SeqSet))
				GFW4.append(tuple(F4GermSet))

			if ReportType == 'Mutation Frequencies':
				i = 0
				j = 1

				NumScored += 1
				if len(MapThisSeq) > len(MapThisGSeq):
					SeqLen = len(MapThisSeq) - 3
				else:
					SeqLen = len(MapThisGSeq) - 3

				for base in range(0, SeqLen):

					MutType = MapThisGSeq[base][1]

					# RScore = ''
					Nuc = MapThisSeq[base][0]
					GNuc = MapThisGSeq[base][0]
					TotType = GNuc + MutType
					Type = GNuc + Nuc + MutType

					if j == 1:
						Gcodon = MapThisGSeq[base][0] + MapThisGSeq[base + 1][0] + MapThisGSeq[base + 2][0]
						codon = Nuc + Gcodon[1] + Gcodon[2]
						GAA = VGenesSeq.Translator(Gcodon, 0)
						AA = VGenesSeq.Translator(codon, 0)
					elif j == 2:
						codon = Gcodon[0] + Nuc + Gcodon[2]
						AA = VGenesSeq.Translator(codon, 0)
					elif j == 3:
						codon = Gcodon[0] + Gcodon[1] + Nuc
						AA = VGenesSeq.Translator(codon, 0)
					# todo need total counts of each type of mutation to get frequency that each is mutated

					if Nuc != GNuc:
						if AA[0] == GAA[0]:
							if TotType in SMutations:
								Smuts += 1
								SMutations[TotType] += 1
							if Type in SMutations:
								SMutations[Type] += 1
						else:
							if TotType in RMutations:
								Rmuts += 1
								RMutations[TotType] += 1
							if Type in RMutations:
								RMutations[Type] += 1

					if AA[0] == GAA[0]:
						# Smuts +=1
						if TotType in SPositions:
							SPositions[TotType] += 1
							# SPositions[Type]+=1
					else:
						# Rmuts +=1
						if TotType in RPositions:
							RPositions[TotType] += 1
							# RPositions[Type]+=1

					if j < 3:
						j += 1
					else:
						j = 1

			if ReportType == 'Standard':
				for i in range(0, int(sequence[3])):
					if i == 0:

						Line = SeqName + ',' + MapThisSeq[i][0] + ',' + 'FW1' + ',' + MapThisSeq[i][1] + ',' + \
							   MapThisSeq[i][2] + ',' + \
							   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(
							MapThisSeq[i][5]) + ',' \
							   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
							   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
							   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(
							MapThisGSeq[i][5]) \
							   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
							MapThisGSeq[i][8]) + '\n'
					else:
						Line = ',' + MapThisSeq[i][0] + ',' + 'FW1' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
							2] + ',' + \
							   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(
							MapThisSeq[i][5]) + ',' \
							   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
							   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
							   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(
							MapThisGSeq[i][5]) \
							   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
							MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[3]) + 1, int(sequence[4])):
					Line = ',' + MapThisSeq[i][0] + ',' + 'CDR1' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[4]) + 1, int(sequence[5])):
					Line = ',' + MapThisSeq[i][0] + ',' + 'FW2' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[5]) + 1, int(sequence[6])):
					Line = ',' + MapThisSeq[i][0] + ',' + 'CDR2' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[6]) + 1, int(sequence[7])):
					Line = ',' + MapThisSeq[i][0] + ',' + 'FW3' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[7]) + 1, int(sequence[8])):
					Line = ',' + MapThisSeq[i][0] + ',' + 'CDR3' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					Document += Line
				for i in range(int(sequence[8]) + 1, len(MapThisSeq)):
					# try:
					Line = ',' + MapThisSeq[i][0] + ',' + 'FW4' + ',' + MapThisSeq[i][1] + ',' + MapThisSeq[i][
						2] + ',' + \
						   str(MapThisSeq[i][3]) + ',' + str(MapThisSeq[i][4]) + ',' + str(MapThisSeq[i][5]) + ',' \
						   + str(MapThisSeq[i][6]) + ',' + str(MapThisSeq[i][7]) + ',' + str(MapThisSeq[i][8]) + \
						   ',,' + MapThisGSeq[i][0] + ',' + MapThisGSeq[i][1] + ',' + MapThisGSeq[i][2] + ',' \
						   + str(MapThisGSeq[i][3]) + ',' + str(MapThisGSeq[i][4]) + ',' + str(MapThisGSeq[i][5]) \
						   + ',' + str(MapThisGSeq[i][6]) + ',' + str(MapThisGSeq[i][7]) + ',' + str(
						MapThisGSeq[i][8]) + '\n'
					if i == len(MapThisSeq):
						Line += '\n'
					# except:
					#     print('stop')

					Document += Line




			elif ReportType == 'Hotspots' or ReportType == 'HS-Summary':

				Header.append(SeqName)
				F1SeqSet.clear()
				F1GermSet.clear()

				for i in range(0, int(sequence[3])):
					F1SeqSet.append(MapThisSeq[i][1])
					F1GermSet.append(MapThisGSeq[i][1])

				if ReportType == 'HS-Summary':
					FW1counts = Counter(F1SeqSet)
					FW1Gcounts = Counter(F1GermSet)

					FW1keys = FW1counts.keys()
					if 'WRCY' not in FW1keys:
						FW1counts.update({'WRCY': 0})
					if 'WRC' not in FW1keys:
						FW1counts.update({'WRC': 0})
					if 'WA' not in FW1keys:
						FW1counts.update({'WA': 0})
					if 'TW' not in FW1keys:
						FW1counts.update({'TW': 0})
					if 'AID-CS' not in FW1keys:
						FW1counts.update({'AID-CS': 0})
					if 'Neutral' not in FW1keys:
						FW1counts.update({'Neutral': 0})

					FW1Gkeys = FW1Gcounts.keys()
					if 'WRCY' not in FW1Gkeys:
						FW1Gcounts.update({'WRCY': 0})
					if 'WRC' not in FW1Gkeys:
						FW1Gcounts.update({'WRC': 0})
					if 'WA' not in FW1Gkeys:
						FW1Gcounts.update({'WA': 0})
					if 'TW' not in FW1Gkeys:
						FW1Gcounts.update({'TW': 0})
					if 'AID-CS' not in FW1Gkeys:
						FW1Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in FW1Gkeys:
						FW1Gcounts.update({'Neutral': 0})

				else:
					if FW1len < len(F1SeqSet):
						FW1len = len(F1SeqSet)

					item = (tuple(F1SeqSet), tuple(F1GermSet))

					FW1.append(item)

				C1SeqSet.clear()
				C1GermSet.clear()
				for i in range(int(sequence[3]) + 1, int(sequence[4])):
					C1SeqSet.append(MapThisSeq[i][1])
					C1GermSet.append(MapThisGSeq[i][1])
				if ReportType == 'HS-Summary':
					CDR1counts = Counter(C1SeqSet)
					CDR1Gcounts = Counter(C1GermSet)

					CDR1keys = CDR1counts.keys()
					if 'WRCY' not in CDR1keys:
						CDR1counts.update({'WRCY': 0})
					if 'WRC' not in CDR1keys:
						CDR1counts.update({'WRC': 0})
					if 'WA' not in CDR1keys:
						CDR1counts.update({'WA': 0})
					if 'TW' not in CDR1keys:
						CDR1counts.update({'TW': 0})
					if 'AID-CS' not in CDR1keys:
						CDR1counts.update({'AID-CS': 0})
					if 'Neutral' not in CDR1keys:
						CDR1counts.update({'Neutral': 0})

					CDR1Gkeys = CDR1Gcounts.keys()
					if 'WRCY' not in CDR1Gkeys:
						CDR1Gcounts.update({'WRCY': 0})
					if 'WRC' not in CDR1Gkeys:
						CDR1Gcounts.update({'WRC': 0})
					if 'WA' not in CDR1Gkeys:
						CDR1Gcounts.update({'WA': 0})
					if 'TW' not in CDR1Gkeys:
						CDR1Gcounts.update({'TW': 0})
					if 'AID-CS' not in CDR1Gkeys:
						CDR1Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in CDR1Gkeys:
						CDR1Gcounts.update({'Neutral': 0})

				else:
					if CDR1len < len(C1SeqSet):
						CDR1len = len(C1SeqSet)

					item = (tuple(C1SeqSet), tuple(C1GermSet))

					CDR1.append(item)

				F2SeqSet.clear()
				F2GermSet.clear()
				for i in range(int(sequence[4]) + 1, int(sequence[5])):
					F2SeqSet.append(MapThisSeq[i][1])
					F2GermSet.append(MapThisGSeq[i][1])
				if ReportType == 'HS-Summary':
					FW2counts = Counter(F2SeqSet)
					FW2Gcounts = Counter(F2GermSet)

					FW2keys = FW2counts.keys()
					if 'WRCY' not in FW2keys:
						FW2counts.update({'WRCY': 0})
					if 'WRC' not in FW2keys:
						FW2counts.update({'WRC': 0})
					if 'WA' not in FW2keys:
						FW2counts.update({'WA': 0})
					if 'TW' not in FW2keys:
						FW2counts.update({'TW': 0})
					if 'AID-CS' not in FW2keys:
						FW2counts.update({'AID-CS': 0})
					if 'Neutral' not in FW2keys:
						FW2counts.update({'Neutral': 0})

					FW2Gkeys = FW2Gcounts.keys()
					if 'WRCY' not in FW2Gkeys:
						FW2Gcounts.update({'WRCY': 0})
					if 'WRC' not in FW2Gkeys:
						FW2Gcounts.update({'WRC': 0})
					if 'WA' not in FW2Gkeys:
						FW2Gcounts.update({'WA': 0})
					if 'TW' not in FW2Gkeys:
						FW2Gcounts.update({'TW': 0})
					if 'AID-CS' not in FW2Gkeys:
						FW2Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in FW2Gkeys:
						FW2Gcounts.update({'Neutral': 0})

				else:
					if FW2len < len(F2SeqSet):
						FW2len = len(F2SeqSet)

					item = (tuple(F2SeqSet), tuple(F2GermSet))

					FW2.append(item)

				C2SeqSet.clear()
				C2GermSet.clear()
				for i in range(int(sequence[5]) + 1, int(sequence[6])):
					C2SeqSet.append(MapThisSeq[i][1])
					C2GermSet.append(MapThisGSeq[i][1])
				if ReportType == 'HS-Summary':
					CDR2counts = Counter(C2SeqSet)
					CDR2Gcounts = Counter(C2GermSet)

					CDR2keys = CDR2counts.keys()
					if 'WRCY' not in CDR2keys:
						CDR2counts.update({'WRCY': 0})
					if 'WRC' not in CDR2keys:
						CDR2counts.update({'WRC': 0})
					if 'WA' not in CDR2keys:
						CDR2counts.update({'WA': 0})
					if 'TW' not in CDR2keys:
						CDR2counts.update({'TW': 0})
					if 'AID-CS' not in CDR2keys:
						CDR2counts.update({'AID-CS': 0})
					if 'Neutral' not in CDR2keys:
						CDR2counts.update({'Neutral': 0})

					CDR2Gkeys = CDR2Gcounts.keys()
					if 'WRCY' not in CDR2Gkeys:
						CDR2Gcounts.update({'WRCY': 0})
					if 'WRC' not in CDR2Gkeys:
						CDR2Gcounts.update({'WRC': 0})
					if 'WA' not in CDR2Gkeys:
						CDR2Gcounts.update({'WA': 0})
					if 'TW' not in CDR2Gkeys:
						CDR2Gcounts.update({'TW': 0})
					if 'AID-CS' not in CDR2Gkeys:
						CDR2Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in CDR2Gkeys:
						CDR2Gcounts.update({'Neutral': 0})

				else:
					if CDR2len < len(C2SeqSet):
						CDR2len = len(C2SeqSet)

					item = (tuple(C2SeqSet), tuple(C2GermSet))

					CDR2.append(item)

				F3SeqSet.clear()
				F3GermSet.clear()
				for i in range(int(sequence[6]) + 1, int(sequence[7])):
					F3SeqSet.append(MapThisSeq[i][1])
					F3GermSet.append(MapThisGSeq[i][1])
				if ReportType == 'HS-Summary':
					FW3counts = Counter(F3SeqSet)
					FW3Gcounts = Counter(F3GermSet)

					FW3keys = FW3counts.keys()
					if 'WRCY' not in FW3keys:
						FW3counts.update({'WRCY': 0})
					if 'WRC' not in FW3keys:
						FW3counts.update({'WRC': 0})
					if 'WA' not in FW3keys:
						FW3counts.update({'WA': 0})
					if 'TW' not in FW3keys:
						FW3counts.update({'TW': 0})
					if 'AID-CS' not in FW3keys:
						FW3counts.update({'AID-CS': 0})
					if 'Neutral' not in FW3keys:
						FW3counts.update({'Neutral': 0})

					FW3Gkeys = FW3Gcounts.keys()
					if 'WRCY' not in FW3Gkeys:
						FW3Gcounts.update({'WRCY': 0})
					if 'WRC' not in FW3Gkeys:
						FW3Gcounts.update({'WRC': 0})
					if 'WA' not in FW3Gkeys:
						FW3Gcounts.update({'WA': 0})
					if 'TW' not in FW3Gkeys:
						FW3Gcounts.update({'TW': 0})
					if 'AID-CS' not in FW3Gkeys:
						FW3Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in FW3Gkeys:
						FW3Gcounts.update({'Neutral': 0})

				else:
					if FW3len < len(F3SeqSet):
						FW3len = len(F3SeqSet)

					item = (tuple(F3SeqSet), tuple(F3GermSet))

					FW3.append(item)

				C3SeqSet.clear()
				C3GermSet.clear()
				for i in range(int(sequence[7]) + 1, int(sequence[8])):
					C3SeqSet.append(MapThisSeq[i][1])
					C3GermSet.append(MapThisGSeq[i][1])
				if ReportType == 'HS-Summary':
					CDR3counts = Counter(C3SeqSet)
					CDR3Gcounts = Counter(C3GermSet)

					CDR3keys = CDR3counts.keys()
					if 'WRCY' not in CDR3keys:
						CDR3counts.update({'WRCY': 0})
					if 'WRC' not in CDR3keys:
						CDR3counts.update({'WRC': 0})
					if 'WA' not in CDR3keys:
						CDR3counts.update({'WA': 0})
					if 'TW' not in CDR3keys:
						CDR3counts.update({'TW': 0})
					if 'AID-CS' not in CDR3keys:
						CDR3counts.update({'AID-CS': 0})
					if 'Neutral' not in CDR3keys:
						CDR3counts.update({'Neutral': 0})

					CDR3Gkeys = CDR3Gcounts.keys()
					if 'WRCY' not in CDR3Gkeys:
						CDR3Gcounts.update({'WRCY': 0})
					if 'WRC' not in CDR3Gkeys:
						CDR3Gcounts.update({'WRC': 0})
					if 'WA' not in CDR3Gkeys:
						CDR3Gcounts.update({'WA': 0})
					if 'TW' not in CDR3Gkeys:
						CDR3Gcounts.update({'TW': 0})
					if 'AID-CS' not in CDR3Gkeys:
						CDR3Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in CDR3Gkeys:
						CDR3Gcounts.update({'Neutral': 0})


				else:
					if CDR3len < len(C3SeqSet):
						CDR3len = len(C3SeqSet)

					item = (tuple(C3SeqSet), tuple(C3GermSet))

					CDR3.append(item)

				F4SeqSet.clear()
				F4GermSet.clear()
				for i in range(int(sequence[8]) + 1, len(MapThisGSeq)):
					# try:
					F4SeqSet.append(MapThisSeq[i][1])
					F4GermSet.append(MapThisGSeq[i][1])
					# except:
					#     print('stop')
				if ReportType == 'HS-Summary':
					FW4counts = Counter(F4SeqSet)
					FW4Gcounts = Counter(F4GermSet)

					FW4keys = FW4counts.keys()
					if 'WRCY' not in FW4keys:
						FW4counts.update({'WRCY': 0})
					if 'WRC' not in FW4keys:
						FW4counts.update({'WRC': 0})
					if 'WA' not in FW4keys:
						FW4counts.update({'WA': 0})
					if 'TW' not in FW4keys:
						FW4counts.update({'TW': 0})
					if 'AID-CS' not in FW4keys:
						FW4counts.update({'AID-CS': 0})
					if 'Neutral' not in FW4keys:
						FW4counts.update({'Neutral': 0})

					FW4Gkeys = FW4Gcounts.keys()
					if 'WRCY' not in FW4Gkeys:
						FW4Gcounts.update({'WRCY': 0})
					if 'WRC' not in FW4Gkeys:
						FW4Gcounts.update({'WRC': 0})
					if 'WA' not in FW4Gkeys:
						FW4Gcounts.update({'WA': 0})
					if 'TW' not in FW4Gkeys:
						FW4Gcounts.update({'TW': 0})
					if 'AID-CS' not in FW4Gkeys:
						FW4Gcounts.update({'AID-CS': 0})
					if 'Neutral' not in FW4Gkeys:
						FW4Gcounts.update({'Neutral': 0})

				else:
					if FW4len < len(F4SeqSet):
						FW4len = len(F4SeqSet)

					item = (tuple(F4SeqSet), tuple(F4GermSet))

					FW4.append(item)

			if ReportType == 'HS-Summary' or ReportType == 'Combined':
				# for key,value in FW1counts:
				FWWRCY = FW1counts.get('WRCY') + FW2counts.get('WRCY') + FW3counts.get('WRCY') + FW4counts.get(
					'WRCY')
				FWWRC = FW1counts.get('WRC') + FW2counts.get('WRC') + FW3counts.get('WRC') + FW4counts.get('WRC')
				FWWA = FW1counts.get('WA') + FW2counts.get('WA') + FW3counts.get('WA') + FW4counts.get('WA')
				FWTW = FW1counts.get('TW') + FW2counts.get('TW') + FW3counts.get('TW') + FW4counts.get('TW')
				FWAIDCS = FW1counts.get('AID-CS') + FW2counts.get('AID-CS') + FW3counts.get(
					'AID-CS') + FW4counts.get('AID-CS')
				FWNeut = FW1counts.get('Neutral') + FW2counts.get('Neutral') + FW3counts.get(
					'Neutral') + FW4counts.get('Neutral')

				CDRWRCY = CDR1counts.get('WRCY') + CDR2counts.get('WRCY') + CDR3counts.get('WRCY')
				CDRWRC = CDR1counts.get('WRC') + CDR2counts.get('WRC') + CDR3counts.get('WRC')
				CDRWA = CDR1counts.get('WA') + CDR2counts.get('WA') + CDR3counts.get('WA')
				CDRTW = CDR1counts.get('TW') + CDR2counts.get('TW') + CDR3counts.get('TW')
				CDRAIDCS = CDR1counts.get('AID-CS') + CDR2counts.get('AID-CS') + CDR3counts.get('AID-CS')
				CDRNeut = CDR1counts.get('Neutral') + CDR2counts.get('Neutral') + CDR3counts.get('Neutral')

				TotWRCY = FWWRCY + CDRWRCY
				TotWRC = FWWRC + CDRWRC
				TotWA = FWWA + CDRWA
				TotTW = FWTW + CDRTW
				TotAIDCS = FWAIDCS + CDRAIDCS
				TotNeut = FWNeut + CDRNeut

				GFWWRCY = FW1Gcounts.get('WRCY') + FW2Gcounts.get('WRCY') + FW3Gcounts.get('WRCY') + FW4Gcounts.get(
					'WRCY')
				GFWWRC = FW1Gcounts.get('WRC') + FW2Gcounts.get('WRC') + FW3Gcounts.get('WRC') + FW4Gcounts.get(
					'WRC')
				GFWWA = FW1Gcounts.get('WA') + FW2Gcounts.get('WA') + FW3Gcounts.get('WA') + FW4Gcounts.get('WA')
				GFWTW = FW1Gcounts.get('TW') + FW2Gcounts.get('TW') + FW3Gcounts.get('TW') + FW4Gcounts.get('TW')
				GFWAIDCS = FW1Gcounts.get('AID-CS') + FW2Gcounts.get('AID-CS') + FW3Gcounts.get(
					'AID-CS') + FW4Gcounts.get('AID-CS')
				GFWNeut = FW1Gcounts.get('Neutral') + FW2Gcounts.get('Neutral') + FW3Gcounts.get(
					'Neutral') + FW4Gcounts.get('Neutral')

				GCDRWRCY = CDR1Gcounts.get('WRCY') + CDR2Gcounts.get('WRCY') + CDR3Gcounts.get('WRCY')
				GCDRWRC = CDR1Gcounts.get('WRC') + CDR2Gcounts.get('WRC') + CDR3Gcounts.get('WRC')
				GCDRWA = CDR1Gcounts.get('WA') + CDR2Gcounts.get('WA') + CDR3Gcounts.get('WA')
				GCDRTW = CDR1Gcounts.get('TW') + CDR2Gcounts.get('TW') + CDR3Gcounts.get('TW')
				GCDRAIDCS = CDR1Gcounts.get('AID-CS') + CDR2Gcounts.get('AID-CS') + CDR3Gcounts.get('AID-CS')
				GCDRNeut = CDR1Gcounts.get('Neutral') + CDR2Gcounts.get('Neutral') + CDR3Gcounts.get('Neutral')
				GTotWRCY = GFWWRCY + GCDRWRCY
				GTotWRC = GFWWRC + GCDRWRC
				GTotWA = GFWWA + GCDRWA
				GTotTW = GFWTW + GCDRTW
				GTotAIDCS = GFWAIDCS + GCDRAIDCS
				GTotNeut = GFWNeut + GCDRNeut

				SumTup = (SeqName, str(mutcount), str(TotWRCY), str(TotWRC), str(TotWA), str(TotTW), str(TotAIDCS),
						  str(TotNeut)
						  , str(FWWRCY), str(FWWRC), str(FWWA), str(FWTW), str(FWAIDCS), str(FWNeut)
						  , str(CDRWRCY), str(CDRWRC), str(CDRWA), str(CDRTW), str(CDRAIDCS), str(CDRNeut)
						  , FW1counts.get('WRCY'), FW1counts.get('WRC'), FW1counts.get('WA'), FW1counts.get('TW'),
						  FW1counts.get('AID-CS'), FW1counts.get('Neutral')
						  , CDR1counts.get('WRCY'), CDR1counts.get('WRC'), CDR1counts.get('WA'),
						  CDR1counts.get('TW'), CDR1counts.get('AID-CS'), CDR1counts.get('Neutral')
						  , FW2counts.get('WRCY'), FW2counts.get('WRC'), FW2counts.get('WA'), FW2counts.get('TW'),
						  FW2counts.get('AID-CS'), FW2counts.get('Neutral')
						  , CDR2counts.get('WRCY'), CDR2counts.get('WRC'), CDR2counts.get('WA'),
						  CDR2counts.get('TW'), CDR2counts.get('AID-CS'), CDR2counts.get('Neutral')
						  , FW3counts.get('WRCY'), FW3counts.get('WRC'), FW3counts.get('WA'), FW3counts.get('TW'),
						  FW3counts.get('AID-CS'), FW3counts.get('Neutral')
						  , CDR3counts.get('WRCY'), CDR3counts.get('WRC'), CDR3counts.get('WA'),
						  CDR3counts.get('TW'), CDR3counts.get('AID-CS'), CDR3counts.get('Neutral')
						  , FW4counts.get('WRCY'), FW4counts.get('WRC'), FW4counts.get('WA'), FW4counts.get('TW'),
						  FW4counts.get('AID-CS'), FW4counts.get('Neutral')
						  , str(GTotWRCY), str(GTotWRC), str(GTotWA), str(GTotTW), str(GTotAIDCS), str(GTotNeut)
						  , str(GFWWRCY), str(GFWWRC), str(GFWWA), str(GFWTW), str(GFWAIDCS), str(GFWNeut)
						  , str(GCDRWRCY), str(GCDRWRC), str(GCDRWA), str(GCDRTW), str(GCDRAIDCS), str(GCDRNeut)
						  , FW1Gcounts.get('WRCY'), FW1Gcounts.get('WRC'), FW1Gcounts.get('WA'),
						  FW1Gcounts.get('TW'), FW1Gcounts.get('AID-CS'), FW1Gcounts.get('Neutral')
						  , CDR1Gcounts.get('WRCY'), CDR1Gcounts.get('WRC'), CDR1Gcounts.get('WA'),
						  CDR1Gcounts.get('TW'), CDR1Gcounts.get('AID-CS'), CDR1Gcounts.get('Neutral')
						  , FW2Gcounts.get('WRCY'), FW2Gcounts.get('WRC'), FW2Gcounts.get('WA'),
						  FW2Gcounts.get('TW'), FW2Gcounts.get('AID-CS'), FW2Gcounts.get('Neutral')
						  , CDR2Gcounts.get('WRCY'), CDR2Gcounts.get('WRC'), CDR2Gcounts.get('WA'),
						  CDR2Gcounts.get('TW'), CDR2Gcounts.get('AID-CS'), CDR2Gcounts.get('Neutral')
						  , FW3Gcounts.get('WRCY'), FW3Gcounts.get('WRC'), FW3Gcounts.get('WA'),
						  FW3Gcounts.get('TW'), FW3Gcounts.get('AID-CS'), FW3Gcounts.get('Neutral')
						  , CDR3Gcounts.get('WRCY'), CDR3Gcounts.get('WRC'), CDR3Gcounts.get('WA'),
						  CDR3Gcounts.get('TW'), CDR3Gcounts.get('AID-CS'), CDR3Gcounts.get('Neutral')
						  , FW4Gcounts.get('WRCY'), FW4Gcounts.get('WRC'), FW4Gcounts.get('WA'),
						  FW4Gcounts.get('TW'), FW4Gcounts.get('AID-CS'), FW4Gcounts.get('Neutral'))

				Summary.append(SumTup)

		except:
			error = sys.exc_info()[0]
			errrep += SeqName + ' caused: ' + str(error)

	if ReportType == 'Heat Map'  or ReportType == 'GC-Heat Map' or ReportType == 'AnalysisDB Heat Map':
		if ReportType == 'AnalysisDB Heat Map':
			item = 'No germline'
		else:
			items = ('Stacked', 'Paired', 'No germline')
			title = 'Choose germline analysis placement:'
			item = setItem(self, items, title)
		if item == "Cancel":
			return

		Document = ''
		FW1len = FW1len + 5
		for j in range(0, FW1len):  # each row or base
			line = 'FW1,'
			for i in range(0, len(FW1)):  # i is each column or each sequence
				try:
					line += (str(FW1[i][j]) + ',')
					if item == 'Paired':
						# if i == 0:
						# Line+= 'Germline '
						line += (str(GFW1[i][j]) + ',')
				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'
			Document += (line + '\n')
		for j in range(0, CDR1len):  # each row or base
			line = 'CDR1,'
			for i in range(0, len(CDR1)):  # i is each column or each sequence
				try:
					line += (str(CDR1[i][j]) + ',')
					if item == 'Paired':
						line += (str(GCDR1[i][j]) + ',')
				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')
		FW2len = FW2len + 5
		for j in range(0, FW2len):  # each row or base
			line = 'FW2,'
			for i in range(0, len(FW2)):  # i is each column or each sequence
				try:
					line += (str(FW2[i][j]) + ',')
					if item == 'Paired':
						line += (str(GFW2[i][j]) + ',')
				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')
		for j in range(0, CDR2len):  # each row or base
			line = 'CDR2,'
			for i in range(0, len(CDR2)):  # i is each column or each sequence
				try:
					line += (str(CDR2[i][j]) + ',')
					if item == 'Paired':
						line += (str(GCDR2[i][j]) + ',')
				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')

		for j in range(0, FW3len):  # each row or base
			line = 'FW3,'
			for i in range(0, len(FW3)):  # i is each column or each sequence
				try:
					line += (str(FW3[i][j]) + ',')
					if item == 'Paired':
						line += (str(GFW3[i][j]) + ',')
				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')
		for j in range(0, CDR3len):  # each row or base
			line = 'CDR3,'
			for i in range(0, len(CDR3)):  # i is each column or each sequence
				try:
					line += (str(CDR3[i][j]) + ',')
					if item == 'Paired':
						line += (str(GCDR3[i][j]) + ',')

				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')

		for j in range(0, FW4len):  # each row or base
			line = 'FW4,'
			for i in range(0, len(FW4)):  # i is each column or each sequence
				try:
					line += (str(FW4[i][j]) + ',')
					if item == 'Paired':
						line += (str(GFW4[i][j]) + ',')

				except:
					line += ' ,'
					if item == 'Paired':
						line += ' ,'

			Document += (line + '\n')

		if item == 'Stacked':
			Document += '\nGermline:\n'
			for j in range(0, FW1len):  # each row or base
				line = 'FW1,'
				for i in range(0, len(GFW1)):  # i is each column or each sequence
					try:
						line += (str(GFW1[i][j]) + ',')

					except:
						line += ' ,'

				Document += (line + '\n')
			for j in range(0, CDR1len):  # each row or base
				line = 'CDR1,'
				for i in range(0, len(GCDR1)):  # i is each column or each sequence
					try:
						line += (str(GCDR1[i][j]) + ',')

					except:
						line += ' ,'

				Document += (line + '\n')

			for j in range(0, FW2len):  # each row or base
				line = 'FW2,'
				for i in range(0, len(GFW2)):  # i is each column or each sequence
					try:
						line += (str(GFW2[i][j]) + ',')

					except:
						line += ' ,'

				Document += (line + '\n')
			for j in range(0, CDR2len):  # each row or base
				line = 'CDR2,'
				for i in range(0, len(GCDR2)):  # i is each column or each sequence
					try:
						line += (str(GCDR2[i][j]) + ',')

					except:
						line += ' ,'

				Document += (line + '\n')

			for j in range(0, FW3len):  # each row or base
				line = 'FW3,'
				for i in range(0, len(GFW3)):  # i is each column or each sequence
					try:
						line += (str(GFW3[i][j]) + ',')

					except:
						line += ' ,'

				Document += (line + '\n')
			for j in range(0, CDR3len):  # each row or base
				line = 'CDR3,'
				for i in range(0, len(GCDR3)):  # i is each column or each sequence
					try:
						line += (str(GCDR3[i][j]) + ',')


					except:
						line += ' ,'

				Document += (line + '\n')

			for j in range(0, FW4len):  # each row or base
				line = 'FW4,'
				for i in range(0, len(GFW4)):  # i is each column or each sequence
					try:
						line += (str(GFW4[i][j]) + ',')


					except:
						line += ' ,'

				Document += (line + '\n')

	if ReportType == 'Mutation Frequencies':
		# NumScored = 0
		# Smuts = 0
		# Rmuts = 0
		# SMutations = {'AWA':0, 'ANeutral':0, 'TTW':0, 'TNeutral':0, 'CWRCY':0
		#     , 'CWRC':0, 'CAID-CS':0, 'CNeutral':0, 'GWRCY':0, 'GWRC':0, 'GAID-CS':0, 'GNeutral':0
		#     , 'ATWA':0, 'ATNeutral':0, 'AGWA':0, 'AGNeutral':0, 'ACWA':0, 'ACNeutral':0
		#     , 'TAWA':0, 'TANeutral':0, 'TGWA':0, 'TGNeutral':0, 'TCWA':0, 'TCNeutral':0
		#     , 'CAWRCY':0, 'CAWRC':0, 'CAAID-CS':0, 'CANeutral':0
		#     , 'CTWRCY':0, 'CTWRC':0, 'CTAID-CS':0, 'CTNeutral':0
		#     , 'CGWRCY':0, 'CGWRC':0, 'CGAID-CS':0, 'CGNeutral':0
		#     , 'GAWRCY':0, 'GAWRC':0, 'GAAID-CS':0, 'GANeutral':0
		#     , 'GTWRCY':0, 'GTWRC':0, 'GTAID-CS':0, 'GTNeutral':0
		#     , 'GCWRCY':0, 'GCWRC':0, 'GCAID-CS':0, 'GCNeutral':0}
		Document = 'Sequences analyzed: ,' + str(NumScored) + '\n'
		Document += 'Silent Mutations: ,' + str(Smuts) + '\n'
		Document += ',A hotspot, A neutral, AT hotspot, AT neutral, AG hotspot, AG neutral, AC hotspot, AC neutral' \
					', T hotspot, T neutral, TA hotspot, TA neutral, TG hotspot, TG neutral, TC hotspot, TC neutral' \
					', C WRCY, C WRC, C coldspot, C neutral, CA WRCY, CA WRC, CA coldspot, CA neutral, CT WRCY, CT WRC' \
					', CT coldspot, CT neutral, CG WRCY, CG WRC, CG coldspot, CG neutral' \
					', G WRCY, G WRC, G coldspot, G neutral, GC WRCY, GC WRC, GC coldspot, GC neutral, GA WRCY, GA WRC' \
					', GA coldspot, GA neutral, GT WRCY, GT WRC, GT coldspot, GT neutral\n'
		Document += (
		'Mutations:,' + str(SMutations['AWA']) + ',' + str(SMutations['ANeutral']) + ',' + str(SMutations['ATWA']) +
		',' + str(SMutations['ATNeutral']) + ',' + str(SMutations['AGWA']) + ',' + str(SMutations['AGNeutral'])
		+ ',' + str(SMutations['ACWA']) + ',' + str(SMutations['ACNeutral']) + ',' + str(
			SMutations['TTW']) + ',' + str(SMutations['TNeutral'])
		+ ',' + str(SMutations['TATW']) + ',' + str(SMutations['TANeutral']) + ',' + str(
			SMutations['TGTW']) + ',' + str(SMutations['TGNeutral'])
		+ ',' + str(SMutations['TCTW']) + ',' + str(SMutations['TCNeutral']) + ','
		+ str(SMutations['CWRCY']) + ',' + str(SMutations['CWRC']) + ',' + str(SMutations['CAID-CS']) + ',' + str(
			SMutations['CNeutral']) + ','
		+ str(SMutations['CAWRCY']) + ',' + str(SMutations['CAWRC']) + ',' + str(
			SMutations['CAAID-CS']) + ',' + str(SMutations['CANeutral']) + ','
		+ str(SMutations['CTWRCY']) + ',' + str(SMutations['CTWRC']) + ',' + str(
			SMutations['CTAID-CS']) + ',' + str(SMutations['CTNeutral']) + ','
		+ str(SMutations['CGWRCY']) + ',' + str(SMutations['CGWRC']) + ',' + str(
			SMutations['CGAID-CS']) + ',' + str(SMutations['CGNeutral']) + ','
		+ str(SMutations['GWRCY']) + ',' + str(SMutations['GWRC']) + ',' + str(SMutations['GAID-CS']) + ',' + str(
			SMutations['GNeutral']) + ','
		+ str(SMutations['GCWRCY']) + ',' + str(SMutations['GCWRC']) + ',' + str(
			SMutations['GCAID-CS']) + ',' + str(SMutations['GCNeutral']) + ','
		+ str(SMutations['GAWRCY']) + ',' + str(SMutations['GAWRC']) + ',' + str(
			SMutations['GAAID-CS']) + ',' + str(SMutations['GANeutral']) + ','
		+ str(SMutations['GTWRCY']) + ',' + str(SMutations['GTWRC']) + ',' + str(
			SMutations['GTAID-CS']) + ',' + str(SMutations['GTNeutral']) + '\n')

		Document += (
		'Positions:,' + str(SPositions['AWA']) + ',' + str(SPositions['ANeutral']) + ',' + str(SPositions['AWA']) +
		',' + str(SPositions['ANeutral']) + ',' + str(SPositions['AWA']) + ',' + str(SPositions['ANeutral'])
		+ ',' + str(SPositions['AWA']) + ',' + str(SPositions['ANeutral']) + ',' + str(
			SPositions['TTW']) + ',' + str(SPositions['TNeutral'])
		+ ',' + str(SPositions['TTW']) + ',' + str(SPositions['TNeutral']) + ',' + str(
			SPositions['TTW']) + ',' + str(SPositions['TNeutral'])
		+ ',' + str(SPositions['TTW']) + ',' + str(SPositions['TNeutral']) + ','
		+ str(SPositions['CWRCY']) + ',' + str(SPositions['CWRC']) + ',' + str(SPositions['CAID-CS']) + ',' + str(
			SPositions['CNeutral']) + ','
		+ str(SPositions['CWRCY']) + ',' + str(SPositions['CWRC']) + ',' + str(SPositions['CAID-CS']) + ',' + str(
			SPositions['CNeutral']) + ','
		+ str(SPositions['CWRCY']) + ',' + str(SPositions['CWRC']) + ',' + str(SPositions['CAID-CS']) + ',' + str(
			SPositions['CNeutral']) + ','
		+ str(SPositions['CWRCY']) + ',' + str(SPositions['CWRC']) + ',' + str(SPositions['CAID-CS']) + ',' + str(
			SPositions['CNeutral']) + ','
		+ str(SPositions['GWRCY']) + ',' + str(SPositions['GWRC']) + ',' + str(SPositions['GAID-CS']) + ',' + str(
			SPositions['GNeutral']) + ','
		+ str(SPositions['GWRCY']) + ',' + str(SPositions['GWRC']) + ',' + str(SPositions['GAID-CS']) + ',' + str(
			SPositions['GNeutral']) + ','
		+ str(SPositions['GWRCY']) + ',' + str(SPositions['GWRC']) + ',' + str(SPositions['GAID-CS']) + ',' + str(
			SPositions['GNeutral']) + ','
		+ str(SPositions['GWRCY']) + ',' + str(SPositions['GWRC']) + ',' + str(SPositions['GAID-CS']) + ',' + str(
			SPositions['GNeutral']) + '\n\n')

		Document += 'Replacment Mutations: ,' + str(Rmuts) + '\n'
		Document += ' , A hotspot, A neutral, AT hotspot, AT neutral, AG hotspot, AG neutral, AC hotspot, AC neutral' \
					', T hotspot, T neutral, TA hotspot, TA neutral, TG hotspot, TG neutral, TC hotspot, TC neutral' \
					', C WRCY, C WRC, C coldspot, C neutral, CA WRCY, CA WRC, CA coldspot, CA neutral, CT WRCY, CT WRC' \
					', CT coldspot, CT neutral, CG WRCY, CG WRC, CG coldspot, CG neutral' \
					', G WRCY, G WRC, G coldspot, G neutral, GC WRCY, GC WRC, GC coldspot, GC neutral, GA WRCY, GA WRC' \
					', GA coldspot, GA neutral, GT WRCY, GT WRC, GT coldspot, GT neutral\n'
		Document += (
		'Mutations:,' + str(RMutations['AWA']) + ',' + str(RMutations['ANeutral']) + ',' + str(RMutations['ATWA']) +
		',' + str(RMutations['ATNeutral']) + ',' + str(RMutations['AGWA']) + ',' + str(RMutations['AGNeutral'])
		+ ',' + str(RMutations['ACWA']) + ',' + str(RMutations['ACNeutral']) + ',' + str(
			RMutations['TTW']) + ',' + str(RMutations['TNeutral'])
		+ ',' + str(RMutations['TATW']) + ',' + str(RMutations['TANeutral']) + ',' + str(
			RMutations['TGTW']) + ',' + str(RMutations['TGNeutral'])
		+ ',' + str(RMutations['TCTW']) + ',' + str(RMutations['TCNeutral']) + ','
		+ str(RMutations['CWRCY']) + ',' + str(RMutations['CWRC']) + ',' + str(RMutations['CAID-CS']) + ',' + str(
			RMutations['CNeutral']) + ','
		+ str(RMutations['CAWRCY']) + ',' + str(RMutations['CAWRC']) + ',' + str(
			RMutations['CAAID-CS']) + ',' + str(RMutations['CANeutral']) + ','
		+ str(RMutations['CTWRCY']) + ',' + str(RMutations['CTWRC']) + ',' + str(
			RMutations['CTAID-CS']) + ',' + str(RMutations['CTNeutral']) + ','
		+ str(RMutations['CGWRCY']) + ',' + str(RMutations['CGWRC']) + ',' + str(
			RMutations['CGAID-CS']) + ',' + str(RMutations['CGNeutral']) + ','
		+ str(RMutations['GWRCY']) + ',' + str(RMutations['GWRC']) + ',' + str(RMutations['GAID-CS']) + ',' + str(
			RMutations['GNeutral']) + ','
		+ str(RMutations['GCWRCY']) + ',' + str(RMutations['GCWRC']) + ',' + str(
			RMutations['GCAID-CS']) + ',' + str(RMutations['GCNeutral']) + ','
		+ str(RMutations['GAWRCY']) + ',' + str(RMutations['GAWRC']) + ',' + str(
			RMutations['GAAID-CS']) + ',' + str(RMutations['GANeutral']) + ','
		+ str(RMutations['GTWRCY']) + ',' + str(RMutations['GTWRC']) + ',' + str(
			RMutations['GTAID-CS']) + ',' + str(RMutations['GTNeutral']) + '\n')

		Document += (
		'Positions:,' + str(RPositions['AWA']) + ',' + str(RPositions['ANeutral']) + ',' + str(RPositions['AWA']) +
		',' + str(RPositions['ANeutral']) + ',' + str(RPositions['AWA']) + ',' + str(RPositions['ANeutral'])
		+ ',' + str(RPositions['AWA']) + ',' + str(RPositions['ANeutral']) + ',' + str(
			RPositions['TTW']) + ',' + str(RPositions['TNeutral'])
		+ ',' + str(RPositions['TTW']) + ',' + str(RPositions['TNeutral']) + ',' + str(
			RPositions['TTW']) + ',' + str(RPositions['TNeutral'])
		+ ',' + str(RPositions['TTW']) + ',' + str(RPositions['TNeutral']) + ','
		+ str(RPositions['CWRCY']) + ',' + str(RPositions['CWRC']) + ',' + str(RPositions['CAID-CS']) + ',' + str(
			RPositions['CNeutral']) + ','
		+ str(RPositions['CWRCY']) + ',' + str(RPositions['CWRC']) + ',' + str(RPositions['CAID-CS']) + ',' + str(
			RPositions['CNeutral']) + ','
		+ str(RPositions['CWRCY']) + ',' + str(RPositions['CWRC']) + ',' + str(RPositions['CAID-CS']) + ',' + str(
			RPositions['CNeutral']) + ','
		+ str(RPositions['CWRCY']) + ',' + str(RPositions['CWRC']) + ',' + str(RPositions['CAID-CS']) + ',' + str(
			RPositions['CNeutral']) + ','
		+ str(RPositions['GWRCY']) + ',' + str(RPositions['GWRC']) + ',' + str(RPositions['GAID-CS']) + ',' + str(
			RPositions['GNeutral']) + ','
		+ str(RPositions['GWRCY']) + ',' + str(RPositions['GWRC']) + ',' + str(RPositions['GAID-CS']) + ',' + str(
			RPositions['GNeutral']) + ','
		+ str(RPositions['GWRCY']) + ',' + str(RPositions['GWRC']) + ',' + str(RPositions['GAID-CS']) + ',' + str(
			RPositions['GNeutral']) + ','
		+ str(RPositions['GWRCY']) + ',' + str(RPositions['GWRC']) + ',' + str(RPositions['GAID-CS']) + ',' + str(
			RPositions['GNeutral']) + '\n')

	if ReportType == 'HS-Summary' or ReportType == 'Combined':
		for item in Summary:
			for i in range(0, len(item)):
				if i != 121:
					Document += (str(item[i]) + ',')
				else:
					Document += (str(item[i]) + '\n')

	if ReportType == 'Hotspots' or ReportType == 'Combined':
		# counts = {'WRCY':0, 'WRC':0, 'WA':0, 'TW':0, 'AID-CS':0, 'Neutral':0}

		line = 'FW1 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, FW1len + 2):
			line = ','
			j = 0
			for seq in FW1:
				try:
					entry = seq[0][i]
					# counts.append(entry)
				except:
					entry = ''
				if j == len(FW1) - 1:

					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW1 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, FW1len + 2):
			line = ','
			j = 0
			for seq in FW1:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(FW1) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR1 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, CDR1len + 2):
			line = ','
			j = 0
			for seq in CDR1:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(CDR1) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR1 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, CDR1len + 2):
			line = ','
			j = 0
			for seq in CDR1:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(CDR1) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW2 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, FW2len + 2):
			line = ','
			j = 0
			for seq in FW2:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(FW2) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW2 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, FW2len + 2):
			line = ','
			j = 0
			for seq in FW2:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(FW2) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR2 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, CDR2len + 2):
			line = ','
			j = 0
			for seq in CDR2:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(CDR2) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR2 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, CDR2len + 2):
			line = ','
			j = 0
			for seq in CDR2:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(CDR2) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW3 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, FW3len + 2):
			line = ','
			j = 0
			for seq in FW3:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(FW3) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW3 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, FW3len + 2):
			line = ','
			j = 0
			for seq in FW3:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(FW3) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR3 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, CDR3len + 2):
			line = ','
			j = 0
			for seq in CDR3:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(CDR3) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'CDR3 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, CDR3len + 2):
			line = ','
			j = 0
			for seq in CDR3:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(CDR3) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW4 Sequence:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0

		for i in range(0, FW4len + 2):
			line = ','
			j = 0
			for seq in FW4:
				try:
					entry = seq[0][i]
				except:
					entry = ''
				if j == len(FW4) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

		line = 'FW4 germline:,'
		NumSeqs = len(Header)
		i = 0
		for item in Header:
			line += item
			if i != NumSeqs - 1:
				line += ','
			else:
				line += '\n'
			i += 1
		Document += line
		i = 0
		for i in range(0, FW4len + 2):
			line = ','
			j = 0
			for seq in FW4:
				try:
					entry = seq[1][i]
				except:
					entry = ''
				if j == len(FW4) - 1:
					line += (entry + '\n')
				else:
					line += (entry + ',')
				j += 1
			Document += line

	# elif ReportType == 'HS-Summary':


	if errrep != '\nThe following errors occured during processing:\n':
		Style = 'standard'

		self.ShowVGenesTextEdit(errrep, Style)
	f = ''
	f = saveFile(self, 'CSV')
	if f != '' and f != None:
		with open(f, 'w') as currentfile:
			currentfile.write(Document)

	return 'done'