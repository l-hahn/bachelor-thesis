from random import choice as SeqChoice

print("Reading Pfam-File ... ",end="")
Data = open("/data/PFam/Pfam-A.full","r",encoding="latin-1")
print("Done!")
TrainDB = open("TrainDatabase.stockholm","w",encoding="latin-1")
TestDB = open("TestDatabase.stockholm","w",encoding="latin-1")

Pfid = ""
Name = ""
Descr = ""
Seqs = []
TestOffset = 0.25
MinLeng = 10
FamCtr = 0

print("\nSpliting Database ... ")
for Line in Data:
	if Line[0] == '#':
		Input = Line.split()
		Input = list(filter(None,Input))
		if len(Input) > 0:
			if Input[1] == "STOCKHOLM":
				if len(Seqs) >= MinLeng:
					print("\rWriting {}\t{}".format(Pfid, FamCtr),end="")
					TrainDB.write(Line)
					TrainDB.write("#=GF ID   {}\n#=GF AC   {}\n#=GF DE   {}\n".format(Name, Pfid, Descr))
					TestDB.write(Line)
					TestDB.write("#=GF ID   {}\n#=GF AC   {}\n#=GF DE   {}\n".format(Name, Pfid, Descr))

					Offset = int((TestOffset)*len(Seqs))
					SeqCtr = 0
					TestSeq = []
					while SeqCtr < Offset:
						SeqEntry = SeqChoice(Seqs)
						TestSeq.append(SeqEntry)
						Seqs.remove(SeqEntry)
						SeqCtr += 1
					for SeqLine in Seqs:
						TrainDB.write(SeqLine)
					for SeqLine in TestSeq:
						TestDB.write(SeqLine)
				Seqs = []
				FamCtr += 1
			if Input[0] == "#=GF":
				if Input[1] == "ID":
					Name = Input[2]
				elif Input[1] == "AC":
					Pfid = Input[2]
					print("\rReading {}\t{}".format(Pfid, FamCtr),end="")
				elif Input[1] == "DE":
					Descr = " ".join(Input[2:])
	else:
		Seqs.append(Line)
if len(Seq) > 0:
	print("\rWriting {}\t{}".format(Pfid, FamCtr),end="")
	TrainDB.write(Line)
	TrainDB.write("#=GF ID   {}\n#=GF AC   {}\n#=GF DE   {}\n".format(Name, Pfid, Descr))
	TestDB.write(Line)
	TestDB.write("#=GF ID   {}\n#=GF AC   {}\n#=GF DE   {}\n".format(Name, Pfid, Descr))

	Offset = int((TestOffset)*len(Seqs))
	SeqCtr = 0
	TestSeq = []
	while SeqCtr < Offset:
		SeqEntry = SeqChoice(Seqs)
		TestSeq.append(SeqEntry)
		Seqs.remove(SeqEntry)
		SeqCtr += 1
	for SeqLine in Seqs:
		TrainDB.write(SeqLine)
	for SeqLine in TestSeq:
		TestDB.write(SeqLine)
CleanLine = ' '*80;
print(CleanLine,end="")
print("\rDone!")
Data.close()
TrainDB.close()
TestDB.close()