"""
This module are used to parse pdb files, mainly for being imported of other scripts, containing:
	class pdb_atom_record:
		method __init__(self, string): string is a line of "ATOM" or "HETATM" from pdb file
		method update_recordName(self, new_recordName):
		method update_serial(self, new_serial=None, by_shift=0): serial is int
		method update_resName(self, new_resName)
		method update_chainID(self, new_chainID)
		method update_resSeq(self, new_resSeq=None, by_shift=0): resSeq is int
		method update_xyz(self, *new_xyz): x, y, z are floats

	class pdb_ter_record: inherited from pdb_atom_record
		method __init__(self, string):
			only set members string, resSeq and chainID
			must use methods update_serial and update_resSeq for setting serial and resSeq

	function file2rec(inpfile, rec_list): read file to a list of [pdb_atom_record]

	function rec2file(rec_list, outfile, reorder_serial = False): write a list of [pdb_atom_record] to a file
		optional reorder_serial is for reorder the serial for atoms starting from 1
version update:
20240624: accepted revision from Ash related to TER records
"""

class pdb_atom_record:#for ATOM or HETATM
	def __init__(self, string):
		self.string = string
		self.recordName = string[:6].strip()
		self.serial = int(string[6:11])
		self.name = string[12:16].strip()
		self.resName = string[17:20].strip()
		self.chainID = string[21]
		self.resSeq = int(string[22:26])
		self.x =float(string[30:38])
		self.y =float(string[38:46])
		self.z =float(string[46:54])

	def update_recordName(self, new_recordName):
		self.recordName = new_recordName
		self.string = ('%-6s' % self.recordName) + self.string[6:]

	def update_serial(self, new_serial=None, by_shift=0):
		if new_serial is None:
			new_serial = self.serial
		self.serial = new_serial + by_shift
		self.string = self.string[:6] + ('%5d' % self.serial) + self.string[11:]

	def update_resName(self, new_resName):
		self.resName = new_resName
		self.string = self.string[:17] + ('%3s' % self.resName) + self.string[20:]

	def update_chainID(self, new_chainID):
		self.chainID = new_chainID
		self.string = self.string[:21] + ('%1s' % self.chainID) + self.string[22:]

	def update_resSeq(self, new_resSeq=None, by_shift=0):
		if new_resSeq is None:
			new_resSeq = self.resSeq
		self.resSeq = new_resSeq + by_shift
		self.string = self.string[:22] + ('%4d' % self.resSeq) + self.string[26:]
	def update_xyz(self, *new_xyz):
		self.x, self.y, self.z = new_xyz
		self.string = self.string[:30]+ "%8.3f%8.3f%8.3f"%(self.x, self.y, self.z) + self.string[54:]

class pdb_ter_record(pdb_atom_record):#for TER, inherit from ATOM
	def __init__(self, string):
		self.string = string[:-1] + " "*( 27-len(string.rstrip()) ) + string[-1] # "\n" @-1
		# print(self.string)
		self.recordName = string[:6].strip()

		serial_str = string[6:11].strip()
		if serial_str.isdigit():
			self.serial = int(serial_str) # to prevent the error from that some file TER is not numbered
		else:
			self.serial = None
		self.resName = self.string[17:20].strip()
		self.chainID = self.string[21]
		resSeq_str = string[22:26].strip()
		if resSeq_str.isdigit():
			self.resSeq = int(resSeq_str)
		else:
			self.resSeq = None


def file2rec(inpfile, rec_list): # read from file to a list
	last_serial = 0
	last_resSeq = 0
	added_ter_serial = 0 # recording added serial for TER having no serial
	for eachline in inpfile:
		if eachline[:4] == "ATOM" or eachline[:6] == "HETATM":
			curr_atom = pdb_atom_record(eachline)
			curr_atom_serial = curr_atom.serial + added_ter_serial
			curr_atom.update_serial(curr_atom_serial)
			last_serial = curr_atom.serial
			last_resSeq = curr_atom.resSeq
		elif eachline[:3] == "TER":
			curr_atom = pdb_ter_record(eachline)
			if curr_atom.serial is None:
				curr_atom.update_serial(last_serial+1)
				added_ter_serial = added_ter_serial +1
			if curr_atom.resSeq is None:
				curr_atom.update_resSeq(last_resSeq)
		else:
			continue
		rec_list.append(curr_atom)

def rec2file(rec_list, outfile, reorder_serial = False): # write to file
	for x, write_atom in enumerate(rec_list):
		if reorder_serial:
			write_atom.update_serial(x+1)
		outfile.write(write_atom.string)

if __name__ == '__main__':
	test_string = 'ATOM     12  N3    A A   1       4.730  -2.263  -2.712  1.00  0.00           N '
	atom = pdb_atom_record(test_string)
	print('Original:\n', atom.string)
	atom.update_serial(1)
	atom.update_resName('G')
	atom.update_chainID('X')
	atom.update_resSeq(999)
	atom.update_resSeq(by_shift=-9)
	xyzlist=[0.0,1.1,-3.3]
	atom.update_xyz(*xyzlist)
	print('New:\n', atom.string)
	test_string2 = 'TER    1890      G   E  99'
	# test_string2 = 'TER              G   E  99'
	ter = pdb_ter_record(test_string2)
	print('Original:\n', ter.string)
	ter.update_serial(by_shift=10)
	print('New:\n', ter.string)
	test_string3 = "HETATM    2  PB  GTP X  10     -24.259 -10.797 -22.646  1.00 42.59           P "
	hetatm = pdb_atom_record(test_string3)
	print('Original:\n', hetatm.string)
	hetatm.update_serial(by_shift=10)
	print('New:\n', hetatm.string)
	hetatm.update_recordName("HET") # fake record for testing
	print('recordName updated:\n', hetatm.string)
