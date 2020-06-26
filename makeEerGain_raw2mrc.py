#!/software/EMAN2/bin/python2.7
  
# Copyright (c) 2019 by FEI Company

from __future__ import print_function
from __future__ import division



from past.utils import old_div
import struct
from EMAN2 import *
from numpy import *
import sys,os
import getopt
import xml.etree.ElementTree as ET

def printHelp():
	progname = "makeEerGain.py"
	usage = progname + """  Convert eer gain in raw format to mrc
	
	Options:
		-i, --raw \t<input raw file>
		-d, --defect \t<input SensorDefects.xml>
		-o, --gain \t<output gain file (with mrc extenstion)>
	        -h, --help
	
	Examples:
		makeEerGain.py -i gain_post_ec_eer.raw -d SensorDefects.xml -o gain4k.mrc
	
	"""
	print(usage)
	exit(0)

def getFilenames(argv):
	rawFile = ''
	xmlFile = ''
	mrcFile = ''

	try:
		opts, args = getopt.getopt(argv,"hi:d:o:",["help","raw=","defect=","gain="])
	except getopt.GetoptError:
		printHelp()
	for opt, arg in opts:
		if opt in ('-h','--help'):
			printHelp()
			sys.exit()
		elif opt in ("-i", "--raw"):
			rawFile = arg
		elif opt in ("-o", "--gain"):
			mrcFile = arg
		elif opt in ('-d', '--defect'):
			xmlFile = arg
		else:      
			print(opt+'unrecognized option')
	
	#check inputs
	if rawFile=='':
		print("input raw file name is missing")
		printHelp()
	elif xmlFile=='':
		print("input defect file is missing")
		printHelp()
	elif mrcFile=='':
		print("output mrc file name is missing")
		printHelp() 	
		
	return rawFile, xmlFile, mrcFile
	
def main():
	
	if len(sys.argv)<2:
		printHelp()
		
	rawFile, xmlFile, mrcFile = getFilenames(sys.argv[1:])
	
	print("Reading defects from "+xmlFile)
	defectPoint=[]
	defectRow=[]
	defectCol=[]
	#defectArea=[]
	tree = ET.parse(xmlFile)
	root = tree.getroot()
	for elem in root.iter():
		if elem.tag=='point':
			defectPoint.append(map(int,elem.text.split(',')))
		elif elem.tag=='col':
			defectCol.append(map(int,elem.text.split('-')))
		elif elem.tag=='row':
			defectRow.append(map(int,elem.text.split('-')))
	
	
	
	gainList=[]
	print("Reading "+rawFile)
	with open(rawFile,'rb') as a:
		a.seek(49)
		for i in range(0,4096*4096):
			s = struct.unpack('i32',a.read(4))
			gainList.append(s[0])
	gain64=array(gainList).reshape((4096,4096)).astype('float64')

	print("Masking pixels in "+xmlFile)
	for col in defectCol:
		for j in range(col[0]-1,col[1]+2):
			#print('col: %d'%j)
			gain64[:,j]=0
	for row in defectRow:
		for i in range(row[0]-1,row[1]+2):
			#print('row: %d'%i)
			gain64[i,:]=0
	for point in defectPoint:
		gain64[point[1],point[0]]=0
	#print("Masking edge")
	gain64[0,:]=0
	gain64[:,0]=0
	gain64[4095,:]=0
	gain64[:,4095]=0
	if True:
		ind = where(gain64>0)
		mm=gain64[ind].mean()
		ss=gain64[ind].std()
		dd=16
		print("Electrons in the gain: %f"%mm)
		ind2=where(gain64>mm+ss*dd)
		#ind3=where(gain64<mm-ss*dd)
		#print("Extra defects: hot %d cold %d"%(len(ind2[0]),len(ind3[0])))
		#gain64[ind3]=0
		print("Masking %d hot pixels beyond %d sigma"%(len(ind2[0]),dd))
		gain64[ind2]=0
	print("Saving "+mrcFile)
	e = EMNumPy.numpy2em(gain64.astype('float32'))
	e.write_image(mrcFile)	


if __name__ == "__main__":
	main()
