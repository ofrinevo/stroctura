from chimera import runCommand as rc
from chimera.tkgui import saveReplyLog
#from TEMPy.ScoringFunctions import ScoringFunctions
#from TEMPy.MapParser import MapParser
import sys
import os
import math
import numpy
import subprocess

NUMOFRESULTS=10
NUMOFFIELDS=12
DEBUG = 1

def get_principal_axes(map_file):
    rc("open "+map_file)
    rc("measure inertia #0")
    saveReplyLog("log.txt")
    log= open("log.txt", 'r')
    for line in log:
        if "v3" in line:
            symVec=line.split()
        elif "center" in line:
            centerLine=line.split()
            break
    rc("close all")
    for i in range(2,5):
        symVec[i]=float(symVec[i])
        centerLine[i]=float(centerLine[i])
    vecTup= (symVec[2], symVec[3], symVec[4])
    centerTup= (centerLine[2], centerLine[3], centerLine[4])
    
    return [vecTup, centerTup] #array of the two tupples
    #TODO get it into a file or whatever



def powerfit(map_file,res,monomer_file, angle = 10, output_path = '',
             num_of_results = NUMOFRESULTS):
    command = "powerfit "
    command += map_file+" "
    command += str(res)+" "
    command += monomer_file+" "
    if output_path:
        command += "-d "
        command += output_path+" "
    if num_of_results:
        command += "-n"
        command += str(num_of_results)+" "
    command += "-a "
    command += str(angle)
    os.system(command)


"""
    outputs a string with the turn command of the form
    turn symmetry_axis deg center symmetry_center modles #models_index
    deg - turning degree
    symmetry_axis - a tuple of the form (x,y,z)
    symetry_center -  a tuple of the form (x,y,z)
    model_index - the index of the model
"""
def get_turn_command(deg,symetry_axis,symetry_center,model_index):
    command = "turn "
    command += str(symetry_axis).replace(' ','')[1:-1] #get only x,y,z out of (x,y,z)
    command += " "
    command += str(deg) #rotation degree
    command += " center "
    command += str(symetry_center).replace(' ','')[1:-1] #get only x,y,z out of (x,y,z)
    command += " models #"
    command += str(model_index)
    return command

"""
monomer_file - location of the original monomer
fit_file- location of the fitted monomer to the map by powerfit
fits_dir- the directory with the fit-files
N - number of symmetry
symetry_axis - a tuple of the form (x,y,z)
symetry_center -  a tuple of the form (x,y,z)
"""
def cyclic_shift(monomer_file,fit_file,fits_dir,N,symetry_axis = (0,0,1) ,symetry_center = (0,0,0)):
    os.system("mkdir cycledFits")
    rc("open " + monomer_file)
    rc("cd "+fits_dir)
    transArray=[]
    for i in range(N):
        print fit_file
        rc("open " + fit_file)
        if i == 0: 
	    rc("match #" + str(i+1) + " #0 showMatrix true move false")
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i+1) #changed for output
        if DEBUG:
            print command
        rc(command)
	rc("match #" + str(i+1) + " #0 showMatrix true move false")
        #FIX
    filePath= "after"+fit_file
    rc("close 0")
    rc("combine #0,#1 close 1")
    rc("cd ..")
    saveReplyLog("transformsLog.txt")
    transform= findTranform(fit_file) 
	
    return transform
		
def findTranform(fit_file):
	firstMonomer=True
	transFile= open("transformsLog.txt", 'r')
	rotationsMats=[]
	translationArray=[]
	checkLine= fit_file+ " opened"
	for line in transFile:
		if checkLine in line:
			if firstMonomer==False:
        			transFile.next() # the turn command isn't made on the first monomer
			transFile.next()
			transFile.next()
			transFile.next()
			firstRow=transFile.next().split()
			secRow=transFile.next().split()
			thirdRow=transFile.next().split()
			rotationsMats.append([[firstRow[0],firstRow[1],firstRow[2]],[secRow[0],secRow[1],secRow[2]],[thirdRow[0],thirdRow[1],thirdRow[2]]])
			translationArray.append([firstRow[3], secRow[3], thirdRow[3]])
			firstMonomer=False

	return rotationsMats,translationArray

	

def vector2Matrix(vector):
#ONLY for 3x3
    mat = [[0 for x in range(3)] for y in range(3)]
    for i in range(3):
        mat[i]=vector[i*3:i*3+3]
    return mat

def mat2angles(M, cy_thresh=None):
    M = numpy.asarray(M)
    if cy_thresh is None:
        try:
            cy_thresh = numpy.finfo(M.dtype).eps * 4
        except ValueError:
            cy_thresh = _FLOAT_EPS_4
    r11, r12, r13, r21, r22, r23, r31, r32, r33 = M.flat
    cy = math.sqrt(r33*r33 + r23*r23)
    if cy > cy_thresh: # cos(y) not close to zero, standard form
        z = math.atan2(-r12,  r11) # atan2(cos(y)*sin(z), cos(y)*cos(z))
        y = math.atan2(r13,  cy) # atan2(sin(y), cy)
        x = math.atan2(-r23, r33) # atan2(cos(y)*sin(x), cos(x)*cos(y))
    else: # cos(y) (close to) zero, so x -> 0.0 (see above)
        # so r21 -> sin(z), r22 -> cos(z) and
        z = math.atan2(r21,  r22)
        y = math.atan2(r13,  cy) # atan2(sin(y), cy)
        x = 0.0
    return z, y, x 

def get_powerfit_results(map_file,res,monomer_file):
   # powerfit(map_file,res,monomer_file,output_path=monomer_file[:-4]+"PF")
    infile = open(monomer_file[:-4]+"PF"+"/solutions.out","r")
    #infile = open("fit5j40/solutions.out","r")
    res=[[0 for x in range(NUMOFFIELDS)] for y in range(NUMOFRESULTS)]
    resEuler = [[0 for x in range(6)] for y in range(NUMOFRESULTS)]
    i=0
    while True:
        line=infile.readline()
        line=line.rstrip()
        if not line: break
        splitLine=line.split()
        if splitLine[0].isdigit() and int(splitLine[0])<=10:
            for j in range(12):
                res[i][j]=float(splitLine[j+4])
            for k in range(3):
                resEuler[i][k]=res[i][k]
            resEuler[i][3:]=mat2angles(vector2Matrix(res[i][3:]))
            i=i+1
            if i>NUMOFRESULTS+1:
                break
    
    return resEuler

def get_scoreTempy(protein_map, map_file):
	score = ScoringFunctions()
	target_map=MapParser.readMRC(map_file)
	protein_map=MapParser.readMRC(protein_map)
	
	score.CCC(protein_map,target_map)
	return score

def get_score(map_file, N, density):
    #We make a density map out of the structure of the protein we made in cycle_shift
    rc("molmap #"+ str(0) + " " + str(density))
    rc("open " + str(map_file))
    rc("measure correlation #0 #1") #Assumption: the cycled protein map is model #0 and the map we just opened is model #1
    #again, reading from reply log
    saveReplyLog("log.txt")
    log= open("log.txt", 'r')
    for line in log:
        if "Correlation" in line:
            corrLine=line.split()
            ind= corrLine.index("=")
            score= float((corrLine[ind+1])[:-1])
            score*=1000 #in order to make the score nicer than a number between -1 to 1.
    return score
	
def output(resultsArr):
	out= open("Results.txt", 'w')
	for i in range (len(resultsArr)):
		out.write("\nResult "+ str(i+1) +":\t Fit " + str(resultsArr[i][0])+ (":\n") + "score: " + str(resultsArr[i][1])+ "\n" + 
"transformations:\n") 
		for j in range (len(resultsArr[i][2])):
			out.write("monomer " +str(j+1)+ ":\n")
			out.write("rotations: " +str((resultsArr[i][2])[j]) + "\n")
			out.write("translations: " +str((resultsArr[i][3])[j]) + "\n\n")
	 	
def main(monomer_file,N,map_file):
    tupArr=get_principal_axes(map_file,N)
    symAxis= tupArr[0]
    center= tupArr[1]
    powerfit_results=get_powerfit_results(map_file,3,monomer_file)  
    fitsDirPath= monomer_file[:-4]+"PF"

    results=[]
    for nameFile in os.listdir(fitsDirPath):
        if nameFile[:3]=="fit":
            
            numID=nameFile[4:-4]
            rc("close all")
	    rotation,translation= cyclic_shift(monomer_file,nameFile,fits_dir,N,symAxis,center)
            score=get_score(map_file,N,3)
            results.append([numID,score, rotation, translation])
	    #should be xAngle,yAngle,zAngle=mat2angles(rotation)
	    #results.append([numID,score, [xAngle,yAngle,zAngle], translation])
                                  
    #sort the fits according to its score, so that the best score is first(place 0)	
    sortedResults= sorted(results, key= lambda item:item[1], reverse= True)
    rc("close all")
    output(sortedResults)
    print "Done"

def checkAllNoPowerfit(monomer_file,N,map_file, fits_dir):

    tupArr=get_principal_axes(map_file)
    symAxis= tupArr[0]
    center= tupArr[1]
    results=[]
    for nameFile in os.listdir(fits_dir):
    	print nameFile
    for nameFile in os.listdir(fits_dir):
        if nameFile[:3]=="fit":            
            numID=nameFile[4:-4]
            rc("close all")
	    rotation,translation= cyclic_shift(monomer_file,nameFile,fits_dir,N,symAxis,center)
            score=get_score(map_file,N,3)
            results.append([numID,score, rotation, translation])
              
    #sort the fits according to its score, so that the best score is first(place 0)	
    sortedResults= sorted(results, key= lambda item:item[1], reverse= True)
    rc("close all")
    output(sortedResults)
    print "Done"

			
#main("4dyc.pdb",4,"map.mrc")
checkAllNoPowerfit("4dyc.pdb",4,"map.mrc", "highRes")
        
