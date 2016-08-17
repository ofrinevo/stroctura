from chimera import runCommand as rc
from chimera.tkgui import saveReplyLog
#from TEMPy.ScoringFunctions import ScoringFunctions
#from TEMPy.MapParser import MapParser
import sys
import os

NUMOFRESULTS=10
NUMOFFIELDS=12
DEBUG = 0

def get_principal_axes(map_file,num=0):
    rc("open "+mapFile)
   # rc("measure inertia #"+ str(num))
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
    vecTup= (symVec[2], symVec[3], symVec[4])
    centerTup= (centerLine[2], centerLine[3], centerLine[4])
    print [vecTup, centerTup]
	
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
    print command
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
monomer_file - location of the monomer
N - number of symmetry
symetry_axis - a tuple of the form (x,y,z)
symetry_center -  a tuple of the form (x,y,z)
"""
def cyclic_shift(monomer_file,N,symetry_axis = (0,0,1) ,symetry_center = (0,0,0)):
	#os.system("mkdir cycledFits")
    for i in range(N):
        rc("open " + monomer_file)
        if i == 0:
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i)
        if DEBUG:
            print command
        rc(command)
		
		#filePath= "/cycledFits/"+numID+".pdb"
	#	rc("write format pdb 0" + filePath) #saves the protein to the new folder with appropriate name
	#	rc("close #0") #that way the model ID will always be 0 
		
		

def get_powerfit_results(map_file,res,monomer_file):
    powerfit(map_file,res,monomer_file,output_path=monomer_file[:-4]+"PF")
    infile = open(monomer_file[:-4]+"PF"+"/solutions.out","r")
   # infile = open("fit5j40/solutions.out","r")
    #infile=next(infile)
    res=[[0 for x in range(NUMOFFIELDS)] for y in range(NUMOFRESULTS)]
    i=0
    for line in infile:
        splitLine=line.split()
        if DEBUG:
            print splitLine
        if i<=NUMOFRESULTS+1 and splitLine[0]!="#rank" and int(splitLine[0])<=10:
            for j in range(12):
                res[i][j]=splitLine[j+4]
            i=i+1
    return res

def get_scoreTempy(protein_map, map_file):
	score = ScoringFunctions()
	target_map=MapParser.readMRC(map_file)
	protein_map=MapParser.readMRC(protein_map)
	
	score.CCC(protein_map,target_map)
	return score

def get_score(map_file):
	#We make a density map out of the structure of the protein we made in cycle_shift
	rc("molmap #0 3")
	rc("open " + str(map_file))
	rc("measure correlation #1 #0") #Assumption: the cycled protein map is model #0 and the map we just opened is model #1
	#again, reading from reply log
	saveReplyLog("log.txt")
	log= open("log.txt", 'r')
    for line in log:
        if "Correlation" in line:
            corrLine=line.split()
            break
	ind= corrLine.index("=")
	score= float((corrLine[ind+1])[:-1])
	score*=1000 #in order to make the score nicer than a number between -1 to 1.
	
	return score
	
	

def main(monomer_file,N,map_file):
    #(axis,center)=get_principal_axes(map_file,N)
	tupArr=get_principal_axes(map_file,N)
	symAxis= tupArr[0]
	center= tupArr[1]
    #powerfit_results=get_powerfit_results(map_file,3,monomer_file)
	   
	fitsDirPath= monomer_file[:-4]+"PF"
	results=[]
	for nameFile in os.listdir(fitsDirPath):
		if nameFile[:3]=="fit":
			numID=nameFile[4:-4]
			cyclic_shift(nameFile,N,symAxis,center)
		#	fileProtein= "/cycledFits/"+numID+".pdb"
			score=get_score(map_file)
			results.append([numID,score])
			
	#sort the fits according to its score, so that the best score is first(place 0)		
	sorted(results, key= lambda item:item[1], reverse= True) 
			
	
        
