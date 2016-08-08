from chimera import runCommand as rc
from chimera.tkgui import saveReplyLog
import sys
import os

NUMOFRESULTS=10
DEBUG = 0

def get_principal_axes(map_file,num):
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
    for i in range(N):
        rc("open " + monomer_file)
        if i == 0:
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i)
        if DEBUG:
            print command
        rc(command)

def get_powerfit_results(map_file,res,monomer_file):
    powerfit(map_file,res,monomer_file,output_path=monomer_file[:-4]+"PF")
    infile = open(monomer_file[:-4]+"PF"+"/solutions.out","r")
    infile=next(infile)
   # for line in infile:
        #TODO
    
def main(monomer_file,N,map_file):
    #(axis,center)=get_principal_axes(map_file,N)
	tupArr=get_principal_axes(map_file,N)
	symAxis= tupArr[0]
	center= tupArr[1]
    #powerfit_results=get_powerfit_results(map_file,3,monomer_file)
    #for result in powerfit_results:
       # cyclic_shift(monomer_file,N,symAxis,center)
        
