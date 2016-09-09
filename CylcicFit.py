from chimera import runCommand as rc
from chimera.tkgui import saveReplyLog
import sys
import os
import math
import numpy
import subprocess
import VolumeViewer
from Segger import segment_dialog as SD

#Constants
ROTATION_DEGREE = 20
NUMOFRESULTS = 10
NUMOFFIELDS = 12
_FLOAT_EPS_4 = numpy.finfo(float).eps * 4.0




def get_principal_axes(map_file):
    """
    Recives a map file, returns it's three axes and the center of it.
    The return value is: [axis 1, axis 2, axis 3, center]
    Each component is a tuple of three numbers
    """
    rc("open "+map_file)
    rc("measure inertia #0")
    saveReplyLog("log.txt")
    log= open("log.txt", 'r')
    for line in log:
        if "v1" in line:
            lineVec1=line.split()
        elif "v2" in line:
            lineVec2=line.split()
        elif "v3" in line:
            lineVec3=line.split()
        elif "center" in line:
            centerLine=line.split()
            break
    rc("close all")
    for i in range(2,5):
        lineVec1[i]=float(lineVec1[i])
        lineVec2[i]=float(lineVec2[i])
        lineVec3[i]=float(lineVec3[i])
        centerLine[i]=float(centerLine[i])
    vec1= (lineVec1[2], lineVec1[3], lineVec1[4])
    vec2= (lineVec2[2], lineVec2[3], lineVec2[4])
    vec3= (lineVec3[2], lineVec3[3], lineVec3[4])
    centerTup= (centerLine[2], centerLine[3], centerLine[4])
    
    return [vec1, vec2, vec3, centerTup]


def runSegment(map_file):
    """Segments the map from the given path and saves the results in
    the directory segmentationsDir.
    """
    chimeraMap = VolumeViewer.open_volume_file(map_file)[0] #opens the map
    dialog = SD.Volume_Segmentation_Dialog() #creates a segmentor object
    trsh = dialog.Segment() #segments the map

    #code for saving the results
    dmap = dialog.SegmentationMap()
    smod = dialog.CurrentSegmentation()
    regs = smod.selected_regions()
    if len(regs)==0 :
            regs = smod.regions
    dir = os.path.dirname ( dmap.data.path )
    fprefix = os.path.splitext ( dmap.name ) [0]
    fname = fprefix + "_region_%d.mrc"
    path = "segmentationsDir/"+fname
    os.system("rm -r segmentationsDir ; mkdir segmentationsDir")
    for reg in regs : #save each result in a different file
        dialog.SaveRegsToMRC ( [reg], dmap, path % (reg.rid,) )
    #end of code for saving the results
        
    rc("close all")

def fitToSegment(monomer_file):
    """
    saves all of the map's segments
    """
    os.system("rm -r fitDir ; mkdir fitDir")
    i=0
    for segmentedMapFile in os.listdir( "segmentationsDir" ) :
        rc( "open " + monomer_file)
	v= VolumeViewer.open_volume_file("segmentationsDir/" + segmentedMapFile)[0]
	d=v.data
	#We move the monomer close to the segmented-region map to improve powerfit, so we measure the position of the monomer in atomic coordinates
	#measure the position of the segmented-map in grid indices, convert it to atomic corrdinates, and move the monomer in the coordinates 		
	#difference axis, then running fitmap to fit the monomer to the region-map when they are closed to each other
	rc( "measure center #0")
	rc( "measure center #1")
	saveReplyLog("log.txt")
	log= open("log.txt", 'r')
	lines=log.readlines()
	centerLine1=lines[len(lines)-2]
	centerLine0=lines[len(lines)-3] 
	line0Array=centerLine0.split()
	line1Array=centerLine1.split()
	monomerCenter= (float(line0Array[-3][1:-1]), float(line0Array[-2][:-1]), float(line0Array[-1][:-1]))
	segmentCenter= (float(line1Array[-3][1:-1]), float(line1Array[-2][:-1]), float(line1Array[-1][:-1]))
	x,y,z= d.ijk_to_xyz(segmentCenter) 
	moveAxis= (x-monomerCenter[0], y-monomerCenter[1], z-monomerCenter[2])
	rc("move "+str(moveAxis[0]) +","+ str(moveAxis[1]) + "," + str(moveAxis[2]) + " models #0")
        rc( "fitmap #0 #1" )
        rc("write format pdb #0 fitDir/model" + str(i)+".pdb")
        rc("close all")
        i+=1
    
def get_turn_command(deg,symetry_axis,symetry_center,model_index):
    """outputs a string with the turn command of the form
    turn symmetry_axis deg center symmetry_center modles #models_index
    deg - turning degree
    symmetry_axis - a tuple of the form (x,y,z)
    symetry_center -  a tuple of the form (x,y,z)
    model_index - the index of the model
    """
    command = "turn "
    command += str(symetry_axis).replace(' ','')[1:-1] #get only x,y,z out of (x,y,z)
    command += " "
    command += str(deg) #rotation degree
    command += " center "
    command += str(symetry_center).replace(' ','')[1:-1] #get only x,y,z out of (x,y,z)
    command += " models #"
    command += str(model_index)
    return command


def cyclic_shift(fit_file,fits_dir,N,symetry_axis ,symetry_center,numOfAxis):
    """
    fit_file- location of the fitted monomer to the map by powerfit
    fits_dir- the directory with the fit-files
    N - number of symmetry
    symetry_axis - a tuple of the form (x,y,z)
    symetry_center -  a tuple of the form (x,y,z)
    numOfAxis- the index of symmetry axis (1/2/3)
    """
    rc("cd "+fits_dir)
    for i in range(N):
        rc("open " + fit_file)
        if i == 0: 
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i) #changed for output
        rc(command)
    filePath= "after"+fit_file+" axis "+str(numOfAxis)
    os.system("mkdir cycledFits")
    rc("cd cycledFits")
    rc("combine #0,#1 close 1")
    rc("write format pdb #"+ str(N) + " " + filePath+".pdb") 
    rc("cd ..")
    rc("cd ..")


def cyclicForOutput(monomer_file,fit_file,fits_dir,N,symetry_axis,symetry_center):

    rc("open " + monomer_file) 
    rc("cd " + fits_dir)
    for i in range(N):
        rc("open " + fit_file)
        if i == 0: 
            print "Calculating transformations and translations " + fit_file #important! Do not remove
            rc("match #" + str(i+1) + " #0 showMatrix true move false")
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i+1) #changed for output
        rc(command)
        print "Calculating transformations and translations " + fit_file #important! Do not remove
        rc("match #" + str(i+1) + " #0 showMatrix true move false")
    rc("cd ..")
    saveReplyLog("log.txt")
    transform= findTransform(fit_file)
    rc("close all") 

    return transform
     
def findTransform(fit_file):
    """
    Given a fit, returns the set of transformations and torartions
    """
    transFile= open("log.txt", 'r')
    rotationsMats=[]
    translationArray=[]
    checkLine= "Calculating transformations and translations " + fit_file
    for line in transFile:
        if checkLine in line:
            transFile.next()
            transFile.next()
            transFile.next()
            firstRow=transFile.next().split()
            secRow=transFile.next().split()
            thirdRow=transFile.next().split()
            rotationsMats.append([[float(firstRow[0]),float(firstRow[1]),float(firstRow[2])],[float(secRow[0]),float(secRow[1]),float(secRow[2])],[float(thirdRow[0]),float(thirdRow[1]),float(thirdRow[2])]])
            translationArray.append([float(firstRow[3]), float(secRow[3]), float(thirdRow[3])])

    return rotationsMats,translationArray

    

def vector2Matrix(vector):
    """
    Given a vector with 9 numbers, returns a 3x3 matrix of it
    """
    mat = [[0 for x in range(3)] for y in range(3)]
    for i in range(3):
        mat[i]=vector[i*3:i*3+3]
    return mat

def mat2angles(M, cy_thresh=None):
    """
    Converts a rotation matrix into 3 euler angels
    Takned from the interweb
    """
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

def get_score(map_file, N, density):
    """
    Given a map, the number of symm and the density, computes and returns the
    correlation between the original map and the current tested model.
    Assums that the current model is open.
    """
    rc("molmap #"+str(N)+" " + str(density))
    rc("open " + str(map_file))
    rc("measure correlation #0 #"+str(N)) #Assumption: the cycled protein map is model #0 and the map we just opened is model #1
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
    """
    Taking care of creating an output file
    """
    out= open("Results.txt", 'w')
    for i in range (len(resultsArr)):
        out.write("\nResult "+ str(i+1) +":\t Fit " + str(resultsArr[i][0])+ (":\n") + "score: " + str(resultsArr[i][1])+ "\n" + 
        "transformations:\n") 
        for j in range (len(resultsArr[i][2])):
            out.write("monomer " +str(j+1)+ ":\n")
            out.write("rotations: " +str(mat2angles((resultsArr[i][2])[j])) + "\n")
            out.write("translations: " +str((resultsArr[i][3])[j]) + "\n\n")

                
def main(monomer_file,N,map_file,density):
    #calculates symm axes and center
    tupArr=get_principal_axes(map_file)
    axes=tupArr[0:3]
    center= tupArr[3]
    #all closed
    
    #gets segmentation
    runSegment(map_file)
    #all closed
    
    tempScore=[0,0,0]
    results=[]
    fitToSegment(monomer_file)
    fits_dir="fitDir"
    for fits in os.listdir(fits_dir):
        if fits[-4:]!=".pdb":
            continue
        numID=fits
        for i in range (3):
            cyclic_shift(fits,fits_dir,N,axes[i],center,i)
            tempScore[i]=get_score(map_file,N,density)
            rc("close all")
        maxIndex=numpy.argmax(tempScore)
        rotationsMat,translationArray= cyclicForOutput(monomer_file,fits,fits_dir,N,axes[maxIndex],center)
        results.append([numID,tempScore[maxIndex], rotationsMat, translationArray])
          
    #sort the fits according to its score, so that the best score is first(place 0) 
    sortedResults= sorted(results, key= lambda item:item[1], reverse= True)
    rc("close all")
    output(sortedResults)
    print "Done"


######ENTER RUNNING LINE HERE:
#args are: monomer file, symm number, map file, map res
main("4dyc.pdb",4,"4dyc_res6.mrc",6)
