from chimera import runCommand as rc
import os

DEBUG = 0

def powerfit(map_file,res,monomer_file, angle = 10, output_path = ''):
    command = "powerfit "
    command += map_file+" "
    command += str(res)+" "
    command += monomer_file+" "
    if output_path:
        command += "-d "
        command += output_path+" "
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
