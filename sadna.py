from chimera import runCommand as rc

DEBUG = 0

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
def cyclic_shift(monomer_file,N,symetry_axis,symetry_center):
    for i in range(N):
        rc("open " + monomer_file)
        if i == 0:
            continue
        command = get_turn_command(i*360/N,symetry_axis,symetry_center,i)
        if DEBUG:
            print command
        rc(command)
