import numpy as np
import scipy.linalg as lin


def menu():
    welcome_message = """
    *****************
      Huckel Solver
    *****************
   
    Select molecule type by choosing a number from the following
    1. Linear Polyene
    2. Cyclic Polyene
    3. Platonic Solid
    4. Buckminsterfullerene
    """

    print(welcome_message)

    # Ensure one of he options is selected, or let user try again
    mol_type = get_mol_type()

    # Buckminsterfullerene does not have any choice in n
    if mol_type in ("1", "2"):
        n_atoms = get_n_atoms()
        huckel_mat = np.zeros([n_atoms, n_atoms])

    # Call the right funciton to generate Huckel Matrix
    if mol_type == "1":
        huckel_mat = linear_polyene_mat(n_atoms, huckel_mat)
    if mol_type == "2":
        huckel_mat = cyclic_polyene_mat(n_atoms, huckel_mat)
    if mol_type == "3":
        n_atoms, huckel_mat = platonic()
    if mol_type == "4":
        n_atoms, huckel_mat = buckminsterfullerene()

    ## Fill in the remainder of the matrix (in case it has not been done already)
    symmetrize_matrix(n_atoms, huckel_mat)

    ## Diagonalize huckel matrix
    evals, evects = diagonalize(huckel_mat)

    ## Round values to 3 d.p
    evals_rounded = round_values(evals)
    printresults(n_atoms, evals_rounded)


def platonic():

    platonic_question = """
    
    Which platonic solid would you like?

     4. Tetrahedron
     6. Octahedron
     8. Hexahedron
    12. Icosahedron
    20. Dodecahedron

    """

    print(platonic_question)

    ## Limit choice to acceptable responses
    platonic_accepted = False
    while platonic_accepted == False:
        platonic_type = input()
        if platonic_type in ("4", "6", "8", "12", "20"):
            platonic_accepted = True
            n_atoms = int(platonic_type)
        else:
            print("The input must be 4,6,8,12 or 20 Try again")

    ## Call the corresponsing huckel matrix generator
    if platonic_type == "4":
        huckel_mat = tetrahedron()
    if platonic_type == "6":
        huckel_mat = octahedron()
    if platonic_type == "8":
        huckel_mat = hexahedron()
    if platonic_type == "12":
        huckel_mat = icosahedron()
    if platonic_type == "20":
        huckel_mat = dodecahedron()

    return n_atoms, huckel_mat


def tetrahedron():

    huckel_mat = np.zeros([4, 4])
    # All vertices are connected in a tetrahedron, so set all except diagonal values to 1
    for i in range(4):
        for j in range(i + 1, 4):
            huckel_mat[i, j] = 1

    return huckel_mat



def hexahedron():

    ## Having set beta to 1, the huckel matrix becomes the mathemetical quantity known as 
    ## an adjacency matrix. These are available online, or can be calculated with algorithms
    # Adjacency matrix calcuated from openai.com, which iterated over all pairs
    # and calculated if their x, y and z differences summed to 0 or 1. This makes sense,
    # since vertices which are not connected have differences 2 or 3

    huckel_mat = np.matrix([
        [0, 1, 1, 1, 1, 0, 0, 0],
        [1, 0, 1, 1, 0, 1, 0, 0],
        [1, 1, 0, 1, 0, 0, 1, 0],
        [1, 1, 1, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 0]])

    return huckel_mat


def icosahedron():

    # Adjacency matrix text file was downloaded from distanceregular.org
    huckel_mat = np.zeros([12, 12])
    file_name = "icosahedron.am.txt"

    f = open(file_name, "r")
    lines = f.readlines()

    for i in range(12):
        row_string = lines[i].strip()
        for j in range(12):
            huckel_mat[i, j] = row_string[j]

    return huckel_mat


def dodecahedron():

    # Adjacency matrix text file was downloaded from distanceregular.org

    huckel_mat = np.zeros([20, 20])

    file_name = "dodecahedron.am.txt"

    f = open(file_name, "r")
    lines = f.readlines()

    for i in range(20):
        row_string = lines[i].strip()
        for j in range(20):
            huckel_mat[i, j] = row_string[j]

    return huckel_mat


def octahedron():

    # Adjacency matrix text file was downloaded from distanceregular.org
    # This one was small enough to manually enter
    huckel_mat = np.matrix([
        [0, 1, 1, 1, 1, 0],
        [1, 0, 1, 1, 0, 1],
        [1, 1, 0, 0, 1, 1],
        [1, 1, 0, 0, 1, 1],
        [1, 0, 1, 1, 0, 1],
        [0, 1, 1, 1, 1, 0]])

    return huckel_mat


def get_n_atoms():

    ## Ensure only physical values for n accepted
    n_atoms_accepted = False
    while n_atoms_accepted == False:
        print("Choose the nuber of atoms in your molecule")

        try:
            # Try converting input to integer
            n_atoms = int(input())
            if n_atoms > 0:
                n_atoms_accepted = True
            else:
                print("Number of atoms must be a positive integer")
        except:
            "Input is not a number"

    return n_atoms


def get_mol_type():


    # Ensure only allowable values of molecule type are entered
    mol_type_accepted = False

    while mol_type_accepted == False:
        mol_type = input()
        if mol_type in ("1", "2", "3", "4"):
            mol_type_accepted = True
        else:
            print("The input must be 1, 2, 3 or 4. Try again")

    return mol_type


def round_values(evals):
    decimal_places = 3
    evals_rounded = []
    # Iterate over e-values and round them, then append to new set of rounded e-values
    for i in range(len(evals)):
        evals_rounded.append(round(evals[i], decimal_places))

    return(evals_rounded)


def symmetrize_matrix(n_atoms, huckel_mat):
    # Takes upper right half of the matrix and duplicates it
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            huckel_mat[j, i] = huckel_mat[i, j]
    return huckel_mat


def diagonalize(huckel_mat):

    # Scipy handles this very quickly
    evals, evects = lin.eigh(huckel_mat)

    return evals, evects


def linear_polyene_mat(n_atoms, huckel_mat):

    # Diagonal values of alpha already filled in
    # Set RHS adjacent to diagonal to beta

    for i in range(n_atoms - 1):
        huckel_mat[i, i+1] = 1

    return huckel_mat


def cyclic_polyene_mat(n_atoms, huckel_mat):

    # Same as linear except for additional 0, n-1 coupling
    huckel_mat = linear_polyene_mat(n_atoms, huckel_mat)

    huckel_mat[0, n_atoms - 1] = 1

    return huckel_mat


def printresults(n_atoms, evals_rounded):

    ## Define Unicode characters
    alpha = '\u03B1'
    beta = '\u03B2'
    dash = '\U00002500'
    print("The Molecular Orbital Energies are as follows:")
    print()

    # Create a new tuple, which removes any duplicate energies
    # This is made more complicated by the rounding, which treats
    # 0 and -0 as separate values

    collapsed_energies = []
    degeneracies = []

    # Set default degeneracy to 1
    dcount = 1
    i = 0
    while i < n_atoms - 1:
        # Check if e-value is the same as the next one [they are ordered]
        ## If so, don't add it yet, but add one to the degeneracy counter
        e_diff = round(evals_rounded[i]-evals_rounded[i+1], 3)
        if abs(e_diff) == 0.000:
            dcount = dcount + 1
        else:
            collapsed_energies.append(evals_rounded[i])
            degeneracies.append(dcount)
            dcount = 1

        i = i + 1

    ## Ensure top level gets added on if it is non-degenerate
    if dcount == 1:
        collapsed_energies.append(evals_rounded[n_atoms - 1])
        degeneracies.append(1)

  
    ## Diagram is created by defining an overall width based on the largest
    ## Degeneracy, then adding the number of levels, then adding in the
    ## same amount of spaces either side to make them centrally aligned

    max_degen = max(degeneracies)
    diagram_width = 4 * max_degen - 1
    level = dash + dash + dash + " "
    diagram = ""

    for i in range(len(collapsed_energies)):
        energy = alpha + " + " + str(collapsed_energies[i]) + beta

        for j in range(degeneracies[i]):
            diagram = diagram + level
        spaces_required = 0.5 * (diagram_width - len(diagram))
        diagram = nspaces(spaces_required) + diagram + nspaces(spaces_required)
        line = diagram + energy
        print(line)
        diagram = ""







def buckminsterfullerene():

    # MATLAB has a bucky function, which gives an adjacency matrix but displaying it leads to 
    # an unfortunate output format (see bucky.txt). The following code takes
    # the file and places it into the matrix accordingly
    # MATLAB famously indexes from 1 which also needs to be dealt with

    huckel_mat = np.zeros([60, 60])
    matlab_mat = np.zeros([61, 61])

    file_name = "bucky.txt"

    f = open(file_name, "r")
    lines = f.readlines()

    separator = "    "
    for i in range(len(lines)):

        # Firstly, remove the whitespace and 1 digits from the file
        coord = lines[i].split(separator, 1)[0]
        coord = coord.strip()

        # Change the brackets to square brackets
        coord = coord.replace("(", "[")
        coord = coord.replace(")", "]")

        # Assign the relevvant matrix value by executing string as line of code
        command = "matlab_mat" + coord + " = 1"
        exec(command)

    for i in range(1, 61):
        for j in range(1, 61):
            huckel_mat[i-1, j-1] = matlab_mat[i, j]

    n_atoms = 60
    return n_atoms, huckel_mat


def nspaces(n):

    ## Creates a string of spaces of the requested length
    nspaces = ""
    for i in range(round(n)):
        nspaces = nspaces + " "
    return nspaces


def main():

    ## Allows program to be looped until user wishes to quit
    finished = False
    while finished == False:
        menu()
        print()
        print("Enter q to quit or any button to continue")
        if input() == "q":
            finished = True


main()
