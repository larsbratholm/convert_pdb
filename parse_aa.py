import sys

data = {}
residue = ""
for line in open(sys.argv[1]).readlines():
    if "RESIDUE" in line:
        residue = line.split()[1]
        data[residue] = {}
    elif "CONECT" in line:
        tokens = line.split()
        if tokens[1] in ["OXT","HXT", "H2"] or (residue == "PRO" and tokens[1] == "H"):
            continue
        #if len(tokens[1]) < 4:
        #    tokens[1] = " " + tokens[1]
        ## do this twice instead of while loop
        #if len(tokens[1]) < 4:
        #    tokens[1] += " "
        #if len(tokens[1]) < 4:
        #    tokens[1] += " "
        data[residue][tokens[1]] = []
        for atom in tokens[3:]:
            if atom == "H2" or (residue == "PRO" and tokens[1] == "H"):
                atom = "pC"
            elif atom == "OXT":
                atom = "nN"
            elif atom == "HXT":
                atom = "nH"
            #if len(atom) < 4:
            #    atom = " " + atom
            ## do this twice instead of while loop
            #if len(atom) < 4:
            #    atom += " "
            #if len(atom) < 4:
            #    atom += " "
            data[residue][tokens[1]].append(atom)

# create 1 and 2 bond atoms to predict the atom position
estimation_data = {}
for residue in data.keys():
    estimation_data[residue] = {}
    for main_atom in data[residue].keys():
        estimation_data[residue][main_atom] = []
        neighbours = data[residue][main_atom]
        for n_index, n in enumerate(neighbours):
            needed_atoms = [n]
            if n[0] in ["p","n"]:
                m = n[1:]
                prefix = n[0]
            else:
                m = n
                prefix = ""
            # hack due to pro having no H
            if residue == "PRO" and m == "H":
                neighbour_atoms = data["ALA"][m][:]
            else:
                neighbour_atoms = data[residue][m][:]
            #if residue == "GLY" and main_atom == "HA2":
            #    print m, neighbour_atoms
            #    quit()
            for index_2, atom_2 in enumerate(neighbour_atoms):
                if prefix != "" and atom_2[0] in ["p","n"]:
                    if main_atom == atom_2[1:]:
                        main_index = index_2
                    neighbour_atoms[index_2] = atom_2[1:]
                else:
                    if main_atom == atom_2:
                        main_index = index_2
                    neighbour_atoms[index_2] = prefix+atom_2
                    #    quit()
            # To be consistent in stereo assignment we have to keep track of odd and even indices
            #if len(main_atom) == 4 and main_atom[:2] in "HH":
            #    print neighbour_atoms, main_index, n_index, needed_atoms
            # Arginine special case
            if main_index % 2 == 1 or (len(main_atom) == 4 and main_atom[:3] == "HH1"):
                needed_atoms.extend(neighbour_atoms[:main_index]+neighbour_atoms[main_index+1:])
            else:
                tmp_list = neighbour_atoms[:main_index]+neighbour_atoms[main_index+1:]
                tmp_list.reverse()
                needed_atoms.extend(tmp_list[:])

            if len(needed_atoms) > 1:
                estimation_data[residue][main_atom].append(needed_atoms)
print estimation_data

quit()
for i in estimation_data.keys():
    for j in estimation_data[i]:
        for k in estimation_data[i][j]:
            print i, j, k



