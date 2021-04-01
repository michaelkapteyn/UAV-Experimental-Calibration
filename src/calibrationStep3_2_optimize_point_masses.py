### NOTE: This final will not run!
### The following functions are not included as our model is proprietary.
### The following (self explanatory) functions would need to be implemented in order for this script to interact with a given structural model.

# modify_material_properties_in_structural_FEA_model(Emultiplier)
# modify_point_masses_in_structural_FEA_model(point_mass_dict)
# frequencies = run_structural_FEA_model_and_return_natural_frequencies()

import os
import json
import numpy as np
import scipy.optimize
import csv

TARGET_FREQS = [5.6492, 29.9699] # computed in STEP3_1
ZETA = [0.0281, 0.0717]   # computed in STEP3_1
def modify_masses_and_do_solve(mass):
    print(">>\n>>\n>>\n>>\n>>\n>>\n")
    print("Running solve with masses:")

    if mass[0] < 0:
        mass[0] = 0.0

    mass3 = 0.472-2.0*mass[0]
    if  mass3 < 0:
        mass3 = 0.0

    point_mass_dict = { "servo_1" : mass[0],
                        "servo_2" : mass[0],
                        "pitot_tube" : mass3}

    print(point_mass_dict)
    modify_point_masses_in_structural_FEA_model(point_mass_dict)
    frequencies = run_structural_FEA_model_and_return_natural_frequencies()
    frequencies = [frequencies[0], frequencies[2]]

    print("Target frequencies: ", TARGET_FREQS)
    print("Model frequencies: ", frequencies)


    error = 0.5*(np.abs(TARGET_FREQS[0] - frequencies[0])/TARGET_FREQS[0] + np.abs(TARGET_FREQS[1] - frequencies[1])/TARGET_FREQS[1])

    print("Mean Relative Error: ", error)
    return error

def final_modify_masses_and_do_solve(mass):
    print(">>\n>>\n>>\n>>\n>>\n>>\n")
    print("Running solve with masses:")

    if mass[0] < 0:
        mass[0] = 0.0

    mass3 = 0.472-2.0*mass[0]
    if  mass3 < 0:
        mass3 = 0.0

    point_mass_dict = { "servo_1" : mass[0],
                        "servo_2" : mass[0],
                        "pitot_tube" : mass3}
    print(point_mass_dict)

    modify_point_masses_in_structural_FEA_model(point_mass_dict)
    frequencies = run_structural_FEA_model_and_return_natural_frequencies()
    frequencies = [frequencies[0], frequencies[2]]
    return frequencies


if __name__ == '__main__':
    # read in the e samples generated in Step 3_1 and saved in e_samples.csv
    Emultiplierslist = []
    with open('datafiles/step3/e_samples.csv', newline='') as inputfile:
        for row in csv.reader(inputfile):
            Emultiplierslist.append(row[0])

    # loop over the e samples and run the mass optimization routine
    for i in range(100):
        print("\n")
        print("Starting iteration ", i)

        # get the Young's modulus scale factor from the data saved in e_samples.csv
        Emultiplier = Emultiplierslist[i]

        print("E multiplier= ", Emultiplier)
        print("\n")

        # modify Young's modulus and hg commit+push
        modify_material_properties_in_structural_FEA_model(Emultiplier)
        mass_opt = scipy.optimize.fmin(func=modify_masses_and_do_solve, x0=[0.169], xtol=0.001)

        # scipy.optimize.fmin returns optimal masses, but has not necessarily computed the corresponding frequencies. do that here
        f_opt = final_modify_masses_and_do_solve(mass_opt)
        mass3 = 0.472-2.0*mass_opt[0];

        # write data to file
        with open('datafiles/step3/massoptimizationresults_new.csv', 'a') as the_file:
            the_file.write(str(Emultiplier)+","+str(TARGET_FREQS[0])+","+str(TARGET_FREQS[1])+","+str(ZETA[0])+","+str(ZETA[1])+","+str(mass_opt[0])+","+str(mass_opt[0])+","+str(mass3)+","+str(f_opt[0])+","+str(f_opt[1])+'\n')

        print("Optimized masses:")
        point_mass_dict_opt = { "servo_1" : mass_opt[0],
                                "servo_2" : mass_opt[0],
                                "pitot_tube" : mass3}
        print(point_mass_dict_opt)
        print("\n")
        print("Optimized frequency mode 1: ", f_opt[0])
        print("Optimized frequency mode 2: ", f_opt[0])
        print("\n")
        print("End iteration ",i)
