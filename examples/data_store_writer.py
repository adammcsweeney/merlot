# returns two csv files for collection of data from trajectory optimisation runs
# merlot segments earth return mission concepts. headers related to outbound and inbound leg trajectory data
# note: this module will clear the contents of the csv files, removing data previously stored

# TODO ==== create means to store run data from separate sessions in individual files, so data is not overwritten
# TODO ==== add thrust profile, and spacecraft ephemris data to collection
# TODO ==== add conjunction/planetary angle data
# TODO ==== create title for run data based on inputs

import csv

# OUTBOUND DATA

def outbound_datafile_make(planet,departure_year,engine,no_engines,launcher):

    # create file name sequence
    s = "-"
    name = ('outbound',str(planet),str(departure_year),str(engine),str(no_engines),str(launcher))
    filename = s.join(name)+'.csv'

    # create path to data store location
    z = '\\'
    path = r"C:\Users\andrew mcsweeney\Documents\projects\merlot\run_data"
    file = (path,filename)
    path = z.join(file)

    # create file
    with open(path, 'w+', newline='') as csvfile:

        # write headers
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(  # RUN DATA
            ['ID', 'Date', 'No_segments', 'Run_time',
             # INPUTS
             'Origin', 'Destination', 'Power', 'Ion_eng', 'No_eng', 'Launch_vehicle',
             # TIMELINE DATA
             'Launch_date', 'helio_tof', 'target_arr_date', 'SpiralIn_tof', 'orb_arr_date',
             # DEPARTURE AND ARRIVAL EXCESS VELOCITY
             'C3', 'Dep_ex_v', 'Arr_ex_v',
             # MASS DATA
             'Launch_mass', 'helio_Xe', 'target_arr_mass', 'SpiralIn_Xe', 'orb_arr_mass',
             # DELTA-V
             'helio_delta-v', 'SpiralIn_delta-v',
             # SPIRAL THRUST PARAMETERS
             'thr_at_target', 'isp_at_target',
             # HELIO DISTANCE
             'max_helio_dist', 'min_helio_dist'])

def return_datafile_make(planet,departure_year,engine,no_engines):

    # create file name sequence
    s = "-"
    name = ('return',str(planet),str(departure_year),str(engine),str(no_engines))
    filename = s.join(name)+'.csv'

    # create path to data store location
    z = '\\'
    path = r"C:\Users\andrew mcsweeney\Documents\projects\merlot\run_data"
    file = (path,filename)
    path = z.join(file)

    with open(path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(  # RUN DATA
            ['ID', 'Date', 'No_segments', 'Run_time',
             # INPUTS
             'Origin', 'Destination', 'Power', 'Ion_eng', 'No_eng', 'Return_mass'
             # TIMELINE DATA
            'orb_dep_date', 'SpiralOut_tof', 'target_dep_date',
             'helio_tof', 'Earth_arr_date',
             # DEPARTURE AND ARRIVAL EXCESS VELOCITY
             'Dep_ex_v', 'Arr_ex_v',
             # MASS DATA
             'mass_in_orbit', 'SpiralOut_Xe', 'target_dep_mass', 'helio_Xe', 'arr_mass',
             # DELTA-V
             'SpiralOut_delta-v', 'helio_delta-v',
             # SPIRAL THRUST PARAMETERS
             'thr_at_target', 'isp_at_target',
             # HELIO DISTANCE
             'max_helio_dist', 'min_helio_dist'])






# planet = 'mars'
# depyear = 2026
# engine = 'T6'
# no = 1
# launcher = 'falcon9'
#
# outbound_data_write(planet,depyear,engine,no,launcher)


# with open('merlot_outbound_data.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile, delimiter = ',')
#     writer.writerow( # RUN DATA
#                     ['ID', 'Date', 'No_segments', 'Run_time',
#                      # INPUTS
#                      'Origin', 'Destination', 'Power', 'Ion_eng', 'No_eng','Launch_vehicle',
#                      # TIMELINE DATA
#                      'Launch_date','helio_tof', 'target_arr_date', 'SpiralIn_tof', 'orb_arr_date',
#                      # DEPARTURE AND ARRIVAL EXCESS VELOCITY
#                      'C3','Dep_ex_v', 'Arr_ex_v',
#                      # MASS DATA
#                      'Launch_mass', 'helio_Xe', 'target_arr_mass', 'SpiralIn_Xe', 'orb_arr_mass',
#                      # DELTA-V
#                      'helio_delta-v','SpiralIn_delta-v',
#                      # SPIRAL THRUST PARAMETERS
#                      'thr_at_target', 'isp_at_target',
#                      # HELIO DISTANCE
#                      'max_helio_dist', 'min_helio_dist'])
#
# # RETURN DATA
#
with open('merlot_return_data.csv', 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    writer.writerow(  # RUN DATA
                    ['ID', 'Date', 'No_segments', 'Run_time',
                     # INPUTS
                     'Origin', 'Destination','Power', 'Ion_eng', 'No_eng','Return_mass'
                     # TIMELINE DATA
                     'orb_dep_date', 'SpiralOut_tof', 'target_dep_date', 'helio_tof', 'Earth_arr_date',
                     # DEPARTURE AND ARRIVAL EXCESS VELOCITY
                     'Dep_ex_v', 'Arr_ex_v',
                     # MASS DATA
                     'mass_in_orbit', 'SpiralOut_Xe', 'target_dep_mass', 'helio_Xe', 'arr_mass',
                     # DELTA-V
                     'SpiralOut_delta-v', 'helio_delta-v',
                     # SPIRAL THRUST PARAMETERS
                     'thr_at_target', 'isp_at_target',
                     # HELIO DISTANCE
                     'max_helio_dist', 'min_helio_dist'])