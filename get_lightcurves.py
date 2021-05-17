import psycopg2 as pg2
import pandas as pd
import numpy as np
import os
import argparse
import signal
import time
import shutil
from datetime import datetime, timezone
from tqdm import tqdm
from subprocess import check_output

now = datetime.now(timezone.utc)

startTime = time.time()

# Fetch the input file using argument parser

ap = argparse.ArgumentParser()

ap.add_argument("--input", type=str, required=True, help="Input target file with coordinates for GOTO query or a directory of existing light curve files.")
ap.add_argument("--output", type=str, required=False, help="Output directory of light curves (ignored if an existing input directory is used).")
ap.add_argument("--existing", required=False, default=False, action='store_true', help="Instead of a target input file, use an existing directory of light curves.")
ap.add_argument("--nondet", required=False, default=False, action='store_true', help="Get the non-detection frames (slow).")
ap.add_argument("--hashhead", required=False, default=False, action='store_true', help="Add a hash to the 'column names' row of the output table files.")
ap.add_argument("--startdate", type=float, required=False, default=0.0, help="Start date of any light curve in MJD.")
ap.add_argument("--enddate", type=float, required=False, default=66154.0, help="Start date of any light curve in MJD.")
ap.add_argument("--quality", type=int, required=False, default=256, help="Quality limit. Set 0 for high quality, 128 for medium quality.")
ap.add_argument("--imgid", required=False, default=False, action='store_true', help="Return the unique GOTO ID of the images.")
ap.add_argument("--filepath", required=False, default=False, action='store_true', help="Return the file paths of the images.")
ap.add_argument("--noowchk", required=False, default=False, action='store_true', help="Skip the overwrite check for the output light curves directory.")
ap.add_argument("--getframes", required=False, default=False, action='store_true', help="Download Science fits frames into a subdirectory in the output light curves directory.")
ap.add_argument("--noszchk", required=False, default=False, action='store_true', help="Skip the file size check for the download of science frames.")

args = vars(ap.parse_args())

input_file = args['input']
output_folder = args['output']
old_file = args['existing']
get_non_det = args['nondet']
hash_head = args['hashhead']
start_date = args['startdate']
end_date = args['enddate']
quality = args['quality']
imgid = args['imgid']
filepath = args['filepath']
noowchk = args['noowchk']
getframes = args['getframes']
noszchk = args['noszchk']

start_date = start_date + 2400000.5

end_date = end_date + 2400000.5

# Set the current directory as the work directory

cwd = os.getcwd()

# Grab the GOTO account information from a local file

try:

    pwd_file = open(os.path.join(cwd, "GOTOinfo.txt"), "r")

    try:

        GOTO_dbname = pwd_file.readline().rstrip()

        GOTO_port = pwd_file.readline().rstrip()

        GOTO_user = pwd_file.readline().rstrip()

        GOTO_pwd = pwd_file.readline().rstrip()

        GOTO_host = pwd_file.readline().rstrip()

    except:

        pwd_file.close()

        print("GOTO account information file has incorrect format. Exiting...")

        quit()

    pwd_file.close()

except:

    print("GOTO account information file missing. Exiting...")

    quit()

# Get the GOTO server download information from a local file if getframes is True


if getframes == True:

    try:

        rsync_file = open(os.path.join(cwd, "GOTOrsync.txt"), "r")

        try:

            GOTO_rsync_user = rsync_file.readline().rstrip()

            GOTO_rsync_host = rsync_file.readline().rstrip()

            GOTO_rsync_port = rsync_file.readline().rstrip()

        except:

            rsync_file.close()

            print("GOTO RSYNC information file has incorrect format.")

            print("Correct this file or re-run this script without the --getframes argument. Exiting...")

            quit()

        rsync_file.close()

    except:

        print("GOTO RSYNC information file missing.")

        print("Create this file or re-run this script without the --getframes argument. Exiting...")

        quit()

    GOTO_rsync_address = GOTO_rsync_user + '@' + GOTO_rsync_host + ':' + str(GOTO_rsync_port)

# Check the arguments for any issues

if old_file == False:

    input_present = os.path.isfile(os.path.join(cwd, input_file))

else:

    input_present = os.path.isdir(os.path.join(cwd, input_file))

if input_present == False:
    print("Input target file/directory not found. Exiting...")
    quit()

else:

    if old_file == False:

        print("Input target file accepted. Beginning queries now...")

        if output_folder is not None:

            print("Output directory selected. Light curve data will be saved into this folder.")

        else:

            print("No output directory defined. Light curve data will be saved into a 'light_curves' directory in the current work directory.")

            output_folder = "light_curves"

        print("Arguments will be added to the headers of the output data files.")

        if get_non_det == True:

            print("Non-detection frames selected. This will take longer as a polygon search is required.")

        if hash_head == True:

            print("The column names row will have a # added to the start for machine readability.")

        if imgid == True and filepath == True:

            print("Returning the Image ID and Filepath columns of the light curves.")

        elif imgid == False and filepath == True:

            print("Returning the Filepath column of the light curves.")

        elif imgid == True and filepath == False:

            print("Returning the Image ID column of the light curves.")

        else:

            print("The Image ID and Filepath columns will not be returned. If these are required, use the appropriate flags.")

        print("Commencing light curves from the MJD: " + str(start_date - 2400000.5))

        print("Final MJD for light curves: " + str(end_date - 2400000.5))

        print("Only observations with a quality flag <= " + str(quality) + " will be selected.")

    else:

        print("Using existing light curve directory for input.")

        print("Arguments will be determined from the light curves in this directory.")

        if output_folder is not None:

            print("Output directory selected. Light curve data will be saved into this folder.")

        else:

            print("No output directory defined. Light curve data will be saved into the existing light curve directory.")

            output_folder = input_file

        print("No other input arguments from the script call will be used except for ignoring the overwrite check.")

# Check to see if science frames are requested. If they are, read out to the command line and set filepath to True if needed.

if getframes == True:

    print("Science frames requested. These frames will be download into a subdirectory named 'frames' in the output folder.")

    if filepath == False:

        print("The --filepath argument must be set to access the science frames. Setting this argument automatically...")

        filepath = True

# Make a dialogue for overwriting existing files.

def overwrite():

    print("Output folder already exists. If light curves with the same target name are present they will be overwritten.")

    while True:
            
        response = input("Continue? Y/N:")

        if response == "Y" or response == "y":

            print("Program will continue with this output directory.")

            break

        elif response == "N" or response == "n":

            print("Program will exit allowing a new output directory to be supplied.")

            quit()

            break

        else:

            print("Incorrect input supplied, please type either Y for yes or N for no.")

# Make a dialogue to shut down science frame download if there are likely to be more than 10 files.

def check_if_size_okay(size_val):

    print("There are " + str(size_val) + " files in the download list. Each file is around 100 MB in size.")

    print("If you select Yes, the download will commence. If you select No, the download will skip.")

    while True:
            
        response = input("Continue? Y/N:")

        if response == "Y" or response == "y":

            print("Program will continue with the download of these files.")

            return True

        elif response == "N" or response == "n":

            print("Program will skip these file downloads.")

            return False

        else:

            print("Incorrect input supplied, please type either Y for yes or N for no.")

# Make folder named light_curves for the output data files if it does not exist and no output folder is supplied. Also make output folder if it does not exist.

if old_file == False or output_folder != input_file:

    if noowchk == False:

        if output_folder is None:

            if not os.path.exists(os.path.join(cwd, "light_curves")):
                os.makedirs(os.path.join(cwd, "light_curves"))
            else:
                status = overwrite()

        else:

            if not os.path.exists(os.path.join(cwd, output_folder)):
                os.makedirs(os.path.join(cwd, output_folder))
            else:
                status = overwrite()

    else:

        print("Skipping output directory overwrite check...")

        if not os.path.exists(os.path.join(cwd, output_folder)):
                os.makedirs(os.path.join(cwd, output_folder))

# Establish connection to the GOTO database

print("Connecting to the GOTO photometry database...")

print("If this attempt takes longer than 10 seconds it will result in a 'timeout expired' error.")

connectTime = time.time()

try:
    
    dbconn = pg2.connect(dbname=GOTO_dbname, port=GOTO_port, user=GOTO_user, password=GOTO_pwd, host=GOTO_host, connect_timeout=10)

except pg2.DatabaseError as err:

    print("The connection to the GOTO Database has resulted in an error: " + str(err))

    print("Check your GOTOinfo.txt file. Alternatively, the database server might be down. Exiting...")

    quit()

finconnectTime = str(time.time() - connectTime)

print("Connection successful in " + finconnectTime + " seconds.")

cursor = dbconn.cursor()

# Read in the input file list of sources or scan the input folder for light curves to build an input file.

if old_file == False:

    tar_list = pd.read_csv(input_file, sep = " ", header = None).values

else:

    # Check the input directory for light curve files.

    lcfiles = []

    for lcfile in os.listdir(os.path.join(cwd, input_file)):
        if lcfile.endswith(".dat"):
            lcfiles.append(lcfile)

    # Read the argument line in from the light curves. If it fails... skip that light curve.

    lcargs = []

    lc_header_hash = []

    for lc in lcfiles:

        try:

            lcargs.append(pd.read_csv(os.path.join(cwd, input_file, lc), header=None, skiprows=1, nrows=1, delimiter=" ").iloc[:,1:].values[0].tolist())

            header_row_0 = pd.read_csv(os.path.join(cwd, input_file, lc), header=None, skiprows=18, nrows = 1, delimiter=" ").values[0].tolist()[0]

            if header_row_0 == "#":

                lc_header_hash.append(True)

            else:

                lc_header_hash.append(False)

        except:

            continue

        tar_list = np.array(lcargs, dtype=object)

# Run a loop to query each row in the target list

for i, target in enumerate(tar_list[:,0]):

    print("Processing object #" + str(i+1) + " named: " + str(target) + ".")

    # If using an existing set of light curves, get the query arguments from the tar_list array, not the script arguments.

    if old_file == True:

        # start_date = tar_list[i,7] + 2400000.5 (Commented out because of new method of getting start_date from light curve below)

        end_date = 2466154.5

        quality = int(tar_list[i,8])

        get_non_det = tar_list[i,5]

        imgid = tar_list[i,9]

        filepath = tar_list[i,10]

        hash_head = lc_header_hash[i]

        last_data_row = pd.read_csv(os.path.join(cwd, input_file, target + ".dat"), skiprows=18, delimiter=" ").tail(1)

        try:

            if imgid == False:

                start_date = last_data_row.values[0][0] + 0.00001

            else:

                start_date = last_data_row.values[0][1] + 0.00001

        except:

            start_date = tar_list[i,6]

            print("No observations detected in existing light curve. Setting start date to " + str(start_date) + " MJD.")

            start_date = start_date + 2400000.5

    # Get the detections joined onto the images
    
    cursor.execute("SELECT img.id, img.jd, det.mag, det.mag_err, det.mag_instr, det.mag_err_instr, det.mag_err_calib, img.exptime, img.filter, img.instrument, img.limmag5, img.fwhm, img.ncoadds, img.ndets, det.ra, det.dec, img.ccdid, img.runid, img.quality, img.quality_template, img.filepath FROM image AS img JOIN detection AS det ON img.id = det.image_id WHERE img.quality <= " + str(quality) + " and img.jd BETWEEN " + str(start_date) + " and " + str(end_date) + " and q3c_radial_query(det.ra, det.dec, " + str(tar_list[i,1]) + ", " + str(tar_list[i,2]) + ", " + str(tar_list[i,3]/3600.) + ")")
    
    data = cursor.fetchall()

    if data == []:

        empty_data = True

    else:

        empty_data = False

        data = pd.DataFrame(data)

        data.rename(columns={0: "id", 1: "jd", 2: "mag", 3: "mag_err", 4: "mag_instr", 5: "mag_err_instr", 6: "mag_err_calib", 7: "exptime", 8: "filter", 9: "instrument", 10: "limmag5", 11: "fwhm", 12: "ncoadds", 13: "ndets", 14: "ra", 15: "dec", 16: "ccdid", 17: "runid", 18: "quality", 19: "quality_template", 20: "filepath"}, inplace = True)

        data.reset_index(drop = True, inplace = True)

    if get_non_det == True:

        # Now we need to get the images where our ra, dec would have been in the frame polygon (detections and non-detections)

        cursor.execute("SELECT id, jd, exptime, filter, instrument, limmag5, fwhm, ncoadds, ndets, ccdid, runid, quality, quality_template, filepath FROM image WHERE quality <= " + str(quality) + " and jd BETWEEN " + str(start_date) + " and " + str(end_date) + " and q3c_poly_query(" + str(tar_list[i,1]) + " , " + str(tar_list[i,2]) + ", ('{' || regexp_replace(fov::TEXT, '[()]', '', 'g') || '}')::DOUBLE PRECISION[])")

        datan = cursor.fetchall()

        if datan == []:

            empty_datan = True

        else:

            empty_datan = False

            datan = pd.DataFrame(datan)

            datan.rename(columns={0: "id", 1: "jd", 2: "exptime", 3: "filter", 4: "instrument", 5: "limmag5", 6: "fwhm", 7: "ncoadds", 8: "ndets", 9: "ccdid", 10: "runid", 11: "quality", 12: "quality_template", 13: "filepath"}, inplace = True)

            datan.reset_index(drop = True, inplace = True)

        if empty_data == False:

            if empty_datan == False:

                # Remove the detections from the non-detections

                datann = datan[~datan['id'].isin(data['id'].values)]

                # Concatenate the detections and non-detections leaving NaNs in the empty cells

                result = pd.concat([data, datann])

                data = data.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                print("Number of detections: " + str(len(data.index)) + ". Number of non-detections: " + str(len(result.index) - len(data.index)) + ".")

            else:

                result = data

                data = data.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                print("Number of detections: " + str(len(data.index)) + ". Number of non-detections: 0.")

        else:

            if empty_datan == False:

                data = pd.DataFrame(columns=["id", "jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template", "filepath"])

                result = pd.concat([data, datan])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                print("Number of detections: 0. Number of non-detections: " + str(len(result.index)) + ".")

            else:

                print("Number of detections: 0. Number of non-detections: 0.")

    else:

        if empty_data == True:

            print("Number of detections: 0.")

        else:

            result = data

            result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

            print("Number of detections: " + str(len(result.index)) + ".")

    try:

        no_new_results = False

        result = result.sort_values(by = ["jd", "id"])

        result.reset_index(drop = True, inplace = True)

        # Finally drop the image id column as it is no longer needed (This operation is not needed, preserved here for reference)

        # result.drop(columns=['id'], inplace = True)

        # Add the MJD column

        result.insert(2, 'mjd', result['jd'] - 2400000.5)

        # Grab and store the unique values in the filepath column

        all_filepaths = result.filepath.unique()

        # Prepare the data for writing

        result = result.values

    except:

        no_new_results = True

    # Fix the filepath information for both writing the output light curves.

    updated_filepaths = []

    if no_new_results == False:

        for path_num, path_value in enumerate(result[:,21]):

            path_split = path_value.split("/")

            # Check to see if the target file is in the old or new pipeline

            new_pipeline_date = datetime.fromisoformat('2020-07-01')

            target_file_date = datetime.fromisoformat(path_split[5])

            time_diff_days = (new_pipeline_date - target_file_date).days

            rebuilt_path = path_split[5] + "/final/" + path_split[7]

            if time_diff_days < 0:

                updated_filepaths.append('/pipeline/' + rebuilt_path)

            else:

                updated_filepaths.append('/pipeline1/' + rebuilt_path)
    
    data_filepath = os.path.join(cwd, output_folder, target + ".dat")

    # If using an existing light curve, read it in now. Also get the existing light curve start time for the file writing.

    if old_file == True:

        input_data_path = os.path.join(cwd, input_file, target + ".dat")

        try:

            old_lc = pd.read_csv(input_data_path, header=None, skiprows=19, delimiter=" ")

            old_lc_parse = True

        except:

            old_lc_parse = False

        old_lc_header = pd.read_csv(input_data_path, header=None, skiprows=18, nrows = 1, delimiter=" ").values[0]

        if hash_head == True:

            old_lc_header = old_lc_header[1:]

        if old_lc_parse == True:

            old_lc.columns = old_lc_header

            old_lc = old_lc.sort_values(by = ["jd", "id"]).values

        else:

            old_lc = pd.DataFrame(columns=old_lc_header)

        # Get the start date from the header of the file if there are no observations in the light curve.

        start_date = tar_list[i,6] + 2400000.5
    
    datafile = open(data_filepath, "w")

    header_full = "id jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template filepath\n"

    header_noimgid = "jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template filepath\n"

    header_nofilepath = "id jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template\n"

    header_min = "jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template\n"

    # Put the arguments of the run in as comments

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("# " + target + " " + str(tar_list[i,1]) + " " + str(tar_list[i,2]) + " " + str(tar_list[i,3]) + " " + str(old_file) + " " + str(get_non_det) + " " + str(start_date - 2400000.5) + " " + str(end_date - 2400000.5) + " " + str(int(quality)) + " " + str(imgid) + " " + str(filepath) + "\n")

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("# List of important arguments to GOTO get_lightcurves script for this output:\n")

    datafile.write("# UTC of most recent script call: " + now.strftime("%d/%m/%Y %H:%M:%S") + "\n")

    datafile.write("# Source name: " + target + "\n")

    datafile.write("# Source RA: " + str(tar_list[i,1]) + "\n")

    datafile.write("# Source DEC: " + str(tar_list[i,2]) + "\n")

    datafile.write("# Cone search radius in arcsecs: " + str(tar_list[i,3]) + "\n")

    datafile.write("# Use existing light curve directory: " + str(old_file) + "\n")

    datafile.write("# Get non-detections: " + str(get_non_det) + "\n")

    datafile.write("# Starting MJD of light curves: " + str(start_date - 2400000.5) + "\n")

    datafile.write("# Ending MJD of light curves: " + str(end_date - 2400000.5) + "\n")

    datafile.write("# Quality flag limit for GOTO images: " + str(int(quality)) + "\n")

    datafile.write("# Return image ids: " + str(imgid) + "\n")

    datafile.write("# Return image file paths: " + str(filepath) + "\n")

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("#\n")

    if hash_head == True:

        if ((imgid == True) and (filepath == True)):

            datafile.write("# " + header_full)

        elif ((imgid == False) and (filepath == True)):

            datafile.write("# " + header_noimgid)

        elif ((imgid == True) and (filepath == False)):

            datafile.write("# " + header_nofilepath)

        else:

            datafile.write("# " + header_min)

    else:

        if ((imgid == True) and (filepath == True)):

            datafile.write(header_full)

        elif ((imgid == False) and (filepath == True)):

            datafile.write(header_noimgid)

        elif ((imgid == True) and (filepath == False)):

            datafile.write(header_nofilepath)

        else:

            datafile.write(header_min)

    # Put in the existing light curve if it is being used as input.

    if old_file == True:

        for obs in range(old_lc.shape[0]):

            try:

                quality_out = str(int(old_lc[obs,19]))

            except:

                quality_out = str(old_lc[obs,19])

            try:

                quality_template_out = str(int(old_lc[obs,20]))

            except:

                quality_template_out = str(old_lc[obs,20])

            if imgid == True and filepath == True:

                datafile.write(str(old_lc[obs,0]) + " " + str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(old_lc[obs,21]) + "\n")

            elif imgid == False and filepath == True:

                datafile.write(str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(old_lc[obs,21]) + "\n")

            elif imgid == True and filepath == False:

                datafile.write(str(old_lc[obs,0]) + " " + str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

            else:

                datafile.write(str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

    if no_new_results == False:
    
        for obs in range(result.shape[0]):

            try:

                quality_out = str(int(result[obs,19]))

            except:

                quality_out = str(result[obs,19])

            try:

                quality_template_out = str(int(result[obs,20]))

            except:

                quality_template_out = str(result[obs,20])

            if imgid == True and filepath == True:

                datafile.write(str(result[obs,0]) + " " + str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(updated_filepaths[obs]) + "\n")

            elif imgid == False and filepath == True:

                datafile.write(str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(updated_filepaths[obs]) + "\n")

            elif imgid == True and filepath == False:

                datafile.write(str(result[obs,0]) + " " + str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

            else:

                datafile.write(str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")
        
    datafile.close()

    # Use RSYNC to download the science frames if getframes is True and there are new frames to get

    if getframes == True and no_new_results == False:

        # Check to see if it is okay to proceed when noszchk is False

        if noszchk == False:

            proceed_with_download = check_if_size_okay(len(all_filepaths))

        else:

            proceed_with_download = True

        if proceed_with_download == True:

            # Create an output directory for the fits files

            output_frames_path = os.path.join(cwd, output_folder, "frames")

            if not os.path.exists(output_frames_path):
                os.makedirs(output_frames_path)

            # Create a subdirectory for the given target

            output_frames_tarpath = os.path.join(output_frames_path, target)

            if not os.path.exists(output_frames_tarpath):
                os.makedirs(output_frames_tarpath)

            print("Commencing download of requested science frames...")

            for path_num, path_value in tqdm(enumerate(all_filepaths)):

                path_split = path_value.split("/")

                # Check to see if the target file is in the old or new pipeline

                new_pipeline_date = datetime.fromisoformat('2020-07-01')

                target_file_date = datetime.fromisoformat(path_split[5])

                time_diff_days = (new_pipeline_date - target_file_date).days

                rebuilt_path = path_split[5] + "/final/" + path_split[7]

                if time_diff_days < 0:

                    # Try to run the rsync to download it to a given directory

                    try:

                        check_output(['rsync -vrdt rsync://' + GOTO_rsync_address + '/pipeline/' + rebuilt_path + ' --password-file rsync.pwd ' + output_frames_tarpath], shell = True)

                    except:

                        continue

                else:

                    try:

                        check_output(['rsync -vrdt rsync://' + GOTO_rsync_address + '/pipeline1/' + rebuilt_path + ' --password-file rsync.pwd ' + output_frames_tarpath], shell = True)

                    except:

                        continue

    # Clear out the junk

    if get_non_det == True:

        if empty_data == False:
            del data
        
        if empty_datan == False:
            del datan

        try:
            del datann
        except:
            continue

        try:
            del result
        except:
            continue

    else:

        if empty_data == False:
            del data

        try:
            del result
        except:
            continue

cursor.close()
    
dbconn.close()

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
