###############################################################
####      GOTO DATABASE ACCESS COMMAND LINE PROGRAM        ####
#### FOR LIGHT CURVE GENERATION AND SCIENCE FRAME DOWNLOAD ####
####  WRITTEN BY ROSS MCWHIRTER, LATEST VERSION: MAY 2021  ####
###############################################################

# Import the required packages. Some are base packages, others need to be pip installed (see .sh scripts in this directory).

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

# Get and store the current time at script call to display the final runtime at the end of the program run.

now = datetime.now(timezone.utc)

startTime = time.time()

# Fetch the input file using argument parser.

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

# Parse the arguments and assign to variables.

args = vars(ap.parse_args())

input_file = args['input']
output_folder = args['output']
existing = args['existing']
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

# As the dates must be in JD for GOTO but called as MJD, convert the dates to JD.

start_date = start_date + 2400000.5

end_date = end_date + 2400000.5

# Set the current directory as the work directory.

cwd = os.getcwd()

# Grab the GOTO account information from a local file.
# This is in a try loop so if it fails the user can be informed in the exception.

try:

    pwd_file = open(os.path.join(cwd, "GOTOinfo.txt"), "r")

    # The dbname, port, username, password and host are all stored in separate lines in GOTOinfo.txt.
    # If they cannot be read in correctly, return an error. We don't bother specifically identifying the fault.

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

# Get the GOTO server download information from a local file if getframes is True.
# This follows a similar ruleset as the GOTOinfo.txt file above.

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

# Check the arguments for any issues.
# Start with checking to see if the --input argument is valid.
# If an existing light curve, it should be a directory otherwise it should be a file.

if existing == False:

    input_present = os.path.isfile(os.path.join(cwd, input_file))

else:

    input_present = os.path.isdir(os.path.join(cwd, input_file))

# If the input file/directory doesn't exist, we exit with an error message.

if input_present == False:
    print("Input target file/directory not found. Exiting...")
    quit()

# Otherwise, we continue with the program.

else:

    # This section contains some if statements to print information to the console on the user selections.

    if existing == False:

        print("Input target file accepted. Beginning queries now...")

        if output_folder is not None:

            print("Output directory selected. Light curve data will be saved into this directory.")

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

            print("Output directory selected. Light curve data will be saved into this directory.")

        else:

            print("No output directory defined. Light curve data will be saved into the existing light curve directory.")

            output_folder = input_file

        print("No other input arguments from the script call will be used except for ignoring the overwrite check.")

# Check to see if science frames download is requested. If it is, read out to the command line and set filepath to True if needed.

if getframes == True:

    print("Science frames requested. These frames will be download into a subdirectory named 'frames' in the output directory.")

    if filepath == False:

        print("The --filepath argument must be set to access the science frames. Setting this argument automatically...")

        filepath = True

# Make a function containing a dialogue for overwriting existing files.
# This function will cause the program to quit if a user does not reply Y or y to a command line input.

def overwrite():

    print("Output directory already exists. If light curves with the same target name are present they will be overwritten.")

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
# It takes size_val as an argument, a number of potential download frames used only for a print statement.
# It returns a Boolean value indicating if the user consents to the download or not based on a Y or y user input.

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

# Make directory named light_curves for the output data files if it does not exist and no output directory is supplied. Also make output directory if it does not exist.
# If the output directory is already present, the overwrite() function above is triggered to confirm that the user consents to the old directory being overwritten.
# If noowchk (no overwrite check) is True from the input arguments, the overwrite() function is bypassed and the directory is automatically overwritten.

if existing == False or output_folder != input_file:

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

# Establish connection to the GOTO database using the Psycopg2 package.
# We record the start time of this connection attempt to resolve issues with connections later.

print("Connecting to the GOTO photometry database...")

print("If this attempt takes longer than 10 seconds it will result in a 'timeout expired' error.")

connectTime = time.time()

# Try to use the variables read in earlier from the GOTOinfo.txt file to create a connection to the GOTO database.
# Put in a try-except so errors can be read out to the console for user debugging.
# We assume that the errors are due to the GOTOinfo.txt file rather than bugs in this script or an issue with the database server.

try:
    
    dbconn = pg2.connect(dbname=GOTO_dbname, port=GOTO_port, user=GOTO_user, password=GOTO_pwd, host=GOTO_host, connect_timeout=10)

except pg2.DatabaseError as err:

    print("The connection to the GOTO Database has resulted in an error: " + str(err))

    print("Check your GOTOinfo.txt file. Alternatively, the database server might be down. Exiting...")

    quit()

# Display the time required for the database connection to monitor consistency.

finconnectTime = str(time.time() - connectTime)

print("Connection successful in " + finconnectTime + " seconds.")

# Create a new 'cursor' an interface to pass SQL commands to the database.

cursor = dbconn.cursor()

# Read in the input file list of sources or scan the input directory for light curves to build an input file.

if existing == False:

    tar_list = pd.read_csv(input_file, sep = " ", header = None).values

else:

    # Check the input directory for light curve files.

    lcfiles = []

    for lcfile in os.listdir(os.path.join(cwd, input_file)):
        if lcfile.endswith(".dat"):
            lcfiles.append(lcfile)

    # Read the argument line in from the existing light curves in the input directory. If it fails... skip that light curve.

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

# Run a loop to query each row in the target list, be it from an input file, or the previous directory light curve read operation.

for i, target in enumerate(tar_list[:,0]):

    # Report to the user which object is currently being processed.

    print("Processing object #" + str(i+1) + " named: " + str(target) + ".")

    # If using an existing set of light curves, get the query arguments from the tar_list array, not the script arguments.

    if existing == True:

        # start_date = tar_list[i,7] + 2400000.5 (Commented out because of new method of getting start_date from light curve below)

        end_date = 2466154.5

        quality = int(tar_list[i,8])

        get_non_det = tar_list[i,5]

        imgid = tar_list[i,9]

        filepath = tar_list[i,10]

        hash_head = lc_header_hash[i]

        last_data_row = pd.read_csv(os.path.join(cwd, input_file, target + ".dat"), skiprows=18, delimiter=" ").tail(1)

        # Grab the start date by inspecting the existing light curve JD column (which can be either column 0 or 1 depending on the imgid argument).
        # We add a tiny time instant onto the most recent observation JD value of the existing light curve.
        # This is shorter than the exposure time of the next frame so no observations should be missed.

        try:

            if imgid == False:

                start_date = last_data_row.values[0][0] + 0.00001

            else:

                start_date = last_data_row.values[0][1] + 0.00001

        except:

            # If the above operation fails (usually because the input light curve had no observations), just take the start date from the light curve header.

            start_date = tar_list[i,6]

            print("No observations detected in existing light curve. Setting start date to " + str(start_date) + " MJD.")

            start_date = start_date + 2400000.5

    # Use the Psycopg2 cursor created earlier to send a select query to the GOTO database to get the details of all images with detections within a pre-defined radius of the input coordinates.
    # To do this the detections table must be joined to the images table using the image_id which we perform using the aliases in 'FROM image AS img JOIN detection AS det ON img.id = det.image_id'.
    # The conditions are then defined after the WHERE by requesting better quality than the input argument, between the start_date and end_date and using a custom cone-search function on the
    # Right Ascension, Declination and cone-search radius as defined in either the input file or the existing light curves in the input directory.
    # We return a set of table columns as shown at the start of the query aftet the SELECT command. The format is <table>.<column>, for example img.jd is the jd (time) column from the image table.
    
    cursor.execute("SELECT img.id, img.jd, det.mag, det.mag_err, det.mag_instr, det.mag_err_instr, det.mag_err_calib, img.exptime, img.filter, img.instrument, img.limmag5, img.fwhm, img.ncoadds, img.ndets, det.ra, det.dec, img.ccdid, img.runid, img.quality, img.quality_template, img.filepath FROM image AS img JOIN detection AS det ON img.id = det.image_id WHERE img.quality <= " + str(quality) + " and img.jd BETWEEN " + str(start_date) + " and " + str(end_date) + " and q3c_radial_query(det.ra, det.dec, " + str(tar_list[i,1]) + ", " + str(tar_list[i,2]) + ", " + str(tar_list[i,3]/3600.) + ")")

    # We use the cursor.fetchall() function to grab the data returned by this query which we add to the data variable.
    
    data = cursor.fetchall()

    # If no observations are present in the 'detections' variable (data), a Boolean is created indicating this for future concatenation of dataframes.

    if data == []:

        empty_data = True

    else:

        empty_data = False

        # If observations are present in the data variable, turn data into a Pandas dataframe and assign column names to each of the columns (as called in the query above).

        data = pd.DataFrame(data)

        data.rename(columns={0: "id", 1: "jd", 2: "mag", 3: "mag_err", 4: "mag_instr", 5: "mag_err_instr", 6: "mag_err_calib", 7: "exptime", 8: "filter", 9: "instrument", 10: "limmag5", 11: "fwhm", 12: "ncoadds", 13: "ndets", 14: "ra", 15: "dec", 16: "ccdid", 17: "runid", 18: "quality", 19: "quality_template", 20: "filepath"}, inplace = True)

        # Reset the index of the data dataframe.

        data.reset_index(drop = True, inplace = True)

    # We now check to see if the user has requested non-detections.

    if get_non_det == True:

        # In this case, we need to get the images where our RA, DEC would have been in the frame polygon (detections and non-detections).
        # This requires a second SQL query which we again send with cursor.execute().
        # As only one table is queried (image), we do not need to give the table name before the columns as in the previous query.
        # We use a similar set of conditions as in the previous query but we now use the q3c_poly_query function.
        # This function will determine if the RA and DEC are within the field of view of any image in the image table by comparing them with the FoV polygon stored for each image.

        cursor.execute("SELECT id, jd, exptime, filter, instrument, limmag5, fwhm, ncoadds, ndets, ccdid, runid, quality, quality_template, filepath FROM image WHERE quality <= " + str(quality) + " and jd BETWEEN " + str(start_date) + " and " + str(end_date) + " and q3c_poly_query(" + str(tar_list[i,1]) + " , " + str(tar_list[i,2]) + ", ('{' || regexp_replace(fov::TEXT, '[()]', '', 'g') || '}')::DOUBLE PRECISION[])")

        # Cursor.fetchall() returns the data from this query which we assign to the variable datan (short for non-detection data).

        datan = cursor.fetchall()

        # Again we check to see if datan contains no observations and assign a Boolean indicating this state.

        if datan == []:

            empty_datan = True

        else:

            empty_datan = False

            # Turn datan into a Pandas dataframe and assign appropriate columns as seen in the second query.

            datan = pd.DataFrame(datan)

            datan.rename(columns={0: "id", 1: "jd", 2: "exptime", 3: "filter", 4: "instrument", 5: "limmag5", 6: "fwhm", 7: "ncoadds", 8: "ndets", 9: "ccdid", 10: "runid", 11: "quality", 12: "quality_template", 13: "filepath"}, inplace = True)

            # Reset the index of the datan dataframe.

            datan.reset_index(drop = True, inplace = True)

        # With data (detections data) and datan (non-detections data) collected, these two dataframes must be concatenated.

        if empty_data == False:

            if empty_datan == False:

                # If empty_data and empty_datan are false, we have both detections and non-detections and must concatenate them together.

                # First, we must remove the detections from the non-detections as all detected observations will also be included in the non-detections query.
                # We can do this by selecting all rows from datan where the row in datan is NOT present in data. These are non-detection images.
                # We assign this set of pure non-detections to a new dataframe datann.

                datann = datan[~datan['id'].isin(data['id'].values)]

                # Concatenate the detections (data) and non-detections (datann) leaving NaNs in the empty cells (performed automatically in the concatenation).
                # Name this combined dataframe of detections and non-detections 'result'.

                result = pd.concat([data, datann])

                # In some situations a detection pass by the pipeline may be done multiple times on the same image resulting in identical rows. Drop the duplicates.
                # We apply this to the data dataframe as well as the result dataframe purely to maintain a count of the detections and non-detections for the user.

                data = data.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                # Read out the number of detections (given by the number of rows in data) and the number of non-detections (number of rows in result minus the number of detections).

                print("Number of detections: " + str(len(data.index)) + ". Number of non-detections: " + str(len(result.index) - len(data.index)) + ".")

            else:

                # In this section of the if statements, there were detections and no non-detections. Therefore just set result to be equal to data and drop the duplicates as above.

                result = data

                data = data.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                # Read out the number of detections (number of rows in data) and non-detections is automatically zero in this case.

                print("Number of detections: " + str(len(data.index)) + ". Number of non-detections: 0.")

        else:

            if empty_datan == False:

                # Here empty_data is True, and empty_datan is False: we only have non-detections.
                # In this case things are a bit more complicated, we cannot just assign result = datan as datan has less columns than the expected data/result dataframe.
                # The solution is to create an empty dataframe named data with the appropriate column names expected from the data dataframe.
                # This empty dataframe can then be concatenated with the datan dataframe giving the expected results.

                data = pd.DataFrame(columns=["id", "jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template", "filepath"])

                result = pd.concat([data, datan])

                result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

                # Read out the number of detections (in this case zero) and the number of non-detections which is the number of rows in the result dataframe.

                print("Number of detections: 0. Number of non-detections: " + str(len(result.index)) + ".")

            else:

                # In this case there were no detections or non-detections so read this out to console and continue.

                print("Number of detections: 0. Number of non-detections: 0.")

    else:

        # This section involves the operation of the script when non-detections are not requested.
        # In this case, the only dataframe is data which is either empty or not denoted by empty_data.

        if empty_data == True:

            # There were no detections and non-detections were not queried. Read this out to the console.

            print("Number of detections: 0.")

        else:

            # As there is no datan in this mode, result = data, drop the duplicates and read out the number of rows in result to the user as detections.

            result = data

            result = result.drop_duplicates(subset = ["jd", "mag", "mag_err", "mag_instr", "mag_err_instr", "mag_err_calib", "exptime", "filter", "instrument", "limmag5", "fwhm", "ncoadds", "ndets", "ra", "dec", "ccdid", "runid", "quality", "quality_template"])

            print("Number of detections: " + str(len(result.index)) + ".")

    # We have now completed the collection of the detections and (if requested) non-detections. Now they must be prepared for writing to the output light curves files.
    # First, we use a try-except to perform some processing steps on the result dataframe. In the case were there were no results, this will force an error.
    # If an error is forced, we continue the script but we define a Boolean variable named no_new_results to inform later parts of the program that result is empty (or isn't defined).

    try:

        no_new_results = False

        # Sort the values by the Julian Date of the observation (as is typical in light curves). In the case of two identical JDs, sort by the image id.

        result = result.sort_values(by = ["jd", "id"])

        result.reset_index(drop = True, inplace = True)

        # Add the MJD column into position 2 (the third column after image id and JD) for easier readability to the user.

        result.insert(2, 'mjd', result['jd'] - 2400000.5)

        # Grab and store the unique values in the filepath column. These will be used for processing the filepath column and for potential science frame downloading later.

        all_filepaths = result.filepath.unique()

        # Prepare the data for writing by converting the Pandas dataframe into a numpy array using the df.values function.

        result = result.values

    except:

        # The 'result' dataframe is empty or doesn't exist causing the above code to fail. Set no_new_results = True.

        no_new_results = True

    # Update the filepath information for both writing the output light curves for use in downloading the science frames.
    # This is required as the filepath column contains a filepath for the internal database, not the one needed for external access.
    # This is only necessary if there are no new results (obviously because there are no filepaths in an empty dataframe).

    updated_filepaths = []

    if no_new_results == False:

        # Loop through all the values in result[:,21] (the filepath column). We need to call it result[:,21] as it is now a numpy array with no column names.

        for path_num, path_value in enumerate(result[:,21]):

            # Split the filepath by the directories in the path.

            path_split = path_value.split("/")

            # Check to see if the target file is in the old or new pipeline.
            # This is performed by seeing if the image date (grabbed from the filepath as nightly data is assigned into a folder with the format YYYY-MM-DD)
            # subtracted from the new_pipeline_date, the date the pipeline swapped from the old to the new version on 2020-07-01
            # after converting the string-based date into the iso format using the datetime package retaining the number of integer days.

            new_pipeline_date = datetime.fromisoformat('2020-07-01')

            target_file_date = datetime.fromisoformat(path_split[5])

            time_diff_days = (new_pipeline_date - target_file_date).days

            # Using the filepath format required for external download, create the external filepath and append to a list of filepaths for later light curve writing.

            rebuilt_path = path_split[5] + "/final/" + path_split[7]

            # Use either pipeline for the newer images (younger than 2020-07-01) and pipeline1 for the older images.

            if time_diff_days < 0:

                updated_filepaths.append('/pipeline/' + rebuilt_path)

            else:

                updated_filepaths.append('/pipeline1/' + rebuilt_path)

    # Get the filepath for the output light curve for the currently processed target by combining the output_folder with <target name>.dat.
    
    data_filepath = os.path.join(cwd, output_folder, target + ".dat")

    # If using an existing light curve, read the data for that light curve now. Also get the existing light curve start time for the file writing.

    if existing == True:

        # Existing light curve filepath.

        input_data_path = os.path.join(cwd, input_file, target + ".dat")

        # Use a try-except to attempt to read in the existing light curve. If we fail to parse the light curve, set old_lc_parse = False and continue.
        # The targeted failure mode is in the event the existing light curve file contains a header but no observations: old_lc_parse informs the program of this case.
        # Currently this is NOT logged for the user anywhere and is a target for a future iteration on this program.

        try:

            old_lc = pd.read_csv(input_data_path, header=None, skiprows=19, delimiter=" ")

            old_lc_parse = True

        except:

            old_lc_parse = False

        # SPECIAL NOTE: As I add in these additional descriptive comments I notice that the code as-is could potentially lead to crashes here.
        # The existing try-except is designed to bypass read issue caused by EMPTY input light curves (header present, but no rows).
        # In the event the user manually places .dat files in the existing light curve directory which do not conform to the expected format,
        # the code will fail, possibly before this point but I note it here to make future developers aware of this potential issue.

        # Read in the column names from the exisiting light curve (despite being called old_lc_header (Pandas terminology) this differs from the
        # description of header used elsewhere here which refers to the full block of text above the light curves in the output data files.

        old_lc_header = pd.read_csv(input_data_path, header=None, skiprows=18, nrows = 1, delimiter=" ").values[0]

        # With an existing light curve, hash_head has been set from the existing light curve header. We need to clip off this hash to set the column names.

        if hash_head == True:

            old_lc_header = old_lc_header[1:]

        # Use Pandas to format the existing light curve dataframe old_lc by setting column names and sorting by JD and then image id.
        # Then immediately use the df.values function to turn the dataframe into a numpy array.

        if old_lc_parse == True:

            old_lc.columns = old_lc_header

            old_lc = old_lc.sort_values(by = ["jd", "id"]).values

        else:

            # If the existing light curve is empty, we just create an empty dataframe with column names from the existing light curve file.

            old_lc = pd.DataFrame(columns=old_lc_header).values

        # Get the start date from the header of the file if there are no observations in the light curve.

        start_date = tar_list[i,6] + 2400000.5

    # Open the output light curve data file with the filepath defined earlier for writing.
    
    datafile = open(data_filepath, "w")

    # Rather ugly way of defining the column names for writing the output file column names depending on whether the image id, file path, both or neither columns are requested in the arguments.

    header_full = "id jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template filepath\n"

    header_noimgid = "jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template filepath\n"

    header_nofilepath = "id jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template\n"

    header_min = "jd mjd mag mag_err mag_instr mag_err_instr mag_err_calib exptime filter instrument limmag5 fwhm ncoadds ndets ra dec ccdid runid quality quality_template\n"

    # Put the arguments of the run in as comments creating the header of the light curve file. This format is required for the 'existing' argument to function.

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("# " + target + " " + str(tar_list[i,1]) + " " + str(tar_list[i,2]) + " " + str(tar_list[i,3]) + " " + str(existing) + " " + str(get_non_det) + " " + str(start_date - 2400000.5) + " " + str(end_date - 2400000.5) + " " + str(int(quality)) + " " + str(imgid) + " " + str(filepath) + "\n")

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("# List of important arguments to GOTO get_lightcurves script for this output:\n")

    datafile.write("# UTC of most recent script call: " + now.strftime("%d/%m/%Y %H:%M:%S") + "\n")

    datafile.write("# Source name: " + target + "\n")

    datafile.write("# Source RA: " + str(tar_list[i,1]) + "\n")

    datafile.write("# Source DEC: " + str(tar_list[i,2]) + "\n")

    datafile.write("# Cone search radius in arcsecs: " + str(tar_list[i,3]) + "\n")

    datafile.write("# Use existing light curve directory: " + str(existing) + "\n")

    datafile.write("# Get non-detections: " + str(get_non_det) + "\n")

    datafile.write("# Starting MJD of light curves: " + str(start_date - 2400000.5) + "\n")

    datafile.write("# Ending MJD of light curves: " + str(end_date - 2400000.5) + "\n")

    datafile.write("# Quality flag limit for GOTO images: " + str(int(quality)) + "\n")

    datafile.write("# Return image ids: " + str(imgid) + "\n")

    datafile.write("# Return image file paths: " + str(filepath) + "\n")

    datafile.write("# ----------------------------------------------------------------------------------------------------------------------------------------------\n")

    datafile.write("#\n")

    # This block of if statements is used to correctly format the column names based on whether the user wants the column names 'commented out' with a hash (hash_head = True)
    # and which of the columns are being returned (image id, file path, both or neither).

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

    # Write in the existing light curve to the output light curve if it is being used as the input for this program run.

    if existing == True:

        # Iterate through the observations rows in the existing light curve, format them correctly and write each row as a new line in the output file.
        # Whilst iterating through the observation rows, make sure the format of quality_out, quality_template_out are consistent.
        # This is required as the output must be a string, but the input can be either integer or string.
        # Without this segment, the integers would be transformed into floats before becoming strings which would be wrong in the output file.

        for obs in range(old_lc.shape[0]):

            try:

                quality_out = str(int(old_lc[obs,19]))

            except:

                quality_out = str(old_lc[obs,19])

            try:

                quality_template_out = str(int(old_lc[obs,20]))

            except:

                quality_template_out = str(old_lc[obs,20])

            # This is where the observations are written as new lines to the output file taking care to write the correct elements depending on the requested columns.

            if imgid == True and filepath == True:

                datafile.write(str(old_lc[obs,0]) + " " + str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(old_lc[obs,21]) + "\n")

            elif imgid == False and filepath == True:

                datafile.write(str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(old_lc[obs,21]) + "\n")

            elif imgid == True and filepath == False:

                datafile.write(str(old_lc[obs,0]) + " " + str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

            else:

                datafile.write(str(old_lc[obs,1]) + " " + str(old_lc[obs,2]) + " " + str(old_lc[obs,3]) + " " + str(old_lc[obs,4]) + " " + str(old_lc[obs,5]) + " " + str(old_lc[obs,6]) + " " + str(old_lc[obs,7]) + " " + str(old_lc[obs,8]) + " " + str(old_lc[obs,9]) + " " + str(old_lc[obs,10]) + " " + str(old_lc[obs,11]) + " " + str(old_lc[obs,12]) + " " + str(old_lc[obs,13]) + " " + str(old_lc[obs,14]) + " " + str(old_lc[obs,15]) + " " + str(old_lc[obs,16]) + " " + str(old_lc[obs,17]) + " " + str(old_lc[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

    # At this point the light curve contains a header, column names and (if using existing = True) the existing light curve observations.
    # Now the new observations queried in this program run are written to the output light curve files.

    # If there are new results to write to the output file as denoted by the no_new_results = False, write them to the file.

    if no_new_results == False:

        # As with the existing light curve, iterate through the rows in the 'result' numpy array containing the newly queried results.
        # Correct the quality_out and quality_template_out as shown above to maintain a correct integer to string conversion.
    
        for obs in range(result.shape[0]):

            try:

                quality_out = str(int(result[obs,19]))

            except:

                quality_out = str(result[obs,19])

            try:

                quality_template_out = str(int(result[obs,20]))

            except:

                quality_template_out = str(result[obs,20])


            # Write the requested columns from the result numpy array to the output light curve file.

            if imgid == True and filepath == True:

                datafile.write(str(result[obs,0]) + " " + str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(updated_filepaths[obs]) + "\n")

            elif imgid == False and filepath == True:

                datafile.write(str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + " " + str(updated_filepaths[obs]) + "\n")

            elif imgid == True and filepath == False:

                datafile.write(str(result[obs,0]) + " " + str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")

            else:

                datafile.write(str(result[obs,1]) + " " + str(result[obs,2]) + " " + str(result[obs,3]) + " " + str(result[obs,4]) + " " + str(result[obs,5]) + " " + str(result[obs,6]) + " " + str(result[obs,7]) + " " + str(result[obs,8]) + " " + str(result[obs,9]) + " " + str(result[obs,10]) + " " + str(result[obs,11]) + " " + str(result[obs,12]) + " " + str(result[obs,13]) + " " + str(result[obs,14]) + " " + str(result[obs,15]) + " " + str(result[obs,16]) + " " + str(result[obs,17]) + " " + str(result[obs,18]) + " " + quality_out + " " + quality_template_out + "\n")
        
    # With all the writing operations complete for this target, close the datafile preventing additional writing by the program this call.

    datafile.close()

    # With this loops target light curve output file written, the program now looks to see if the user has requested to download the science frames of this target.
    # Use RSYNC to download the science frames if getframes is True and there are new frames to get.

    if getframes == True and no_new_results == False:

        # Check to see if it is okay to proceed when noszchk is False (the size check function has not been skipped).
        # If the user decides to proceed with the download, the check_if_size_okay function returns True which is assigned to the variable proceed_with_download.

        if noszchk == False:

            proceed_with_download = check_if_size_okay(len(all_filepaths))

        else:

            proceed_with_download = True

        # If the user has selected to download the files, commence the next block of code in the if statement.

        if proceed_with_download == True:

            # Create an output directory for the fits files
            # This output directory is named frames within the output_folder directory.
            # Currently the user cannot provide a separate output folder for science frames. This may be a target for a future program iteration.

            output_frames_path = os.path.join(cwd, output_folder, "frames")

            if not os.path.exists(output_frames_path):
                os.makedirs(output_frames_path)

            # Create a subdirectory for the given target. The target frames will be stored in a directory named <target name> in the frames directory.

            output_frames_tarpath = os.path.join(output_frames_path, target)

            if not os.path.exists(output_frames_tarpath):
                os.makedirs(output_frames_tarpath)

            # Inform the user that the download is commencing.

            print("Commencing download of requested science frames...")

            # For each unique filepath from the queried detections and (if requested) non-detections, stored in the all_filepaths list
            # we again perform the operations on the filepath to transform it from an internal to an external filepath.
            # NOTE: This is a prime target for code refactoring, there is NO NEED for this to be performed twice, despite it being relatively
            # fast compared to most other operations in this program. This has occurred because the correction of the filepaths to the external
            # version for the output light curve data files occurred AFTER this section for downloading the science frames.
            # The output light curves filepaths were changed afterward to allow a manual http download of the files from the GOTO file server.

            for path_num, path_value in tqdm(enumerate(all_filepaths)):

                path_split = path_value.split("/")

                # Check to see if the target file is in the old or new pipeline

                new_pipeline_date = datetime.fromisoformat('2020-07-01')

                target_file_date = datetime.fromisoformat(path_split[5])

                time_diff_days = (new_pipeline_date - target_file_date).days

                rebuilt_path = path_split[5] + "/final/" + path_split[7]

                if time_diff_days < 0:

                    # Try to run the rsync to download it to the given directory. This is performed by passing a shell command out of the Python runtime.

                    try:

                        check_output(['rsync -vrdt rsync://' + GOTO_rsync_address + '/pipeline/' + rebuilt_path + ' --password-file rsync.pwd ' + output_frames_tarpath], shell = True)

                    except:

                        continue

                else:

                    # Try to run the rsync to download it to the given directory. This is performed by passing a shell command out of the Python runtime.

                    try:

                        check_output(['rsync -vrdt rsync://' + GOTO_rsync_address + '/pipeline1/' + rebuilt_path + ' --password-file rsync.pwd ' + output_frames_tarpath], shell = True)

                    except:

                        continue

    # Prior to the next target, delete the larger dataframes and numpy arrays created during this iteration.
    # As these variables are not always created, a set of if statements and try-excepts determine which variables are to be deleted.

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

# Once the loop has iterated through every target in the input target list (from either file or directory), the Psycopg2 cursor and database connection are closed.

cursor.close()
    
dbconn.close()

# Finally, the elapsed time is calculated by getting the current clock time and subtracting the time recorded at program start.
# This is written out to the command line and the program exits.

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
