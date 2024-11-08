import configparser
import subprocess
from datetime import datetime
import sys

# List of dates and corresponding save_data_loc for each iteration
# start_dates = [
#     '2023-04-07',
#     '2023-04-15',
#     '2023-02-18',
#     '2023-01-25',
#     '2023-04-28',
# ]
start_dates = ['2023-04-07','2023-04-15',
    '2023-02-18',
    '2023-01-25',
    '2023-04-23',
]

config_path = '../Config_files/Script_Config.ini'
main_script_path = '../Codes/SEBAL_ET_Canals_Framework_BWDB_CMD.py'

for start_date in start_dates:
    # Read the existing config file
    config = configparser.ConfigParser()
    config.read(config_path)

    # Update the config file with the current start date and formatted save_data_loc
    formatted_date = datetime.strptime(start_date, '%Y-%m-%d').strftime('%Y_%B_%d')
    config['Date_Running']['start_date'] = start_date
    config['Save_Data_Location']['save_data_loc'] = f'../Data_{formatted_date}_Resistance_Correction3/'

    # Write the updated config file
    with open(config_path, 'w') as configfile:
        config.write(configfile)

    # Run the main script with the updated config file and capture the output
    process = subprocess.Popen(['python', main_script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)

    # Display the output in real time
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            sys.stdout.write(output)
            sys.stdout.flush()
    
    # Print any errors in real time
    while True:
        error = process.stderr.readline()
        if error == '' and process.poll() is not None:
            break
        if error:
            sys.stderr.write(error)
            sys.stderr.flush()

    # Ensure the process has completed
    process.wait()
    print(f'Finished for {start_date}')

print('All tasks completed.')
