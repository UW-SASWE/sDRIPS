import configparser
import subprocess

# Define your lists of start_date and save_data_loc values
start_dates = [ '2022-04-07', '2022-04-14','2022-04-23'] 
save_data_locs = [ '../Data_2022_April_07/', '../Data_2022_April_14/','../Data_2022_April_23/']

# Path to the config file and main script
config_file_path = '../Config_files/Script_Config.ini'
main_script_path = '../Codes/SEBAL_ET_Canals_Framework_BWDB_CMD.py'
conda_env_name = 'canals_env' 
working_directory = 'C:/Users/skhan7/OneDrive - UW/Desktop/Research/PhD/Chapter2/Codes//'
# Create a function to update the config file
def update_config(start_date, save_data_loc):
    # Read the config file
    config = configparser.ConfigParser()
    config.read(config_file_path)
    
    # Update the start_date and save_data_loc
    config['Date_Running']['start_date'] = start_date
    config['Save_Data_Location']['save_data_loc'] = save_data_loc
    
    # Write the updated config back to the file
    with open(config_file_path, 'w') as configfile:
        config.write(configfile, space_around_delimiters=False)

def run_main_script():
    # Command to activate the conda environment and run the main script
    command = f'conda run -n {conda_env_name} python {main_script_path}'
    result = subprocess.run(command, shell=True, check=True, cwd=working_directory,stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)

for start_date, save_data_loc in zip(start_dates, save_data_locs):
    print(f'Running script for start date: {start_date}')
    update_config(start_date, save_data_loc)
    run_main_script()