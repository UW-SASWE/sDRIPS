import os
import logging


def clear_data(
    save_data_loc: str,
    clear_et: bool = True,
    clear_precip: bool = True,
    clear_weather: bool = True,
    clear_uploads: bool = True,
) -> None:
    """
    Clear generated files from specified directories.

    Parameters
    ----------
    save_data_loc : str
        Root directory where the project is running.
    clear_et : bool, optional
        Whether to clear evapotranspiration (ET) data such as Penman, SEBAL, and irrigation. Default is True.
    clear_precip : bool, optional
        Whether to clear precipitation data. Default is True.
    clear_weather : bool, optional
        Whether to clear weather data (excluding precipitation). Default is True.
    clear_uploads : bool, optional
        Whether to clear uploaded files from 'uploads/' directory. Default is True.

    Returns
    -------
    None
    """
    logger = logging.getLogger()

    def safe_clear_folder(folder_path: str, filetypes: tuple = ('.zip', '.tif', '')) -> None:
        """Delete files in a folder matching extensions; pass if folder doesn't exist."""
        if not os.path.exists(folder_path):
            return
        try:
            for filename in os.listdir(folder_path):
                if filename.endswith(filetypes) or filetypes == ('',):  # delete all files if empty string passed
                    os.remove(os.path.join(folder_path, filename))
        except Exception as err:
            logger.error(f"Error clearing folder: {folder_path}")
            logger.exception(err)

    if clear_et:
        logger.info("Clearing evapotranspiration data (Penman, SEBAL, Irrigation)...")
        sensors = ['landsat']
        weeks = ['lastweek', 'currentweek']
        params = ['penman', 'sebal', 'irrigation']

        for sensor in sensors:
            for param in params:
                for week in weeks:
                    et_dir = os.path.join(save_data_loc, sensor, param, week)
                    safe_clear_folder(et_dir, filetypes=('.zip', '.tif'))

    if clear_precip:
        precip_dir = os.path.join(save_data_loc, 'precip')
        logger.info("Clearing precipitation data...")
        safe_clear_folder(precip_dir, filetypes=('','.tif','.csv','.txt'))  # adjust as needed

    if clear_weather:
        logger.info("Clearing weather data (excluding precipitation)...")
        weather_params = ['avgt', 'tmax', 'tmin', 'ugrd', 'vgrd', 'wind', 'avgt']
        for param in weather_params:
            weather_dir = os.path.join(save_data_loc, param)
            safe_clear_folder(weather_dir, filetypes=('','.tif','.csv','.txt'))

    if clear_uploads:
        uploads_dir = os.path.join(save_data_loc, 'uploads')
        logger.info("Clearing uploaded files...")
        safe_clear_folder(uploads_dir, filetypes=('','.zip','.json','.csv','.txt'))
