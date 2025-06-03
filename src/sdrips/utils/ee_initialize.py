import ee
import logging

def initialize_earth_engine():
    """
    Authenticate and initialize the Earth Engine API.

    If the user is not authenticated, this function will prompt for authentication.
    Once authenticated, it initializes the Earth Engine library.

    Raises:
        Exception: If initialization fails after authentication.
    """
    logger = logging.getLogger()

    try:
        ee.Initialize()
        logger.info("Earth Engine initialized successfully.")
    except Exception as init_error:
        logger.warning("Earth Engine not initialized. Attempting authentication...")
        try:
            ee.Authenticate()
            ee.Initialize()
            logger.info("Earth Engine authenticated and initialized successfully.")
        except Exception as auth_error:
            logger.error(f"Earth Engine authentication or initialization failed: {auth_error}")
            raise
