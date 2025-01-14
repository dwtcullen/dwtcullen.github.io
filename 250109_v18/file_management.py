"""
File management toolbox | David Cullen | Jones Lab 2024
"""

import pandas as pd
import numpy as np
import os
import shutil

def set_directory(path, user):
    os.chdir(f'/scratch/user/{user}/{path}')
    
def folders(output_masks = False,
            has_KTR = False,
            output_image_stacks = False) -> None:
    ''' Create required folders, removing existing ones if necessary. '''
    if output_masks:
        folder = 'Output_Images/Nuclear_Masks/'
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)
        if has_KTR:
            folder = 'Output_Images/KTR_Masks/'
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder)

    if output_image_stacks:
        folder = 'Output_Images/Stacks/'
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)

import re
def extract_file_info(
    file_name: str, 
    microscope_name: str, 
    ) -> dict[str, str]:
    """
    Extracts the file information from the file name using regular expression.

    Parameters:
        file_name (str): The file name to extract information from.
        microscope_name (str) Either 'Vader', 'Maul' or 'Deltavision'.

    Returns:
        Dict[str, str]: A dictionary with keys as the extracted file information and values 
        as the corresponding information. The regex used is aimed to extract information from 
        Incell Analysier 6500 tifs or live cell imaging data from the Deltavision Ultra.
    """
    microscope_regex = {
        'Vader': r'(?P<Row>[a-zA-Z])_(?P<Column>[0-9]+)_fld_(?P<Field>[0-9]+)_wv_(?P<WV>[0-9]+)_(?P<Filter>[a-zA-Z]+)',
        'Maul': r'(?P<Row>[a-zA-Z]) - (?P<Column>[0-9]+)\(fld (?P<Field>[0-9]+) wv (?P<WV>[0-9]+) - (?P<Filter>[a-zA-Z]+)\)',
        'Deltavision': r'expt_(?P<Row>[A-Za-z])(?P<Column>[0-9]+)_R3D_w(?P<WV>[0-9]+)_p(?P<Field>[0-9]+)\.tif',
        'Operetta': r'r(?P<Row>[0-9]+)c(?P<Column>[0-9]+)f(?P<Field>[0-9]+)p(?P<Plane>[0-9]+)-ch(?P<Channel>[0-9]+)t(?P<Time>[0-9]+)'
    }
    pattern = microscope_regex[microscope_name]
    match = re.match(pattern, file_name)
    
    if match is None:
        return {}
    
    components = match.groupdict()
    if microscope_name == 'Operetta':
        numeric_row = int(components['Row'])
        components['Row'] = chr(64 + numeric_row)
    
    return components


def get_metadata(path: str, 
                 microscope_name: str,
                 groups: list[str]
                ) -> pd.DataFrame:
    """
    Retrieve metadata information for files in a directory and return the grouped metadata.

    Parameters:
    path (str): Path to the directory where the files are located.
    microscope_name (str): Name of the microscope I'm using. 'Vader', 'Maul', or 'Deltavision'.
    groups (List[str]): List of column names to use for grouping the metadata.

    Returns:
    pd.DataFrame: A dataframe containing the grouped metadata information.
    """

    file_list = [entry.name for entry in os.scandir(path) if entry.is_file()]
    metadata = [{'String': file_name, **extract_file_info(file_name, microscope_name=microscope_name)} for file_name in file_list]
    metadata = pd.DataFrame(metadata)
    if microscope_name == 'Vader':
        metadata['Column'] = metadata['Column'].astype(int)
        metadata['Column'] = metadata['Column'].apply(lambda x: f"{x:02d}")
        metadata['Field'] = metadata['Field'].astype(int)
        metadata['Field'] = metadata['Field'].apply(lambda x: f"{x:02d}")
    
    metadata = metadata.sort_values(by='WV')
    return metadata.groupby(groups)

      
def get_DAPI_string(df):
    """
    Function to retrieve the first matching string from the 'String' column
    where 'WV' equals '405' in a given DataFrame.
    
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing the columns 'WV' and 'String'.
    
    Returns:
    str: The first matching string if found, otherwise None.
    """
    # Filter the DataFrame
    filtered = df.loc[df['WV'] == '405', 'String']
    
    # Check if the filtered DataFrame is not empty and return the first entry
    if not filtered.empty:
        return filtered.iloc[0]
    else:
        return None

import pandas as pd    
def write_df(results: list[dict], 
             image: np.ndarray, 
             output_df_name: str,
             green_ratio: list = None, 
             orange_ratio: list = None, 
             green_ring_intensity: list = None,
             orange_ring_intensity: list = None,
             has_KTR: bool = False,
             measure_pearson: bool = False
            ) -> None:
    '''
    Write dataframe, merge it with the plate map and export as a CSV.

    Parameters:
    - results: A list of dictionaries containing measured properties of nuclei.
    - image: numpy ndarray representing the image.
    - green_ratio, orange_ratio, green_ring_intensity, orange_ring_intensity: Optional additional metrics.
    - has_KTR: Flag indicating if KTR measurements should be included.

    Returns:
    - None
    '''
    df = pd.DataFrame(results)
    
    if has_KTR:
        df['CDK2_Activity'] = green_ratio
        df['CDK46_Activity'] = orange_ratio
        df['Green_Ring_Mean'] = green_ring_intensity
        df['Red_Ring_Mean'] = orange_ring_intensity

    df.to_csv(f'{output_df_name}.csv', index=False)
    print(f"Image processing completed.{ ' ' * 100 }", end="\r")
   

import re
def parse_image_name(image_name, pattern):
    """
    Parse the image name using the provided regex pattern.
    
    Args:
    - image_name (str): The image file name to parse.
    - pattern (str): The regex pattern to match against.
    
    Returns:
    - bool: True if the pattern matches, False otherwise.
    """
    match = re.match(pattern, image_name)
    if match:
        components = match.groupdict()
        print(f"File name: {image_name}")
        print(f"Pattern: {pattern}")
        print(f"Row: {components['Row']}")
        print(f"Column: {components['Column']}")
        print(f"Field: {components['Field']}")
        print(f"WV: {components['WV']}")
        print(f"Filter: {components['Filter']}")
        return True
    else:
        print(f"File name: {image_name} - Pattern did not match")
        return False
