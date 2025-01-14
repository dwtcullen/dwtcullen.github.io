"""
Propidium Iodide Cell Death Assay Toolbox | David Cullen | Jones Lab 2024
"""
import numpy as np
def rgb(image: np.ndarray):
    image_height = image.shape[1]
    rgb = np.zeros(shape=(image_height, image_height, 3),dtype=np.uint16)

    rgb[:,:,0] = image[0]
    rgb[:,:,1] = image[1]
    return rgb

import cellpose
from cellpose import models
def segment_pi(image: np.ndarray,
                   progress_report: str,
                   user: str
                  ) -> np.ndarray:
    '''
    Use pre-trained model to segment nuclei in an image. 
    This is the default function, for 10x air.

    Parameters:
        image (np.ndarray): A single channel image.

    Returns:
        np.ndarray: A binary mask of the segmented nuclei.
    '''

    from skimage.segmentation import clear_border
    print(progress_report + "|Segmenting PI...                                                                                                             ", end='\r')
    #model = models.CellposeModel(model_type='nuclei')
    #model = models.CellposeModel(model_type='cellpose7')
    model = models.CellposeModel(model_type=f'/scratch/user/{user}/cellpose/models/241029_PI')
    masks = model.eval(image[-1], channels=[0,0])
    masks = clear_border(masks[0])      
    return masks

import numpy as np
import image_output
def intensity_sum(regionmask: np.ndarray, 
                  intensity_image: np.ndarray
                 ) -> float:
    '''
    Return the total intensity for all pixels inside the mask.

    Parameters
    ----------
    regionmask : np.ndarray
        A binary mask indicating which pixels belong to the region of interest.
    intensity_image : np.ndarray
        The 2D intensity image corresponding to the regionmask.

    Returns:
    float
        The sum of intensities for all pixels inside the regionmask.
    '''
    return np.sum(intensity_image[regionmask])

def intensity_mean_eroded(regionmask: np.ndarray, 
                          intensity_image: np.ndarray
                         ) -> float:
    '''
    Erode the regionmask, then find the mean intensity of each channel within.
    
    Parameters:
    - regionmask (np.ndarray): The binary mask representing the region of interest.
    - intensity_image (np.ndarray): The multi-channel intensity image.
    
    Returns:
    - float: The mean intensity of the eroded regionmask.
    '''
    eroded = erode(regionmask, erosion = 1)
    eroded_image = intensity_image[eroded]
    return np.mean(eroded_image)

def measure_dist_from_centre(centroid: tuple[int, int],
                             image_centre: float):
    '''
    Basic trig to calculate the distance from the centre of the nucleus to the centre of the image.
    
    Parameters:
    - centroid: Tuple of two elements representing the row and column in which the centre of the nucleus is located in the image.
    -image: ndarray of the current image.
    
    Returns:
    - dist_from_centre: a float representing the distance from the centre of the nucleus to the centre of the image.
    '''
    prop_row, prop_col = centroid
    rows_from_centre = np.abs(prop_row-image_centre)
    cols_from_centre = np.abs(prop_col-image_centre)
    dist_from_centre = np.sqrt(rows_from_centre**2 + cols_from_centre**2)
    return dist_from_centre
    
def get_regionprops(masks, image):
    from skimage.measure import regionprops
    regionprops_image = np.moveaxis(image, 0, -1)
    props = regionprops(masks, intensity_image=regionprops_image, extra_properties=(intensity_sum,))
    return props

def measure_properties(name: tuple[str, str, int], 
                       image: np.ndarray, 
                       nuclear_regions: list, 
                       PI_regions: list,
                       results: list,
                       progress_report: str,
                       nuclear_masks,
                       PI_masks,
                       calculate_pearson: bool = False
                      ) -> list[tuple]:
    import matplotlib.pyplot as plt
    
    '''
    Measure properties of nuclei in an image.
    
    Parameters:
    - name: Tuple of three elements representing the well, row, and column of the current image
    - image: numpy ndarray representing the image to be processed
    - regions: list of regionprops objects representing the nuclei in the image
    
    Returns:
    - results: a list of lists where each sublist represents the properties of a single nucleus
    
    The properties included in each sublist are:
    - well: well identifier, a concatenation of the row and column in the form "RowColumn"
    - field: identifier of the field within the well
    - label: identifier of the nucleus within the field
    - area: area of the nucleus, in pixels
    - circularity: measure of circularity, defined as 4 * pi * area / perimeter^2
    - Mean_DAPI: the mean log10 intensity of the nucleus in the DAPI channel
    - Mean_DAPI: the mean log10 intensity of the nucleus in the DAPI channel
    - Mean_Green (optional): the mean log10 intensity of the nucleus in the Green channel, if present
    - Mean_Orange (optional): the mean log10 intensity of the nucleus in the Orange channel, if present
    - Mean_EdU (optional): the mean log10 intensity of the nucleus in the EdU channel, if present
    - Sum_DAPI: the sum log10 intensity of the nucleus in the DAPI channel
    - Sum_Green (optional): the sum log10 intensity of the nucleus in the Green channel, if present
    - Sum_Orange (optional): the sum log10 intensity of the nucleus in the Orange channel, if present
    - Sum_EdU (optional): the sum log10 intensity of the nucleus in the EdU channel, if present
    '''
   
    print(progress_report + f"|Measuring nuclear properties...                                                                                               ", end='\r')
    for props in PI_regions:
        status = 'Dead'
        circularity = 4 * np.pi * props.area / (props.perimeter**2)
        distance_from_centre = measure_dist_from_centre(props.centroid, image.shape[1]/2)
        row, column, field = name
        result_dict = {
            'Column': column,
            'Row': row
        }
        
        result = [status,
                  column,
                  row,
                  row+column,
                  field, 
                  props.label,
                  props.area,
                  circularity,
                  props.centroid,
                  props.centroid[0],
                  props.centroid[1],
                  distance_from_centre] + [np.log10(i) for i in props.intensity_mean] + [np.log10(i) for i in props.intensity_sum]
        results.append(result)
        
    for props in nuclear_regions:
        #if np.log10(props.intensity_mean[1])>2.25:
        #    image_output.output_cropped_images(props, image, name)
        status = 'Live'
        circularity = 4 * np.pi * props.area / (props.perimeter**2)
        distance_from_centre = measure_dist_from_centre(props.centroid, image.shape[1]/2)
        row, column, field = name
        result_dict = {
            'Column': column,
            'Row': row
        }
        result = [status,
                  column,
                  row,
                  row+column,
                  field, 
                  props.label,
                  props.area,
                  circularity,
                  props.centroid,
                  props.centroid[0],
                  props.centroid[1],
                  distance_from_centre] + [np.log10(i) for i in props.intensity_mean] + [np.log10(i) for i in props.intensity_sum]
        results.append(result)
        
        #if calculate_pearson:
        #    pearson = pearson(props)
        #    result.append(pearson)
            
        