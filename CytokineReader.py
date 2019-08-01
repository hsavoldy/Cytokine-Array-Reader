'''
Created on Jun 30, 2019

@author: haels
'''

from PIL import Image
import numpy
import matplotlib.pyplot as plt
import math
import sys
import xlsxwriter 
from win32com.client import Dispatch

space_between_circles = 18;
space_between_clusters = 28;
smaller_space = 19
side_length = 9
#0.03785 
circle_space_horizontal_ratio = .0395 #.0383 #.0398 #.0395 
cluster_space_horizontal_ratio = .0624 #.065 #.06 #.0624  
vertical_space_horizontal_ratio = .063 #0.066 #.07 #.063

circle_space_vertical_ratio = .098 #0.0953 #.098 #.098 
cluster_space_vertical_ratio = .154 #0.17 #.147 #.154 
vertical_space_vertical_ratio = .158 #0.161 #.13 #.158

def find_horizontal_array_length(array):
    '''Use circle centers to find horizontal (wider) length of sample.
    Inputs:
        array: the cytokine image in array format    
    Returns:
        length: the horizontal array length
    '''
    [upper_left_center, upper_right_center, lower_left_center] = find_alignment_circle_centers(array)
    x = upper_right_center[1]-upper_left_center[1]
    y = upper_right_center[0]-upper_left_center[0]
    return math.sqrt(x**2+y**2)
    
def find_vertical_array_length(array):
    '''Use circle centers to find vertical (shorter) length of sample.
    Inputs:
        array: the cytokine image in array format
    Returns:
        length: the vertical array length
    '''
    [upper_left_center, upper_right_center, lower_left_center] = find_alignment_circle_centers(array)
    x = lower_left_center[1]-upper_left_center[1]
    y = lower_left_center[0]-upper_left_center[0]
    return math.sqrt(x**2+y**2)
    
def find_cluster_space(array):
    '''Use set ratios (based on both vertical and horizontal measurements) to 
    determine the distance, in pixels, between horizontal dot clusters.
    Inputs:
        array: the cytokine image in array format
    Returns:
        space length: the space between horizontal clusters
    '''
    x = find_horizontal_array_length(array)
    y = find_vertical_array_length(array)
    cluster_space_horizontal = cluster_space_horizontal_ratio*x
    cluster_space_vertical = cluster_space_vertical_ratio*y
    return math.floor(numpy.average([cluster_space_horizontal, cluster_space_vertical]))

def find_vertical_space(array):
    '''Use set ratios (based on both vertical and horizontal measurements) to 
    determine the distance, in pixels, between vertical dot clusters.
    Inputs:
        array: the cytokine image in array format
    Returns:
        space length: the space between vertical clusters
    '''
    x = find_horizontal_array_length(array)
    y = find_vertical_array_length(array)
    vertical_space_horizontal = vertical_space_horizontal_ratio*x
    vertical_space_vertical = vertical_space_vertical_ratio*y
    return math.floor(numpy.average([vertical_space_horizontal, vertical_space_vertical]))

def find_circle_space(array):
    '''Use set ratios (based on both vertical and horizontal measurements) to 
    determine the distance, in pixels, between individual dots- rounded for 
    accuracy.
    Inputs:
        array: the cytokine image in array format
    Returns:
        space length: the space between individual dots
    '''
    x = find_horizontal_array_length(array)
    y = find_vertical_array_length(array)
    circle_space_horizontal = circle_space_horizontal_ratio*x
    circle_space_vertical = circle_space_vertical_ratio*y
    
    if(numpy.average([circle_space_horizontal, circle_space_vertical])%1 <= .2):
        circle_space_low = math.floor(numpy.average([circle_space_vertical, circle_space_horizontal]))
        circle_space_high = math.floor(numpy.average([circle_space_vertical, circle_space_horizontal]))
    elif(numpy.average([circle_space_horizontal, circle_space_vertical])%1 >= .8):
        circle_space_low = math.ceil(numpy.average([circle_space_vertical, circle_space_horizontal]))
        circle_space_high = math.ceil(numpy.average([circle_space_vertical, circle_space_horizontal]))
    else:
        circle_space_low = math.floor(numpy.average([circle_space_vertical, circle_space_horizontal]))
        circle_space_high = math.ceil(numpy.average([circle_space_vertical, circle_space_horizontal]))
    
    return [circle_space_low, circle_space_high]


def find_upper_left_circle(array):
    '''Find upper left reference dot by minimizing the area of the square made 
    by a pixel of a threshold intensity and the upper left corner.
    Inputs:
        array: the cytokine image in array format
    Returns:
        pixel: upper left reference pixel
    '''
    threshold_intensity = 95 
    pixel_values = numpy.where(array<threshold_intensity)
    min_square = sys.maxsize
    x=0
    y=0
    
    for i in range(len(pixel_values[0])):
        if pixel_values[0][i]*pixel_values[1][i] < min_square:
            min_square = pixel_values[0][i]*pixel_values[1][i]
            x = pixel_values[0][i]
            y = pixel_values[1][i]
             
    return [x,y]

def find_lower_left_circle(array):
    '''Find lower left reference dot by minimizing the area of the square made 
    by a pixel of a threshold intensity and the lower left corner.
    Inputs:
        array: the cytokine image in array format
    Returns:
        pixel: lower left reference pixel
    '''
    threshold_intensity = 95 
    pixel_values = numpy.where(array<threshold_intensity)
    min_square = sys.maxsize
    x=0
    y=0
    
    for i in range(len(pixel_values[0])):
        if (array.shape[0]-pixel_values[0][i])*pixel_values[1][i] <min_square:
            min_square = (array.shape[0]-pixel_values[0][i])*(array.shape[1]-pixel_values[1][i])
            x = pixel_values[0][i]
            y = pixel_values[1][i]
            
    return [x,y]

def find_upper_right_circle(array):
    '''Find upper right reference dot by minimizing the area of the square made 
    by a pixel of a threshold intensity and the upper right corner.
    Inputs:
        array: the cytokine image in array format
    Returns:
        pixel: upper right reference pixel
    '''
    threshold_intensity = 95 
    pixel_values = numpy.where(array<threshold_intensity)
    min_square = sys.maxsize
    x=0
    y=0
    
    for i in range(len(pixel_values[0])):
        if pixel_values[0][i]*(array.shape[1]-pixel_values[1][i]) <min_square:
            min_square = pixel_values[0][i]*(array.shape[1]-pixel_values[1][i])
            x = pixel_values[0][i]
            y = pixel_values[1][i]
            
    return [x,y]
    

def find_circle_center_upper_left(array, starting_coordinates):
    '''Find the center of the upper left circle by rastering from left to right,
    finding the widest and darkest section of the dot, then finding the 
    middle of it.
    Inputs:
        array: the cytokine image in array format
        starting_coordinates: the upper left pixel of the dot
    Returns:
        center: the coordinates of the upper left circle center
    '''
    threshold_intensity = 160
    new_coordinates = starting_coordinates
    width = 0
    longest = 0
    averages = []
    row = []
    column_number = 0
    row_number = 0
    global upper_left_diameter
    
    while True: 
        
        while(array[starting_coordinates[0]+row_number, column_number+1]>threshold_intensity):
            column_number+=1
        while(array[starting_coordinates[0]+row_number, column_number+1]<=threshold_intensity):
            column_number += 1
            width +=1
            row.append([starting_coordinates[0]+row_number, column_number])
            averages.append(array[starting_coordinates[0]+row_number, column_number])
        if(len(row)>longest):
            longest = len(row)
            new_coordinates = row[int(round(width/2))]
            avg = numpy.mean(averages)
        elif len(row)==longest:
            if(numpy.mean(averages)<=avg):
                new_coordinates = row[int(round(width/2))]
                avg = numpy.mean(averages)
        else:
            break
        
        upper_left_diameter = longest
        column_number=0
        width=0
        row_number+=1
        row.clear()
        averages.clear()
    
    return new_coordinates    

def find_circle_center_lower_left(array, starting_coordinates):
    '''Find the center of the lower left circle by rastering from left to right,
    finding the widest and darkest section of the dot, then finding the 
    middle of it.
    Inputs:
        array: the cytokine image in array format
        starting_coordinates: the lower left pixel of the dot
    Returns:
        center: the coordinates of the lower left circle center
    '''
    threshold_intensity = 160
    new_coordinates = starting_coordinates
    width = 0
    longest = 0
    averages = []
    row = []
    column_number=0
    row_number=0
    plt.imshow(array)
    global lower_left_diameter 
    
    while True: 
        
        while(array[starting_coordinates[0]-row_number, column_number+1]>threshold_intensity):
            column_number+=1
        while(array[starting_coordinates[0]-row_number, column_number+1]<=threshold_intensity):
            column_number += 1
            width +=1
            row.append([starting_coordinates[0]-row_number, column_number])
            averages.append(array[starting_coordinates[0]-row_number, column_number])
        if(len(row)>longest):
            longest = len(row)
            new_coordinates = row[int(round(width/2))]
            avg = numpy.mean(averages)
        elif len(row)==longest:
            if(numpy.mean(averages)<=avg):
                new_coordinates = row[int(round(width/2))]
                avg = numpy.mean(averages)
        else:
            break
        
        lower_left_diameter = longest
        column_number=0
        width=0
        row_number+=1
        row.clear()
        averages.clear()
        
    return new_coordinates  

def find_circle_center_upper_right(array, starting_coordinates):
    '''Find the center of the lower left circle by rastering from right to left,
    finding the widest and darkest section of the dot, then finding the 
    middle of it.
    Inputs:
        array: the cytokine image in array format
        starting_coordinates: the upper right pixel of the dot
    Returns:
        center: the coordinates of the upper right circle center
    '''
    threshold_intensity = 160
    width = 0
    new_coordinates = starting_coordinates
    longest = 0
    averages = []
    row = []
    column_number=0
    row_number=0
    end = array.shape[1]
    global upper_right_diameter 
    
    while True: 
        while(array[starting_coordinates[0]+row_number, end-(column_number+1)]>threshold_intensity):
            column_number+=1
        while(array[starting_coordinates[0]+row_number, end-(column_number+1)]<=threshold_intensity):
            column_number += 1
            width +=1
            row.append([starting_coordinates[0]+row_number, end-column_number])
            averages.append(array[starting_coordinates[0]+row_number, end-column_number])
        if(len(row)>longest):
            longest = len(row)
            new_coordinates = row[int(round(width/2))]
            avg = numpy.mean(averages)
        elif len(row)==longest:
            if(numpy.mean(averages)<=avg):
                new_coordinates = row[int(round(width/2))]
                avg = numpy.mean(averages)
        else:
            break
        
        upper_right_diameter = longest
        column_number=0
        width=0
        row_number+=1
        row.clear()
        averages.clear()
        
    return new_coordinates 

def find_alignment_circles(array):
    '''Find each reference dot.
    Inputs:
        array: the cytokine image in array format
    Returns:
        coordinates: coordinates of each reference dot- upper_left, upper_right,
            lower_right
    '''
    upper_left = find_upper_left_circle(array)
    upper_right = find_upper_right_circle(array)
    lower_right = find_lower_left_circle(array)
    return [upper_left, upper_right, lower_right]
     
def find_alignment_circle_centers(array):
    '''Find the center of each reference dot.
    Inputs:
        array: the cytokine image in array format
    Returns:
        coordinates: coordinates of the center of each reference dot- 
            upper_left, upper_right, lower_right
    '''
    [upper_left, upper_right, lower_right] = find_alignment_circles(array)
    upper_left_center = find_circle_center_upper_left(array, upper_left)
    upper_right_center = find_circle_center_upper_right(array, upper_right)
    lower_left_center = find_circle_center_lower_left(array, lower_right)
    return [upper_left_center, upper_right_center, lower_left_center]

def find_angle(array):
    '''Find the angle of rotation of the array using the differences between 
    the upper right and upper left dot centers.
    Inputs:
        array: the cytokine image in array format
    Returns:
        angle: the offset angle of the array
    '''
    [upper_left_center, upper_right_center, lower_left_center] = find_alignment_circle_centers(array)
    x = upper_right_center[1] - upper_left_center[1] 
    y = upper_right_center[0]-upper_left_center[0]
    return abs(math.atan(y/x))
    
def rotate_image(im):
    '''Rotate the image to make it straight.
    Inputs:
        im: the cytokine image
        array: the cytokine image in array format
    Returns:
        rotated_image: the straightened image
    '''
    array = numpy.asarray(im)
    [upper_left_center, upper_right_center, lower_left_center] = find_alignment_circle_centers(array)
    
    if (upper_left_center[0] > upper_right_center[0]):
        angle = - math.degrees(find_angle(array))
    else:
        angle = math.degrees(find_angle(array))
        
    return im.rotate(angle, expand=False, center = [upper_left_center[1],upper_left_center[0]])
 
def find_size_of_dot():
    '''Find the best measurement diameter by finding the closest odd number to 
    the average reference dot diameter.
    Returns:
        diameter: the average calculated measurement dot diameter
    '''
    closest_odd_number = 0
    average_diameter = math.floor((upper_right_diameter+lower_left_diameter+upper_left_diameter)/3)
    
    if(upper_right_diameter%2==0):
        closest_odd_number = average_diameter
    else:
        closest_odd_number = average_diameter+1
        
    return math.floor(closest_odd_number*.47)

def find_intensity_of_dot(coordinates, array, testArray):
    '''Find the average intensity of the measurement area after the area is 
    adjusted to minimize intensity.
    Inputs:
        coordinates: the coordinates of the guessed center of the dot
        array: the cytokine image in array format
    Returns:
        average_intensity: the average intensity of the adjusted measurement
        area
    '''
    dot_size = find_size_of_dot()
    dot_radius = dot_size//2
    area = []
    dot_pixel_locations = []
    
    for row_number in range(dot_size):
        for column_number in range(dot_size):
            area.append(array[coordinates[0]-dot_radius+row_number, coordinates[1]-dot_radius+column_number])
            dot_pixel_locations.append([coordinates[0]-dot_radius+row_number, coordinates[1]-dot_radius+column_number])
    average_intensity = center_dot(array, dot_pixel_locations, area, testArray)
    
    return average_intensity

def center_dot(array, dot_pixel_locations, area, testArray):
    '''Find the average intensity of the centered measurement area that will 
    minimize intensity
    Inputs:
        array: the cytokine image in array format
        dot_pixel_locations: list of pixel coordinates in measurement area
    Returns:
        average_intensity: the average intensity of the adjusted measurement
        area
    '''
    counter = 0
    dot_size = find_size_of_dot()
    top_left = area[0]
    top_right = area[dot_size-1]
    bottom_left = area[-dot_size]
    bottom_right = area[-1]
    corners = [top_left, top_right, bottom_left, bottom_right]
    
    while((needs_adjustment(corners)) & (counter<50)):
        [corners, area, dot_pixel_locations] = adjust_dot(array, dot_pixel_locations, testArray, corners)
        counter +=1
    for pixel in dot_pixel_locations:
        testArray[pixel[0], pixel[1]] = 255
    return numpy.mean(area)
        
def adjust_dot(array, dot_pixel_locations, testArray, corners):  
    '''Find the lightest corner and move the dot in the opposite direction; 
    return the corners, area, and pixel locations in a list
    Inputs:
        array: the cytokine image in array format
        dot_pixel_locations: list of pixel coordinates in measurement area
        corners: list of coordinates of the measurement area corners
    Returns:
        [corners, area, new_dot_pixel_locations]: list of coordinates of corners, 
        list of pixel intensities, and list of pixel locations of adjusted dot
    '''
    lightest_corner = corners.index(max(corners))
    if(lightest_corner == 0):
        return move_dot("top left", array, dot_pixel_locations, testArray, corners)
    elif(lightest_corner==1):
        return move_dot("top right", array, dot_pixel_locations, testArray, corners)
    elif(lightest_corner==2):
        return move_dot("bottom left", array, dot_pixel_locations, testArray, corners)
    else:
        return move_dot("bottom right", array, dot_pixel_locations, testArray, corners)
    
def move_dot(location, array, dot_pixel_locations, testArray, corners): 
    '''Move the dot in the opposite direction of the given corner
    Inputs:
        array: the cytokine image in array format
        dot_pixel_locations: list of pixel coordinates in measurement area
        corners: list of coordinates of the measurement area corners
    Returns:
        [corners, area, new_dot_pixel_locations]: list of coordinates of corners, 
        list of pixel intensities, and list of pixel locations of adjusted dot
    '''
    dot_size = find_size_of_dot()
    area = []
    new_dot_pixel_locations = []
    i = 0
    j = 0
    if(location == "top left"):
        i=1
        j=1
    elif(location == "top right"):
        i=1
        j=-1
    elif(location == "bottom left"):
        i=-1
        j=1
    else:
        i=-1
        j=-1
    
    for pixel in dot_pixel_locations:
        pixel = [pixel[0]+i, pixel[1]+j]
        area.append(array[pixel[0], pixel[1]])
        new_dot_pixel_locations.append(pixel)

    top_left = area[0]
    top_right = area[dot_size-1]
    bottom_left = area[-dot_size]
    bottom_right = area[-1]
    corners = [top_left, top_right, bottom_left, bottom_right]
    
    return[corners, area, new_dot_pixel_locations]
    
def needs_adjustment(corners): 
    '''Determines whether the greatest difference in corner intensity warrants
    an adjustment
    Inputs:
        corners: list of coordinates of the measurement area corners
    Returns:
        needs_adjustment: true if needs adjustment
    '''
    epsilon = 10
    max_difference = max(corners)-min(corners)
    if(max_difference>epsilon):
        return True
    else:
        return False

def analyze_dots(im, array):
    '''Searches for dots on the array in guessed locations; stores average 
    intensity of each dot in a list.
    Inputs:
        array: the cytokine image in array format
        im: the cytokine image
    Returns:
        data: list of each dot's average intensity
    '''
    [upper_left_center, upper_right_center, lower_left_center] = find_alignment_circle_centers(array)
    image = rotate_image(im)
    image_array = numpy.asarray(image)
    image_marked = image_array.copy()
    data = []
    horizontal_distance = 0
    circle_spaces = find_circle_space(array)
    cluster_space = find_cluster_space(array)
    vertical_space = find_vertical_space(array)
    vertical_distance = 0

    for j in range(10):
        for i in range(24):
            data.append(find_intensity_of_dot([upper_left_center[0]+vertical_distance,upper_left_center[1]+horizontal_distance], image_array, image_marked)) 
            next = i+1
            if (next%4==0) & (next!=0):
                horizontal_distance += cluster_space
            elif (next%2==0):
                horizontal_distance += circle_spaces[0]
            else:
                horizontal_distance += circle_spaces[1]
        
        horizontal_distance = 0
        if((j ==2) | (j ==6)):
            vertical_distance += vertical_space
        else:
            vertical_distance += circle_spaces[0]
              
    #plt.imshow(image_marked)
    #plt.show()
    return data

def filter_cytokines(array_cytokine_intensities):
    '''Find which cytokines appear in any array then filter out each array's 
    list of cytokines to only include those.
    Inputs:
        array_cytokine_intensities: list of dictionaries of each cytokine array's 
        cytokine name and corresponding intensity value 
    Returns:
        array_cytokine_intensities: Input filtered to only include wanted 
        cytokines
    '''
    threshold = 25
    valid_keys= []
    final = {}
    
    for array in array_cytokine_intensities:
        final = {k: v for k, v in array.items() if v>threshold}
        valid_keys = valid_keys + list(final.keys())
        
    valid_keys = [x for x in valid_keys if x != "Reference"]
    
    for i in range(len(array_cytokine_intensities)):
        array_cytokine_intensities[i] = {k: v for k, v in array_cytokine_intensities[i].items() if (k in valid_keys)}
    return array_cytokine_intensities
            
    
def assign_cytokines(data):
    '''Assign cytokines to intensity values.
    Inputs:
        data: list of intensity values
    Returns:
        cytokine_dict: dictionary with cytokine names as keys corresponding
        to intensity values
    '''
    cytokine_dict = {}
    cytokines = mouse_cytokines()
    data_index=0
    cytokine_index = 0
    while (data_index < len(data)-1):
        cytokine_dict[cytokines[cytokine_index]] = round(255-(data[data_index]+data[data_index+1])/2, 2)
        data_index+=2
        cytokine_index += 1
    return cytokine_dict

def mouse_cytokines():
    '''Gives list of cytokines that correlate to data points.
    Returns:
        cytokines: list of cytokines in order of array
    '''
    cytokines = ["Reference", "Adiponectin/Acrp30", "Amphiregulin", "Angiopoietin-1", 
                "Angiopoietin-2", "Angiopoietin-like 3", "BAFF/BLyS/TNFSF13B",
                "C1qR1/CD93", "CCL2/JE/MCP-1", "CCL3/CCL4/MIP-1α/β", 
                "CCL5/RANTES", "Reference", "Space", "CCL6/C10", "CCL11/Eotaxin", "CCL12/MCP-5",
                "CCL17/TARC", "CCL19/MIP-3β", "CCL20/MIP-3α", "CCL21/6Ckine",
                "CCL22/MDC", "CD14", "CD40/TNFRSF5", "Space", "Space", "CD160", "Chemerin",
                "Chitinase 3-like 1", "Coagulation Factor III/Tissue Factor",
                "Complement Component C5/C5a", "Complement Factor D",
                "C-Reactive Protein/CRP", "CX3CL1/Fractalkine", "CXCL1/KC",
                "CXCL2/MIP-2", "Space", "CXCL9/MIG", "CXCL10/IP-10", "CXCL11/I-TAC",
                "CXCL13/BLC/BCA-1", "CXCL16", "Cystatin C", "DKK-1", 
                "DPPIV/CD26", "EGF", "Endoglin/CD105", "Endostatin", 
                "Fetuin A/AHSG", "FGF acidic", "FGF-21", "Flt-3 Ligand", 
                "Gas 6", "G-CSF", "GDF-15", "GM-CSF", "HGF", "ICAM-1/CD54",
                "IFN-γ", "IGFBP-1", "IGFBP-2", "IGFBP-3", "IGFBP-5", "IGFBP-6",
                "IL-1α/IL-1F1", "IL-1β/IL-1F3", "IL-1ra/IL-1F3", "IL-2", "IL-3",
                "IL-4", "IL-5", "IL-6", "IL-7", "IL-10", "IL-11", "IL-12 p40",
                "IL-13", "IL-15", "IL-17A", "IL-22", "IL-23", "IL-27 p28",
                "IL-28A/B", "IL-33", "LDL R", "Leptin", "LIF", "Lipocalin-2/NGAL",
                "LIX", "M-CSF", "MMP-2", "MMP-3", "MMP-9", "Myeloperoxidase",
                "Osteopontin (OPN)", "Osteoprotegerin/TNFRSF11B",
                "PD-ECGF/Thymidine phosphorylase", "PDGF-BB", "Pentraxin 2/SAP",
                "Pentraxin 3/TSG-14", "Periostin/OSF-2", "Pref-1/DLK-1/FA1",
                "Proliferin", "Proprotein Convertase 9/PCSK9", "RAGE", "RBP4",
                "Reg3G", "Resistin", "Space", "Reference", "E-Selectin/CD62E", "P-Selectin/CD62P",
                "Serpin E1/PAI-1", "Serpin F1/PEDF", "Thrombopoietin",
                "TIM-1/KIM-1/HAVCR", "TNF-α", "VCAM-1/CD106", "VEGF", 
                "WISP-1/CCN4", "Space"]
    return cytokines
    
def create_excel_graph(cytokine_dicts, treatments):
    # Workbook() takes one, non-optional, argument   
    # which is the filename that we want to create. 
    workbook = xlsxwriter.Workbook('cytokineChart.xlsx') 
  
    # The workbook object is then used to add new   
    # worksheet via the add_worksheet() method.  
    worksheet = workbook.add_worksheet() 
  
    # Create a new Format object to formats cells 
    # in worksheets using add_format() method . 
  
    # here we create bold format object . 
    bold = workbook.add_format({'bold': 1}) 
  
    # create a data list . 
    headings = [''] + treatments
  
    data = [list(cytokine_dicts[0].keys())]
    for i in range(len(cytokine_dicts)):
        data.append(list(cytokine_dicts[i].values()))
  
    # Write a row of data starting from 'A1' 
    # with bold format. 
    worksheet.write_row('A1', headings, bold) 
  
    # Write a column of data starting from 
    # A2, B2, C2 respectively. 
    rows = ['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2', 'J2', 'K2', "L2"]
    for i in range(len(cytokine_dicts)+1):
        worksheet.write_column(rows[i], data[i]) 
  
    # Create a chart object that can be added 
    # to a worksheet using add_chart() method. 
  
    # here we create a pie chart object . 
    chart1 = workbook.add_chart({'type': 'column'}) 
  
    # Add a data series to a chart 
    # using add_series method. 
    # Configure the first series. 
    #[sheetname, first_row, first_col, last_row, last_col]. 
    for i in range(len(treatments)):
        chart1.add_series({ 
            'name':       ['Sheet1', 0, i+1], 
            'categories': ['Sheet1', 1, 0, len(cytokine_dicts[0]),0], 
            'values':     ['Sheet1', 1, i+1, len(cytokine_dicts[0]), i+1], 
            }) 
  
    # Add a chart title  
    chart1.set_title({'cytokine': 'intensities'}) 
  
    # Set an Excel chart style. Colors with white outline and shadow. 
    chart1.set_style(2) 
    
    chart1.set_size({'x_scale': 2.5, 'y_scale':1.5})
  
    # Insert the chart into the worksheet(with an offset). 
    # the top-left corner of a chart is anchored to cell C2.  
    worksheet.insert_chart(rows[i+3], chart1) 
  
    # Finally, close the Excel file   
    # via the close() method.   
    workbook.close()  
    
def create_graph(cytokine_dicts, names):
    
    #cytokine_dicts = [{'amphiregulin': 45.67, 'angiopiotin':76.98},{'amphiregulin': 189.87, 'angiopiotin':10.87},{'amphiregulin': 84.67, 'angiopiotin':43.90},{'amphiregulin': 12.32, 'angiopiotin':1.92}]
    
#     n_groups = 4
#     means_frank = (90, 55, 40, 65)
#     means_guido = (85, 62, 54, 20)
# 
#     # create plot
#     fig, ax = plt.subplots()
#     index = numpy.arange(n_groups)
#     bar_width = 0.35
#     opacity = 0.8
# 
#     rects1 = plt.bar(index, means_frank, bar_width,
#                      alpha=opacity,
#                      color='b',
#                      label='Frank')
#     
    n_groups = len(cytokine_dicts[0])
    n_cytokines = len(cytokine_dicts)
    data = []
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    while(len(colors) < n_groups):
        colors = colors+colors
    # create plot
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)
    index = numpy.arange(n_groups)
    
    #Create graph dimensions
    bar_width = 0.015
    same_cytokine_distance = .02
    space = .06*n_groups
    
    opacity = 0.8
    plt.xlabel('Cytokine')
    plt.ylabel('Intensity')
    plt.title('Scores by person')
    
    #create cytokine labels
    start_amount = n_cytokines*(same_cytokine_distance)/2
    ticks = [.03+space*i for i in range(n_groups)]
    plt.xticks(ticks, cytokine_dicts[0].keys())
    
    plt.legend()
# 
#     plt.bar((same_cytokine_distance*0, same_cytokine_distance*0+space), tuple(cytokine_dicts[0].values()), bar_width, 
#             alpha=opacity,color=colors[0], label=names[0])
#     plt.show()
#     plt.bar((same_cytokine_distance*1, same_cytokine_distance*1+space), tuple(cytokine_dicts[1].values()), bar_width, 
#             alpha=opacity,color=colors[1], label=names[1])
#     plt.show()
#     plt.bar((same_cytokine_distance*2, same_cytokine_distance*2+space), tuple(cytokine_dicts[2].values()), bar_width, 
#             alpha=opacity,color=colors[2], label=names[2])
#     plt.show()
#     plt.bar((same_cytokine_distance*3, same_cytokine_distance*3+space), tuple(cytokine_dicts[3].values()), bar_width, 
#             alpha=opacity,color=colors[3], label=names[3])
#     plt.show()
    for i in range(len(names)):
        plt.bar((same_cytokine_distance*i, same_cytokine_distance*i+space), tuple(cytokine_dicts[i].values()), bar_width, 
                alpha=opacity,color=colors[i], label=names[i])

    plt.tight_layout()
    plt.show()
    print('hey')
    
#     array_labels = names
#     cytokine_labels = list(cytokine_dicts[0].keys())
#     data_points = []
#     
#     fig = plt.figure()
#     fig, ax = plt.subplots()
#     
#     x = numpy.arange(len(array_labels))  # the label locations
#     width = 0.35  # the width of the bars
#     rects = []
#     for i in range(len(array_labels)):
#         data_points.append(list(cytokine_dicts[i].values()))
#         rects.append(ax.bar(data_points[i], width, label=array_labels[i]))
# 
#     # Add some text for labels, title and custom x-axis tick labels, etc.
#     ax.set_ylabel('Intensity')
#     ax.set_title('Intensity of Cytokines')
#     ax.set_xticks(x)
#     ax.set_xticklabels(cytokine_labels)
#     ax.legend()
# 
# 
#     def autolabel(rects):
#         """Attach a text label above each bar in *rects*, displaying its height."""
#         for rect in rects:
#             height = rect.get_height()
#             ax.annotate('{}'.format(height),
#                         xy=(rect.get_x() + rect.get_width() / 2, height),
#                         xytext=(0, 3),  # 3 points vertical offset
#                         textcoords="offset points",
#                         ha='center', va='bottom')
# 
# 
#     for rect in rects:
#         autolabel(rect)
#     
#     fig.tight_layout()
# 
#     plt.show()
#     print('done')
    
def main():
    first_image = Image.open('CytokineSample1.jpg').convert('L') 
    first_array = numpy.asarray(first_image)
    first_data = analyze_dots(first_image, first_array)
    
    second_image = Image.open('CytokineSample2.jpg').convert('L') 
    second_array = numpy.asarray(second_image)
    second_data = analyze_dots(second_image, second_array)
    
    third_image = Image.open('CytokineSample3.jpg').convert('L') 
    third_array = numpy.asarray(third_image)
    third_data = analyze_dots(third_image, third_array)
    
    fourth_image = Image.open('CytokineSample4.jpg').convert('L') 
    fourth_array = numpy.asarray(fourth_image)
    fourth_data = analyze_dots(fourth_image, fourth_array)
    
    names = ['DMSO 3TC', 'DMSO', 'Rev 3TC', 'Rev']
    
    data = filter_cytokines([assign_cytokines(first_data), assign_cytokines(second_data), assign_cytokines(third_data), assign_cytokines(fourth_data)])
    create_excel_graph(data, names)
if __name__ == '__main__':
    main()