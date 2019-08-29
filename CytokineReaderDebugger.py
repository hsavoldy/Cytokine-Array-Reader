'''
Created on Aug 28, 2019

@author: haels
'''
from PIL import Image
import numpy
import matplotlib.pyplot as plt
import math
import sys
import glob, os

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
        
    return math.floor(closest_odd_number*.45) #.47

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
        [corners, area, dot_pixel_locations] = adjust_dot(array, dot_pixel_locations, corners)
        counter +=1
    
    for pixel in dot_pixel_locations:
        if(testArray[pixel[0], pixel[1]] == 255): #to avoid accidental overlap
            return 255
        else:
            testArray[pixel[0], pixel[1]] = 255
        
    return numpy.mean(area)
        
def adjust_dot(array, dot_pixel_locations, corners):  
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
        return move_dot("top left", array, dot_pixel_locations)
    elif(lightest_corner==1):
        return move_dot("top right", array, dot_pixel_locations)
    elif(lightest_corner==2):
        return move_dot("bottom left", array, dot_pixel_locations)
    else:
        return move_dot("bottom right", array, dot_pixel_locations)
    
def move_dot(location, array, dot_pixel_locations): 
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
              
    plt.imshow(image_marked)
    plt.show()
    return data

def main():

    prompt = '>'

    print("Hello- welcome to Cytokine Reader Debugger!")
    
    files = []
    
    while(True):
        print("Please enter the directory path (eg. C:\\user\\arrays) of your array images: ")
        filename = input(prompt).replace('/', "\\")
        try:
            os.chdir(filename)
            break
        except:
            print('That is not a valid path. Please make sure you have the correct spellings.')
        
    
    # TODO handle if directory not found
    os.chdir(filename)
    for file in glob.glob("*.jpg"):
        files.append(file)
        
    while(len(files) == 0):
        print("No JPGs found. Make sure your files are present and in JPG format. \nPlease enter the location of your array images: ")
        filename = input(prompt)
    
        os.chdir(filename)
        for file in glob.glob("*.jpg"):
            files.append(file)
    
    print('Producing images...')
    
    for i in range(len(files)):
        image = Image.open(files[i]).convert('L') 
        plt.title(files[i])
        array = numpy.asarray(image)
        analyze_dots(image, array)
    
if __name__ == '__main__':
    main()