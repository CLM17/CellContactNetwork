def make_spiral_grid(width, height):
    '''
    The high content microscope images a well in spiral shape.
    This function makes a spiral-shaped grid, that we can use to find the locations of the image.
    '''
    NORTH, S, W, E = (0, -1), (0, 1), (-1, 0), (1, 0) # directions
    turn_right = {NORTH: E, E: S, S: W, W: NORTH} # old -> new direction
    
    if width < 1 or height < 1:
        raise ValueError
    
    #x, y = int(np.ceil(width/2)-1), int(np.ceil(height/2)-1) #- 1 # start near the center

	x = 0
	y = 0
	if width%2 == 0:
		x = width // 2 - 1
	else:
		x = width // 2

	if height%2 == 0:
		y = height // 2 - 1
	else:
		y = height // 2
    
    dx, dy = NORTH # initial direction
    matrix = [[None] * width for _ in range(height)]
    count = 0
    while True:
        matrix[y][x] = count # visit
        count += 1
        # try to turn right
        new_dx, new_dy = turn_right[dx,dy]
        new_x, new_y = x + new_dx, y + new_dy
        if (0 <= new_x < width and 0 <= new_y < height and
            matrix[new_y][new_x] is None): # can turn right
            x, y = new_x, new_y
            dx, dy = new_dx, new_dy
        else: # try to move straight
            x, y = x + dx, y + dy
            if not (0 <= x < width and 0 <= y < height):
                return matrix # nowhere to go

print(make_spiral_grid(8, 8))

def make_column_grid(width, height):
    '''
    The FIJI stitching algorithm can stitch images in a column-by-column grid.
    This function makes a column-shaped grid.
    '''
    matrix = np.arange(width * height).reshape(height,width).transpose()
    return matrix
