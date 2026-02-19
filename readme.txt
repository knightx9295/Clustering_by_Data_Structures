******************** INTRO ********************
The above code is for clustering of nearest points in 2D-plane until
K-groups are formed. It uses two data structures to display the difference.
    Part-1 uses standard 2D matrix (adjacenecy lists) to store and locate the nearest points
    and merge to form groups
    Part-2 uses min heap to keep the minimum distance point near the top and
    then extract to merge.
You can run one by one each file as each uses 2 different data structures

******************** File Input ****************
The code uses file.txt as input to store info about point, no of points and
their x-y coordinates to locate in 2D-plane

******************** Constraints ***************
Use the file.txt to change the input.
First row has no of points on left and the final k groups to form on right
Then all the lines below have (x,y) coordinate pair separated by commas
Use the already present format and above constriants to check the grouping of points